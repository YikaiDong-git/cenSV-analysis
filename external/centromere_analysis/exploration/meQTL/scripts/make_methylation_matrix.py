#!/usr/bin/env python3
"""Create a sample-by-CpG methylation matrix for one genomic region."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


BED_COVERAGE_COL = 9
BED_METHYLATION_COL = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--coverage-min", type=float, default=10)
    parser.add_argument("--min-sample-fraction", type=float, default=0.75)
    parser.add_argument("--max-sites", type=int, default=0)
    parser.add_argument(
        "--no-assume-sorted",
        action="store_true",
        help="Scan entire BED files instead of stopping after the requested sorted interval.",
    )
    parser.add_argument(
        "--backend",
        choices=["auto", "duckhts", "tabix", "python"],
        default="auto",
        help="Region extraction backend. auto prefers DuckHTS for indexed files, then tabix, then gzip streaming.",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def read_samples(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def parse_methylation_fraction(raw_value: str) -> float:
    value = float(raw_value)
    if value > 1.0:
        return value / 100.0
    return value


def find_tabix_index(path: Path, index_path: Path | None = None) -> Path | None:
    if index_path and index_path.exists():
        return index_path
    for candidate in [Path(f"{path}.tbi"), Path(f"{path}.csi")]:
        if candidate.exists():
            return candidate
    return None


def load_duckhts_connection():
    try:
        import duckdb
    except ImportError:
        return None
    con = duckdb.connect()
    try:
        con.execute("LOAD duckhts")
    except Exception:
        con.execute("INSTALL duckhts FROM community")
        con.execute("LOAD duckhts")
    return con


def iter_duckhts_rows(
    con,
    path: Path,
    index_path: Path,
    chrom: str,
    start: int,
    end: int,
    coverage_min: float,
):
    region = f"{chrom}:{start + 1}-{end}"
    query = """
        SELECT
            column0 AS chrom,
            CAST(column1 AS BIGINT) AS start,
            CAST(column2 AS BIGINT) AS "end",
            CAST(column10 AS DOUBLE) AS methylation,
            CAST(column9 AS DOUBLE) AS coverage
        FROM read_tabix(?, region := ?, index_path := ?, header := false, auto_detect := false)
    """
    for row in con.execute(query, [str(path), region, str(index_path)]).fetchall():
        yield (row[0], int(row[1]), int(row[2]), parse_methylation_fraction(str(row[3])), float(row[4]))


def iter_region_lines(
    path: Path,
    index_path: Path | None,
    chrom: str,
    start: int,
    end: int,
    assume_sorted: bool,
    use_tabix: bool,
) -> Iterable[str]:
    tabix = shutil.which("tabix")
    index = find_tabix_index(path, index_path)
    if use_tabix and tabix and index:
        region = f"{chrom}:{start + 1}-{end}"
        proc = subprocess.run([tabix, str(path), region], text=True, capture_output=True, check=False)
        if proc.returncode not in {0, 1}:
            raise RuntimeError(proc.stderr.strip() or f"tabix failed for {path} {region}")
        yield from proc.stdout.splitlines()
        return

    with gzip.open(path, "rt") as handle:
        seen_target_chrom = False
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= BED_METHYLATION_COL:
                continue
            if fields[0] != chrom:
                if assume_sorted and seen_target_chrom:
                    break
                continue
            seen_target_chrom = True
            row_start = int(fields[1])
            if row_start < start:
                continue
            if row_start >= end:
                if assume_sorted:
                    break
                continue
            yield line.rstrip("\n")


def iter_modkit_rows(
    path: Path,
    index_path: Path | None,
    chrom: str,
    start: int,
    end: int,
    coverage_min: float,
    assume_sorted: bool,
    use_tabix: bool,
):
    for line in iter_region_lines(path, index_path, chrom, start, end, assume_sorted, use_tabix):
        fields = line.split("\t")
        if len(fields) <= BED_METHYLATION_COL:
            continue
        coverage = float(fields[BED_COVERAGE_COL])
        methylation = parse_methylation_fraction(fields[BED_METHYLATION_COL])
        yield (fields[0], int(fields[1]), int(fields[2]), methylation, coverage)


def manifest_paths(manifest: pd.DataFrame, sample: str) -> list[tuple[Path, Path | None]]:
    row = manifest.loc[sample]
    bed_value = row["methylation_beds"] if "methylation_beds" in row else row["methylation_bed"]
    beds = [Path(value) for value in str(bed_value).split(";") if value]
    index_value = row["methylation_indexes"] if "methylation_indexes" in row else row.get("methylation_index", "")
    indexes = [Path(value) for value in str(index_value).split(";") if value and value != "nan"]
    indexes.extend([None] * (len(beds) - len(indexes)))
    return list(zip(beds, indexes, strict=True))


def deterministic_downsample(df: pd.DataFrame, max_sites: int) -> pd.DataFrame:
    if max_sites <= 0 or len(df) <= max_sites:
        return df
    keep = np.linspace(0, len(df) - 1, num=max_sites, dtype=int)
    return df.iloc[np.unique(keep)].copy()


def main() -> int:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t").set_index("sample_id")
    samples = read_samples(args.samples)
    missing = [sample for sample in samples if sample not in manifest.index]
    if missing:
        raise ValueError(f"{len(missing)} samples are missing from methylation manifest")

    values: dict[tuple[str, int, int], dict[str, float]] = defaultdict(dict)
    sample_counts = {}
    duckhts_con = None
    backend_counts = {"duckhts": 0, "tabix_or_python": 0}
    if args.backend in {"auto", "duckhts"}:
        try:
            duckhts_con = load_duckhts_connection()
        except Exception as exc:
            if args.backend == "duckhts":
                raise RuntimeError("DuckDB/DuckHTS backend requested, but duckhts could not be loaded") from exc
            duckhts_con = None
        if args.backend == "duckhts" and duckhts_con is None:
            raise RuntimeError("DuckDB/DuckHTS backend requested, but duckdb is not installed")

    for sample in samples:
        weighted: dict[tuple[str, int, int], list[float]] = defaultdict(lambda: [0.0, 0.0])
        for bed_path, requested_index in manifest_paths(manifest, sample):
            index_path = find_tabix_index(bed_path, requested_index)
            if args.backend == "duckhts" and index_path is None:
                raise RuntimeError(f"DuckHTS backend requested, but no tabix index found for {bed_path}")
            use_duckhts = args.backend in {"auto", "duckhts"} and duckhts_con is not None and index_path is not None
            if use_duckhts:
                row_iter = iter_duckhts_rows(duckhts_con, bed_path, index_path, args.chrom, args.start, args.end, args.coverage_min)
                backend_counts["duckhts"] += 1
            else:
                if args.backend == "duckhts":
                    raise RuntimeError(f"DuckHTS backend requested, but no tabix index found for {bed_path}")
                row_iter = iter_modkit_rows(
                    bed_path, index_path, args.chrom, args.start, args.end, args.coverage_min,
                    assume_sorted=not args.no_assume_sorted, use_tabix=args.backend != "python",
                )
                backend_counts["tabix_or_python"] += 1
            for chrom, row_start, row_end, methylation, coverage in row_iter:
                key = (chrom, row_start, row_end)
                weighted[key][0] += methylation * coverage
                weighted[key][1] += coverage
        n_rows = 0
        for key, (methylated_weight, coverage) in weighted.items():
            if coverage < args.coverage_min:
                continue
            values[key][sample] = methylated_weight / coverage
            n_rows += 1
        sample_counts[sample] = n_rows

    min_samples = max(1, math.ceil(args.min_sample_fraction * len(samples)))
    rows = []
    for (chrom, row_start, row_end), sample_values in values.items():
        if len(sample_values) < min_samples:
            continue
        row = {
            "chrom": chrom,
            "start": row_start,
            "end": row_end,
            "feature_id": f"{chrom}:{row_start}-{row_end}",
            "n_observed": len(sample_values),
        }
        row.update(sample_values)
        rows.append(row)

    matrix = pd.DataFrame(rows)
    if not matrix.empty:
        matrix = matrix.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
        methyl_values = matrix[samples]
        variance = methyl_values.var(axis=1, skipna=True)
        matrix = matrix.loc[variance > 0].reset_index(drop=True)
        matrix = deterministic_downsample(matrix, args.max_sites)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output, sep="\t", index=False, compression="gzip")

    summary = {
        "chrom": args.chrom,
        "start": args.start,
        "end": args.end,
        "n_samples": len(samples),
        "coverage_min": args.coverage_min,
        "min_sample_fraction": args.min_sample_fraction,
        "min_samples": min_samples,
        "assume_sorted": not args.no_assume_sorted,
        "uses_index_when_available": True,
        "backend_requested": args.backend,
        "backend_counts": backend_counts,
        "n_raw_sites": len(values),
        "n_output_sites": int(len(matrix)),
        "sample_rows_passing_coverage": sample_counts,
        "output": str(args.output),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
