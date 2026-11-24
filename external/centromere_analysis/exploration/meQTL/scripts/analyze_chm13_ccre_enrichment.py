#!/usr/bin/env python3
"""Annotate CHM13 meQTL CpGs and test matched cCRE enrichment."""

from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact


LAYERS = ("SCREEN_registry", "GM12878_active", "hematopoietic_active_any")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    worker = subparsers.add_parser("worker")
    worker.add_argument("--chrom", required=True)
    worker.add_argument("--associations", type=Path, required=True)
    worker.add_argument("--matrix-dir", type=Path, required=True)
    worker.add_argument("--reference", type=Path, required=True)
    worker.add_argument("--screen", type=Path, required=True)
    worker.add_argument("--gm12878", type=Path, required=True)
    worker.add_argument("--activity-dir", type=Path, required=True)
    worker.add_argument("--q-max", type=float, default=0.01)
    worker.add_argument("--n-permutations", type=int, default=10_000)
    worker.add_argument("--output-dir", type=Path, required=True)
    summary = subparsers.add_parser("summary")
    summary.add_argument("--chroms", nargs="+", required=True)
    summary.add_argument("--n-permutations", type=int, default=10_000)
    summary.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def matrix_start(path: Path) -> int:
    return int(path.name.split("_", 2)[1])


def load_intervals(path: Path, chrom: str) -> list[tuple[int, int, str]]:
    intervals = []
    with pysam.TabixFile(str(path)) as tabix:
        try:
            rows = tabix.fetch(chrom)
        except ValueError:
            return intervals
        for line in rows:
            fields = line.split("\t")
            metadata = fields[3] if len(fields) == 4 else fields[5] if len(fields) > 5 else ""
            intervals.append((int(fields[1]), int(fields[2]), metadata))
    return intervals


def merge_intervals(intervals: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    merged: list[list[object]] = []
    for start, end, metadata in sorted(intervals):
        if not merged or start > int(merged[-1][1]):
            merged.append([start, end, {metadata}])
        else:
            merged[-1][1] = max(int(merged[-1][1]), end)
            merged[-1][2].add(metadata)
    return [(int(start), int(end), ";".join(sorted(values))) for start, end, values in merged]


def overlap_metadata(
    intervals: list[tuple[int, int, str]], pointer: int, start: int, end: int
) -> tuple[int, list[str]]:
    while pointer < len(intervals) and intervals[pointer][1] <= start:
        pointer += 1
    metadata = []
    index = pointer
    while index < len(intervals) and intervals[index][0] < end:
        if intervals[index][1] > start:
            metadata.append(intervals[index][2])
        index += 1
    return pointer, metadata


def coverage_bin(value: int) -> str:
    if value < 100:
        return "lt100"
    if value < 110:
        return "100_109"
    if value < 120:
        return "110_119"
    return "ge120"


def sequence_bins(reference: pysam.FastaFile, chrom: str, start: int) -> tuple[int, int]:
    window_start = (start // 10_000) * 10_000
    sequence = reference.fetch(chrom, window_start, min(window_start + 10_000, reference.get_reference_length(chrom))).upper()
    canonical = sum(sequence.count(base) for base in "ACGT")
    gc = (sequence.count("G") + sequence.count("C")) / canonical if canonical else 0.0
    cpgs = sequence.count("CG")
    gc_bin = min(9, int(gc * 10))
    cpg_bin = int(np.digitize(cpgs, [10, 25, 50, 100, 200, 400]))
    return gc_bin, cpg_bin


def circular_shift_distribution(signal: np.ndarray, annotation: np.ndarray) -> np.ndarray:
    if not len(signal):
        return np.zeros(0, dtype=np.int32)
    values = np.fft.ifft(
        np.fft.fft(signal.astype(float)) * np.conj(np.fft.fft(annotation.astype(float)))
    ).real
    return np.rint(values).astype(np.int32)


def run_worker(args: argparse.Namespace) -> None:
    assoc = pd.read_csv(args.associations, sep="\t", compression="gzip")
    significant_table = assoc[
        assoc["methyl_chrom"].eq(args.chrom) & assoc["qvalue_retained_bh"].lt(args.q_max)
    ]
    significant = set(significant_table["feature_id"].astype(str))

    screen = load_intervals(args.screen, args.chrom)
    gm12878 = load_intervals(args.gm12878, args.chrom)
    activity = []
    for path in sorted(args.activity_dir.glob("CHM13v2_*_active_cCREs.bed.gz")):
        activity.extend(load_intervals(path, args.chrom))
    layers = {
        "SCREEN_registry": screen,
        "GM12878_active": gm12878,
        "hematopoietic_active_any": merge_intervals(activity),
    }
    pointers = {layer: 0 for layer in LAYERS}
    patterns: dict[tuple[str, int, str, int, int], tuple[bytearray, bytearray]] = defaultdict(
        lambda: (bytearray(), bytearray())
    )
    totals = {layer: [0, 0, 0, 0] for layer in LAYERS}
    annotations = []
    seen = set()
    sequence_cache: dict[int, tuple[int, int]] = {}
    reference = pysam.FastaFile(str(args.reference))
    try:
        matrices = sorted(args.matrix_dir.glob(f"{args.chrom}_*.methylation.tsv.gz"), key=matrix_start)
        if not matrices:
            raise FileNotFoundError(f"No matrices found for {args.chrom}")
        for matrix in matrices:
            with gzip.open(matrix, "rt") as handle:
                header = next(handle).rstrip("\n").split("\t", 5)
                if header[:5] != ["chrom", "start", "end", "feature_id", "n_observed"]:
                    raise ValueError(f"Unexpected matrix header in {matrix}")
                for line in handle:
                    fields = line.split("\t", 5)
                    start, end = int(fields[1]), int(fields[2])
                    feature_id, n_observed = fields[3], int(fields[4])
                    is_significant = feature_id in significant
                    if is_significant:
                        seen.add(feature_id)
                    sequence_window = (start // 10_000) * 10_000
                    if sequence_window not in sequence_cache:
                        sequence_cache[sequence_window] = sequence_bins(reference, args.chrom, start)
                    gc_bin, cpg_bin = sequence_cache[sequence_window]
                    block = (start // 1_000_000) * 1_000_000
                    cov_bin = coverage_bin(n_observed)
                    feature_annotation = {
                        "chrom": args.chrom,
                        "start": start,
                        "end": end,
                        "feature_id": feature_id,
                    }
                    for layer in LAYERS:
                        pointer, metadata = overlap_metadata(
                            layers[layer], pointers[layer], start, end
                        )
                        pointers[layer] = pointer
                        overlaps = bool(metadata)
                        values = totals[layer]
                        values[0] += 1
                        values[1] += int(is_significant)
                        values[2] += int(overlaps)
                        values[3] += int(is_significant and overlaps)
                        signal, annotation = patterns[(layer, block, cov_bin, gc_bin, cpg_bin)]
                        signal.append(int(is_significant))
                        annotation.append(int(overlaps))
                        if is_significant:
                            feature_annotation[f"{layer}_overlap"] = int(overlaps)
                            feature_annotation[f"{layer}_metadata"] = ";".join(metadata)
                    if is_significant:
                        annotations.append(feature_annotation)
    finally:
        reference.close()
    if seen != significant:
        raise ValueError(f"{args.chrom}: {len(significant - seen)} significant CpGs absent from matrices")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(annotations).to_csv(
        args.output_dir / f"{args.chrom}.significant_cpg_annotations.tsv", sep="\t", index=False
    )
    stats_rows = []
    simulations = {}
    for layer_index, layer in enumerate(LAYERS):
        rng = np.random.default_rng(20260713 + int(args.chrom[3:]) * 10 + layer_index)
        simulated = np.zeros(args.n_permutations, dtype=np.int32)
        expected = 0.0
        informative = 0
        strata = 0
        for key, (signal_bytes, annotation_bytes) in patterns.items():
            if key[0] != layer:
                continue
            strata += 1
            signal = np.frombuffer(signal_bytes, dtype=np.uint8)
            annotation = np.frombuffer(annotation_bytes, dtype=np.uint8)
            if not signal.any() or not annotation.any():
                continue
            shifts = circular_shift_distribution(signal, annotation)
            expected += float(shifts.mean())
            simulated += shifts[rng.integers(0, len(shifts), args.n_permutations)]
            informative += 1
        simulations[layer] = simulated
        total, n_sig, n_annotated, observed = totals[layer]
        stats_rows.append(
            {
                "chrom": args.chrom,
                "annotation_layer": layer,
                "n_tested_cpgs": total,
                "n_significant_cpgs": n_sig,
                "n_annotation_cpgs": n_annotated,
                "observed_overlap": observed,
                "matched_expected_overlap": expected,
                "n_strata": strata,
                "n_informative_strata": informative,
            }
        )
    pd.DataFrame(stats_rows).to_csv(
        args.output_dir / f"{args.chrom}.stats.tsv", sep="\t", index=False
    )
    np.savez_compressed(args.output_dir / f"{args.chrom}.null.npz", **simulations)


def run_summary(args: argparse.Namespace) -> None:
    stats = pd.concat(
        [pd.read_csv(args.output_dir / f"{chrom}.stats.tsv", sep="\t") for chrom in args.chroms],
        ignore_index=True,
    )
    annotations = pd.concat(
        [
            pd.read_csv(args.output_dir / f"{chrom}.significant_cpg_annotations.tsv", sep="\t")
            for chrom in args.chroms
            if (args.output_dir / f"{chrom}.significant_cpg_annotations.tsv").stat().st_size > 1
        ],
        ignore_index=True,
    )
    results = []
    for layer in LAYERS:
        rows = stats[stats["annotation_layer"].eq(layer)]
        total = int(rows["n_tested_cpgs"].sum())
        significant = int(rows["n_significant_cpgs"].sum())
        annotated = int(rows["n_annotation_cpgs"].sum())
        observed = int(rows["observed_overlap"].sum())
        expected = float(rows["matched_expected_overlap"].sum())
        simulation = np.zeros(args.n_permutations, dtype=np.int64)
        for chrom in args.chroms:
            with np.load(args.output_dir / f"{chrom}.null.npz") as values:
                simulation += values[layer]
        contingency = [
            [observed, significant - observed],
            [annotated - observed, total - significant - annotated + observed],
        ]
        odds_ratio, fisher_p = fisher_exact(contingency, alternative="greater")
        matched_p = (1 + np.count_nonzero(simulation >= observed)) / (args.n_permutations + 1)
        results.append(
            {
                "annotation_layer": layer,
                "n_tested_cpgs": total,
                "n_significant_cpgs": significant,
                "n_annotation_cpgs": annotated,
                "n_significant_annotation_cpgs": observed,
                "background_fraction": annotated / total,
                "significant_fraction": observed / significant if significant else math.nan,
                "naive_fisher_odds_ratio": odds_ratio,
                "naive_fisher_p_enrichment": fisher_p,
                "matched_expected_overlap": expected,
                "matched_fold_enrichment": observed / expected if expected else math.nan,
                "matched_p_enrichment": matched_p,
                "matching": "chromosome+1Mb_block+n_observed_bin+10kb_GC_bin+10kb_reference_CpG_density_bin+circular_shift",
                "n_permutations": args.n_permutations,
            }
        )
    pd.DataFrame(results).to_csv(args.output_dir / "chm13_pc30_ccre_enrichment.tsv", sep="\t", index=False)
    annotations.to_csv(args.output_dir / "chm13_pc30_q001_cpg_annotations.tsv", sep="\t", index=False)
    (args.output_dir / "summary.json").write_text(
        json.dumps({"results": results, "n_annotated_significant_cpgs": len(annotations)}, indent=2) + "\n"
    )


def main() -> int:
    args = parse_args()
    if args.command == "worker":
        run_worker(args)
    else:
        run_summary(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
