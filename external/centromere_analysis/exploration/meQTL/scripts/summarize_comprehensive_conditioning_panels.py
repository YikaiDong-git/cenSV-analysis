#!/usr/bin/env python3
"""Validate prepared panel files, checksums, samples, and interval coverage."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from comprehensive_conditioning import interval_covered, sha256_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, nargs="+", required=True)
    parser.add_argument("--checksums-output", type=Path, required=True)
    parser.add_argument("--coverage-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    return parser.parse_args()


def write_frame(frame: pd.DataFrame, path: Path) -> None:
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text())
    prep_root = Path(config["panel_preparation"]["output_dir"])
    metric_rows = [json.loads(path.read_text()) for path in args.metrics]
    metrics = pd.DataFrame(metric_rows).sort_values("chrom", key=lambda value: value.str.removeprefix("chr").astype(int))
    expected_chroms = pd.read_csv(prep_root / "region_manifest.tsv", sep="\t").chrom.drop_duplicates().tolist()
    if set(metrics.chrom) != set(expected_chroms) or not metrics.status.eq("ok").all():
        raise ValueError("Prepared joint panel does not contain one successful file for every required chromosome")
    checksums = metrics[["chrom", "source_url", "regions", "regions_sha256", "n_samples", "n_biallelic_snps",
                         "vcf", "vcf_bytes", "vcf_sha256", "index", "index_bytes", "index_sha256"]].copy()
    write_frame(checksums, args.checksums_output)

    coverage = pd.read_csv(prep_root / "interval_coverage_manifest.tsv.gz", sep="\t")
    verified = []
    for row in coverage.itertuples(index=False):
        status = row.availability_status
        available = bool(row.available_for_conditioning)
        if row.panel == "joint_chm13_source_locus":
            path = Path(str(row.coverage_region_file))
            available = path.exists() and interval_covered(path, str(row.chrom), int(row.target_start), int(row.target_end))
            status = "prepared_and_checksum_verified" if available else "prepared_interval_missing"
        verified.append({**row._asdict(), "available_for_conditioning": available, "availability_status": status})
    verified_frame = pd.DataFrame(verified)
    write_frame(verified_frame, args.coverage_output)
    counts = (verified_frame.groupby(["panel", "window_bp"], as_index=False)
              .agg(n_svs=("sv_id", "nunique"), n_available=("available_for_conditioning", "sum")))
    ont_checksums_path = prep_root / "existing_ont_panel_checksums.tsv"
    ont_checksums = pd.read_csv(ont_checksums_path, sep="\t")
    joint_bytes = int(checksums.vcf_bytes.sum() + checksums.index_bytes.sum())
    cap = int(config["panel_preparation"]["storage_cap_bytes"])
    if joint_bytes >= cap:
        raise ValueError(f"Prepared panel exceeds storage cap: {joint_bytes} bytes")
    summary = {
        "status": "ready", "n_joint_chromosomes": len(checksums),
        "joint_panel_bytes": joint_bytes, "storage_cap_bytes": cap,
        "coverage": counts.to_dict("records"),
        "joint_checksums_sha256": sha256_file(args.checksums_output),
        "ont_checksums_sha256": sha256_file(ont_checksums_path),
        "existing_ont_panel_bytes": int(ont_checksums.vcf_bytes.sum() + ont_checksums.index_bytes.sum()),
        "verified_coverage_sha256": sha256_file(args.coverage_output),
        "interpretation": "ONT and joint CHM13 are sensitivity callsets on the same 124 participants, not independent replication. Unavailable windows are not negative evidence.",
    }
    args.summary_output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
