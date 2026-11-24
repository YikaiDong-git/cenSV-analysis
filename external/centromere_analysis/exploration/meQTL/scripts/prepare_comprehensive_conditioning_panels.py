#!/usr/bin/env python3
"""Prepare merged intervals and source/availability manifests for conditioning panels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from comprehensive_conditioning import sha256_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    return parser.parse_args()


def merged(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    output: list[list[int]] = []
    for start, end in sorted(intervals):
        if not output or start > output[-1][1]:
            output.append([start, end])
        else:
            output[-1][1] = max(output[-1][1], end)
    return [(start, end) for start, end in output]


def deterministic_tsv(frame: pd.DataFrame, path: Path) -> None:
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def ont_sources(listing: Path, samples: list[str], base_url: str) -> pd.DataFrame:
    wanted = set(samples)
    vcfs: dict[str, str] = {}
    indexes: set[str] = set()
    for line in listing.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        key = fields[-1]
        if key.endswith("_PMDV_FINAL.phased.vcf.gz"):
            sample = Path(key).name.split("_PMDV_FINAL", 1)[0]
            if sample in wanted:
                vcfs[sample] = key
        elif key.endswith("_PMDV_FINAL.phased.vcf.gz.tbi"):
            indexes.add(key[:-4])
    if set(vcfs) != wanted:
        raise ValueError(f"ONT listing lacks required samples: {sorted(wanted - set(vcfs))}")
    rows = []
    for sample in samples:
        key = vcfs[sample]
        rows.append({
            "panel": "ont_source_locus", "sample_id": sample, "chrom": "all",
            "source_url": base_url.rstrip("/") + "/" + key,
            "index_url": base_url.rstrip("/") + "/" + key + ".tbi",
            "index_listed": key in indexes,
            "source_status": "remote_object_deleted; existing 1/5-Mb local panel retained",
        })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text())
    prep = config["panel_preparation"]
    root = Path(prep["output_dir"])
    regions_dir = root / "regions"
    regions_dir.mkdir(parents=True, exist_ok=True)
    svs_path = Path(config["workflow_output"]) / "manifests/svs.tsv"
    svs = pd.read_csv(svs_path, sep="\t")
    samples_path = Path(prep["samples"])
    samples = [value for value in samples_path.read_text().splitlines() if value]
    if len(samples) != 124 or len(samples) != len(set(samples)):
        raise ValueError(f"Expected 124 unique conditioning samples; observed {len(samples)}")
    max_window = int(prep["max_window_bp"])
    region_rows = []
    regions_by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, group in svs.groupby("chrom", sort=False):
        intervals = [(max(0, int(row.pos) - max_window), int(row.end) + max_window) for row in group.itertuples()]
        regions_by_chrom[str(chrom)] = merged(intervals)
    for chrom in sorted(regions_by_chrom, key=lambda value: int(value.removeprefix("chr"))):
        path = regions_dir / f"{chrom}.bed"
        text = "".join(f"{chrom}\t{start}\t{end}\n" for start, end in regions_by_chrom[chrom])
        path.write_text(text)
        for index, (start, end) in enumerate(regions_by_chrom[chrom]):
            region_rows.append({"chrom": chrom, "region_index": index, "start": start, "end": end,
                                "bases": end - start, "region_file": str(path), "region_sha256": sha256_file(path)})
    regions = pd.DataFrame(region_rows)
    deterministic_tsv(regions, root / "region_manifest.tsv")

    source_rows = []
    ont_listing = Path(prep["ont_s3_listing"])
    source_rows.extend(ont_sources(ont_listing, samples, str(prep["ont_base_url"])).to_dict("records"))
    joint_template = str(prep["joint_remote_vcf_template"])
    for chrom in regions.chrom.drop_duplicates():
        url = joint_template.format(chrom=chrom)
        source_rows.append({
            "panel": "joint_chm13_source_locus", "sample_id": "2504_source_subset_to_124", "chrom": chrom,
            "source_url": url, "index_url": url + ".tbi", "index_listed": True,
            "source_status": "remote_tabix_source_available; range extraction required",
        })
    sources = pd.DataFrame(source_rows)
    deterministic_tsv(sources, root / "source_manifest.tsv.gz")

    coverage_rows = []
    windows = [int(value) for value in config["windows_bp"]]
    for panel_name, panel in config["panels"].items():
        available_windows = {int(value) for value in panel["available_windows_bp"]}
        for row in svs.itertuples():
            region_path = regions_dir / f"{row.chrom}.bed"
            for window in windows:
                available = window in available_windows
                if panel_name == "ont_source_locus":
                    status = "existing_complete" if available else "unavailable_remote_source_deleted"
                else:
                    status = "planned_remote_range_extraction" if available else "unavailable"
                coverage_rows.append({
                    "panel": panel_name, "sv_id": row.sv_id, "chrom": row.chrom, "window_bp": window,
                    "target_start": max(0, int(row.pos) - window), "target_end": int(row.end) + window,
                    "available_for_conditioning": available, "availability_status": status,
                    "coverage_region_file": str(region_path) if panel_name != "ont_source_locus" else str(panel["coverage_template"]).format(chrom=row.chrom),
                    "prepared_region_sha256": sha256_file(region_path),
                    "unavailable_reason": "ONT PMDV source objects carry S3 delete markers; no retained full local VCF was found" if not available and panel_name == "ont_source_locus" else "",
                })
    coverage = pd.DataFrame(coverage_rows)
    deterministic_tsv(coverage, root / "interval_coverage_manifest.tsv.gz")

    ont_checksums = []
    ont_panel = config["panels"]["ont_source_locus"]
    for chrom in regions.chrom.drop_duplicates():
        vcf = Path(str(ont_panel["vcf_template"]).format(chrom=chrom))
        index = Path(str(vcf) + ".tbi")
        region = Path(str(ont_panel["coverage_template"]).format(chrom=chrom))
        if not (vcf.exists() and index.exists() and region.exists()):
            raise ValueError(f"Existing ONT 1/5-Mb panel input is missing for {chrom}")
        ont_checksums.append({
            "panel": "ont_source_locus", "chrom": chrom, "vcf": str(vcf),
            "vcf_bytes": vcf.stat().st_size, "vcf_sha256": sha256_file(vcf),
            "index": str(index), "index_bytes": index.stat().st_size, "index_sha256": sha256_file(index),
            "coverage": str(region), "coverage_sha256": sha256_file(region),
        })
    deterministic_tsv(pd.DataFrame(ont_checksums), root / "existing_ont_panel_checksums.tsv")

    baseline_root = Path(prep["storage_projection_baseline_regions"])
    baseline_vcf_root = Path(prep["storage_projection_baseline_vcfs"])
    baseline_bp = 0
    for path in baseline_root.glob("*.bed"):
        for line in path.read_text().splitlines():
            fields = line.split("\t")
            baseline_bp += int(fields[2]) - int(fields[1])
    baseline_bytes = sum(path.stat().st_size for path in baseline_vcf_root.glob("*.vcf.gz"))
    projected_bytes = int(regions.bases.sum() * baseline_bytes / baseline_bp) if baseline_bp else 0
    cap = int(prep["storage_cap_bytes"])
    if projected_bytes * float(prep["storage_safety_factor"]) >= cap:
        raise ValueError(f"Projected joint panel storage exceeds cap: {projected_bytes} bytes")
    summary = {
        "n_svs": int(svs.sv_id.nunique()), "n_chromosomes": int(regions.chrom.nunique()),
        "n_merged_regions": len(regions), "merged_region_bases": int(regions.bases.sum()),
        "n_required_samples": len(samples), "ont_sources_listed": 124,
        "ont_remote_status": "unavailable: representative exact source returned HTTP 404 with x-amz-delete-marker=true on 2026-08-26",
        "joint_remote_status": "available: chr1 header range request returned 2504 samples including all 124 required samples on 2026-08-26",
        "baseline_joint_region_bases": baseline_bp, "baseline_joint_vcf_bytes": baseline_bytes,
        "projected_joint_vcf_bytes": projected_bytes,
        "projected_with_safety_factor_bytes": int(projected_bytes * float(prep["storage_safety_factor"])),
        "storage_cap_bytes": cap, "under_storage_cap": True,
        "input_sha256": {str(path): sha256_file(path) for path in [args.config, svs_path, samples_path, ont_listing]},
        "output_sha256": {str(path): sha256_file(path) for path in [root / "region_manifest.tsv", root / "source_manifest.tsv.gz", root / "interval_coverage_manifest.tsv.gz", root / "existing_ont_panel_checksums.tsv"]},
    }
    (root / "preparation_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
