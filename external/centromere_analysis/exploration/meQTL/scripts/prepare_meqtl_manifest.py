#!/usr/bin/env python3
"""Prepare sample, chunk, SV batch, and job manifests for the meQTL pilot."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from string import Formatter

import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    return parser.parse_args()


def read_sv_samples(dosage_path: Path) -> list[str]:
    with dosage_path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if not header or header[0] != "sv_index":
        raise ValueError(f"Expected first column sv_index in {dosage_path}")
    return header[1:]


def render_template(template: str, chrom: str, source: str = "") -> str:
    fields = {field_name for _, field_name, _, _ in Formatter().parse(template) if field_name}
    values = {"chrom": chrom, "source": source}
    missing = fields - values.keys()
    if missing:
        raise ValueError(f"Unsupported template fields in {template}: {sorted(missing)}")
    return template.format(**values)


def resolve_sv_inputs(config: dict) -> pd.DataFrame:
    rows = []
    sv_input_sets = config.get("sv_input_sets")
    if sv_input_sets:
        iterator = sv_input_sets.items() if isinstance(sv_input_sets, dict) else [(item["name"], item) for item in sv_input_sets]
        for source, item in iterator:
            source = str(source)
            sv_chroms = item.get("sv_chroms") or config.get("sv_chroms") or config.get("methylation_chroms")
            if not sv_chroms:
                raise ValueError("sv_chroms or methylation_chroms is required with sv_input_sets")
            for chrom in sv_chroms:
                chrom = str(chrom)
                rows.append(
                    {
                        "sv_source": source,
                        "sv_chrom": chrom,
                        "sv_dosages": render_template(str(item["sv_dosages_template"]), chrom, source),
                        "sv_metadata": render_template(str(item["sv_metadata_template"]), chrom, source),
                    }
                )
        inputs = pd.DataFrame(rows)
        missing = [
            path
            for path in inputs["sv_dosages"].tolist() + inputs["sv_metadata"].tolist()
            if not Path(path).exists()
        ]
        if missing:
            raise FileNotFoundError(f"Missing SV input files: {missing[:5]}")
        return inputs

    sv_inputs = config.get("sv_inputs")
    if sv_inputs:
        for chrom, item in sv_inputs.items():
            rows.append(
                {
                    "sv_source": str(item.get("sv_source", "")),
                    "sv_chrom": str(chrom),
                    "sv_dosages": str(item["sv_dosages"]),
                    "sv_metadata": str(item["sv_metadata"]),
                }
            )
    elif config.get("sv_dosages_template") and config.get("sv_metadata_template"):
        sv_chroms = config.get("sv_chroms") or config.get("methylation_chroms")
        if not sv_chroms:
            raise ValueError("sv_chroms or methylation_chroms is required with SV path templates")
        for chrom in sv_chroms:
            chrom = str(chrom)
            rows.append(
                {
                    "sv_chrom": chrom,
                    "sv_source": str(config.get("sv_source", "")),
                    "sv_dosages": render_template(str(config["sv_dosages_template"]), chrom, str(config.get("sv_source", ""))),
                    "sv_metadata": render_template(str(config["sv_metadata_template"]), chrom, str(config.get("sv_source", ""))),
                }
            )
    else:
        rows.append(
            {
                "sv_source": str(config.get("sv_source", "")),
                "sv_chrom": str((config.get("sv_chroms") or ["sv"])[0]),
                "sv_dosages": str(config["sv_dosages"]),
                "sv_metadata": str(config["sv_metadata"]),
            }
        )
    inputs = pd.DataFrame(rows)
    missing = [
        path
        for path in inputs["sv_dosages"].tolist() + inputs["sv_metadata"].tolist()
        if not Path(path).exists()
    ]
    if missing:
        raise FileNotFoundError(f"Missing SV input files: {missing[:5]}")
    return inputs


def make_chunks(chroms: list[str], chrom_sizes: dict[str, int], chunk_size: int) -> pd.DataFrame:
    rows = []
    for chrom in chroms:
        if chrom not in chrom_sizes:
            raise ValueError(f"Missing chrom_sizes entry for {chrom}")
        chrom_size = int(chrom_sizes[chrom])
        for start in range(0, chrom_size, chunk_size):
            end = min(start + chunk_size, chrom_size)
            rows.append(
                {
                    "chunk_id": f"{chrom}_{start}_{end}",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                }
            )
    return pd.DataFrame(rows)


def filter_sv_metadata(config: dict, sv_metadata: str | Path) -> pd.DataFrame:
    sv_meta = pd.read_csv(sv_metadata, sep="\t")
    keep = pd.Series(True, index=sv_meta.index)
    categories = config.get("sv_categories") or []
    if categories:
        keep &= sv_meta["sv_category"].isin(categories)
    if "min_cv_r2" in config:
        keep &= pd.to_numeric(sv_meta["cv_r2"], errors="coerce") >= float(config["min_cv_r2"])
    if "min_svfreq" in config:
        keep &= pd.to_numeric(sv_meta["svfreq"], errors="coerce") >= float(config["min_svfreq"])
    return sv_meta.loc[keep].sort_values(["chrom", "pos", "sv_index"]).reset_index(drop=True)


def make_sv_batches(
    sv_meta: pd.DataFrame,
    batch_size: int,
    sv_chrom: str,
    sv_dosages: str,
    sv_metadata: str,
    sv_source: str = "",
) -> pd.DataFrame:
    rows = []
    sv_indices = sv_meta["sv_index"].astype(str).tolist()
    prefix = f"{sv_source}_" if sv_source else ""
    for batch_index, start in enumerate(range(0, len(sv_indices), batch_size)):
        batch = sv_indices[start : start + batch_size]
        rows.append(
            {
                "sv_source": sv_source,
                "sv_chrom": sv_chrom,
                "sv_batch_id": f"{prefix}{sv_chrom}_svbatch{batch_index:04d}",
                "sv_batch_index": batch_index,
                "sv_indices": ",".join(batch),
                "n_svs": len(batch),
                "sv_dosages": sv_dosages,
                "sv_metadata": sv_metadata,
            }
        )
    return pd.DataFrame(rows)


def safe_name(sample_id: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in sample_id)


def split_paths(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [part for part in str(value).split(";") if part]


def add_indexed_methylation_paths(methyl_manifest: pd.DataFrame, config: dict, workflow_output: Path) -> pd.DataFrame:
    if not config.get("tabix_methylation", True):
        return methyl_manifest
    index_dir = Path(config.get("indexed_methylation_dir", workflow_output / "indexed_methylation"))
    methyl_manifest = methyl_manifest.copy()
    if "methylation_beds" not in methyl_manifest.columns:
        methyl_manifest["methylation_beds"] = methyl_manifest["methylation_bed"]
    methyl_manifest["original_methylation_beds"] = methyl_manifest["methylation_beds"]
    methyl_manifest["original_methylation_bed"] = methyl_manifest["methylation_bed"]
    local_bed_sets = []
    local_index_sets = []
    for row in methyl_manifest.itertuples(index=False):
        sample_id = safe_name(str(row.sample_id))
        sources = split_paths(row.methylation_beds)
        local_beds = [index_dir / f"{sample_id}__part{i:03d}.methylation.bed.gz" for i in range(1, len(sources) + 1)]
        local_bed_sets.append(";".join(map(str, local_beds)))
        local_index_sets.append(";".join(f"{bed}.tbi" for bed in local_beds))
    methyl_manifest["methylation_beds"] = local_bed_sets
    methyl_manifest["methylation_indexes"] = local_index_sets
    methyl_manifest["methylation_bed"] = methyl_manifest["methylation_beds"].map(lambda value: split_paths(value)[0])
    methyl_manifest["methylation_index"] = methyl_manifest["methylation_indexes"].map(lambda value: split_paths(value)[0])
    methyl_manifest["indexed_mode"] = "pending"
    return methyl_manifest


def main() -> int:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text()) or {}
    workflow_output = Path(config["workflow_output"])
    manifest_dir = workflow_output / "manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    sv_metadata_dir = manifest_dir / "filtered_sv_metadata"
    sv_metadata_dir.mkdir(parents=True, exist_ok=True)
    sv_inputs = resolve_sv_inputs(config)

    methyl_manifest_path = Path(config["methylation_manifest"])
    methyl_manifest = pd.read_csv(methyl_manifest_path, sep="\t")
    sv_sample_sets = [set(read_sv_samples(Path(path))) for path in sv_inputs["sv_dosages"]]
    sv_samples = set.intersection(*sv_sample_sets)
    methyl_manifest = methyl_manifest[methyl_manifest["sample_id"].isin(sv_samples)].copy()
    methyl_manifest = methyl_manifest.sort_values("sample_id").reset_index(drop=True)
    methyl_manifest = add_indexed_methylation_paths(methyl_manifest, config, workflow_output)
    overlap_samples = methyl_manifest["sample_id"].tolist()
    if not overlap_samples:
        raise ValueError("No overlapping samples between methylation manifest and SV dosages")

    samples_path = manifest_dir / "overlap_samples.txt"
    samples_path.write_text("\n".join(overlap_samples) + "\n")
    methyl_manifest.to_csv(manifest_dir / "overlap_methylation_samples.tsv", sep="\t", index=False)

    chunks = make_chunks(
        list(config["methylation_chroms"]),
        {str(k): int(v) for k, v in config["chrom_sizes"].items()},
        int(config["methylation_chunk_size_bp"]),
    )
    chunks.to_csv(manifest_dir / "methylation_chunks.tsv", sep="\t", index=False)

    sv_meta_frames = []
    sv_batch_frames = []
    for sv_input in sv_inputs.itertuples(index=False):
        sv_meta = filter_sv_metadata(config, sv_input.sv_metadata)
        sv_meta = sv_meta.copy()
        sv_source = getattr(sv_input, "sv_source", "")
        sv_meta["sv_source"] = sv_source
        sv_meta["sv_source_chrom"] = sv_input.sv_chrom
        filtered_prefix = f"{sv_source}_" if sv_source else ""
        filtered_path = sv_metadata_dir / f"{filtered_prefix}{sv_input.sv_chrom}.filtered_sv_metadata.tsv"
        sv_meta.to_csv(filtered_path, sep="\t", index=False)
        sv_meta_frames.append(sv_meta)
        if not sv_meta.empty:
            sv_batch_frames.append(
                make_sv_batches(
                    sv_meta,
                    int(config["sv_batch_size"]),
                    sv_input.sv_chrom,
                    sv_input.sv_dosages,
                    str(filtered_path),
                    sv_source,
                )
            )
    sv_meta_all = pd.concat(sv_meta_frames, ignore_index=True) if sv_meta_frames else pd.DataFrame()
    sv_meta_all.to_csv(manifest_dir / "filtered_sv_metadata.tsv", sep="\t", index=False)
    sv_batches = pd.concat(sv_batch_frames, ignore_index=True) if sv_batch_frames else pd.DataFrame(
        columns=["sv_source", "sv_chrom", "sv_batch_id", "sv_batch_index", "sv_indices", "n_svs", "sv_dosages", "sv_metadata"]
    )
    sv_batches.to_csv(manifest_dir / "sv_batches.tsv", sep="\t", index=False)

    job_rows = []
    job_index = 0
    same_chrom_only = bool(config.get("same_chrom_only", False))
    for chunk in chunks.itertuples(index=False):
        eligible_batches = sv_batches
        if same_chrom_only and not sv_batches.empty:
            eligible_batches = sv_batches[sv_batches["sv_chrom"].astype(str) == str(chunk.chrom)]
        for batch in eligible_batches.itertuples(index=False):
            job_rows.append(
                {
                    "job_id": f"meqtl{job_index:06d}",
                    "job_index": job_index,
                    "chunk_id": chunk.chunk_id,
                    "sv_batch_id": batch.sv_batch_id,
                    "chrom": chunk.chrom,
                    "start": chunk.start,
                    "end": chunk.end,
                    "sv_source": getattr(batch, "sv_source", ""),
                    "sv_chrom": batch.sv_chrom,
                    "sv_dosages": batch.sv_dosages,
                    "sv_metadata": batch.sv_metadata,
                    "sv_indices": batch.sv_indices,
                    "n_svs": batch.n_svs,
                }
            )
            job_index += 1
    job_manifest = pd.DataFrame(job_rows)
    job_manifest.to_csv(manifest_dir / "job_manifest.tsv", sep="\t", index=False)

    run_manifest = {
        "config": str(args.config),
        "workflow_output": str(workflow_output),
        "methylation_manifest": str(methyl_manifest_path),
        "n_overlap_samples": len(overlap_samples),
        "n_methylation_chunks": len(chunks),
        "n_sv_inputs": int(len(sv_inputs)),
        "same_chrom_only": same_chrom_only,
        "n_filtered_svs": int(len(sv_meta_all)),
        "n_sv_batches": len(sv_batches),
        "n_jobs": len(job_manifest),
        "methylation_tests_planned_upper_bound": int(len(job_manifest) * int(config["sv_batch_size"])),
    }
    (manifest_dir / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
