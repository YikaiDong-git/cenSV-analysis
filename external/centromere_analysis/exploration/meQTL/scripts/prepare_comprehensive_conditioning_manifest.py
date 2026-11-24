#!/usr/bin/env python3
"""Prepare the exact comprehensive PC10 cis-conditioning manifest and task graph."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from comprehensive_conditioning import (
    METHOD_VERSION, build_pair_union, expand_analysis_tasks, make_batches,
    sha256_file, stable_frame_hash,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    return parser.parse_args()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text())
    output = Path(config["workflow_output"]) / "manifests"
    output.mkdir(parents=True, exist_ok=True)
    hierarchy_path = Path(config["hierarchical_discoveries"])
    global_path = Path(config["global_pair_discoveries"])
    job_path = Path(config["job_manifest"])
    design_path = Path(config["canonical_design"])
    hierarchical = pd.read_csv(hierarchy_path, sep="\t")
    global_pairs = pd.read_csv(global_path, sep="\t")
    pairs = build_pair_union(hierarchical, global_pairs)

    expected = config["expected"]
    observed = {
        "pairs": len(pairs),
        "features": pairs.feature_id.nunique(),
        "svs": pairs.sv_id.nunique(),
        "hierarchical_primary": int(pairs.hierarchical_primary.sum()),
        "global_only_secondary": int((pairs.discovery_scope == "global_only_secondary").sum()),
    }
    if observed != {key: int(value) for key, value in expected.items()}:
        raise ValueError(f"Discovery identity mismatch: observed={observed}, expected={expected}")

    jobs = pd.read_csv(job_path, sep="\t")
    if jobs.job_id.duplicated().any():
        raise ValueError("Source job manifest has duplicate job_id values")
    source = jobs.set_index("job_id")
    unresolved = sorted(set(pairs.job_id.astype(str)) - set(source.index.astype(str)))
    if unresolved:
        raise ValueError(f"Pair job IDs absent from source manifest: {unresolved[:5]}")
    pairs["chunk_id"] = pairs.job_id.map(source.chunk_id)
    pairs["matrix_path"] = pairs.chunk_id.map(lambda value: str(Path(config["matrix_dir"]) / f"{value}.methylation.tsv.gz"))
    pairs["sv_dosages"] = pairs.job_id.map(source.sv_dosages)
    pairs["sv_metadata"] = pairs.job_id.map(source.sv_metadata)

    sv_rows = []
    for sv_id, group in pairs.groupby("sv_id", sort=False):
        first = group.iloc[0]
        metadata = pd.read_csv(first.sv_metadata, sep="\t")
        hit = metadata.loc[metadata.id.astype(str).eq(str(sv_id))]
        if len(hit) != 1:
            raise ValueError(f"Expected exactly one metadata row for {sv_id}; found {len(hit)}")
        meta = hit.iloc[0]
        sv_rows.append({
            "sv_id": str(sv_id), "chrom": str(meta.chrom), "pos": int(meta.pos),
            "end": int(meta.end), "svtype": str(meta.svtype), "svlen": int(meta.svlen),
            "sv_index": int(meta.sv_index), "sv_dosages": str(first.sv_dosages),
            "sv_metadata": str(first.sv_metadata), "n_pairs": len(group),
        })
    svs = pd.DataFrame(sv_rows)
    batches = make_batches(pairs, int(config["outcome_batch_size"]))
    panels = list(config["panels"])
    windows = [int(value) for value in config["windows_bp"]]
    snp_tasks, analysis_tasks = expand_analysis_tasks(batches, panels, windows)

    pair_path = output / "pairs.tsv.gz"
    sv_path = output / "svs.tsv"
    batch_path = output / "outcome_batches.tsv"
    snp_path = output / "snp_tasks.tsv"
    analysis_path = output / "analysis_tasks.tsv"
    write_tsv(pairs, pair_path)
    write_tsv(svs, sv_path)
    write_tsv(batches, batch_path)
    write_tsv(snp_tasks, snp_path)
    write_tsv(analysis_tasks, analysis_path)
    design = np.load(design_path, allow_pickle=True)
    pc_counts = design["pc_counts"].astype(int).tolist()
    covariate_counts = design["covariate_counts"].astype(int).tolist()
    covariate_names = design["covariate_names"].astype(str).tolist()
    if pc_counts != [5, 10, 20, 30] or covariate_counts != [11, 16, 26, 36]:
        raise ValueError(f"Canonical nested design contract changed: PCs={pc_counts}, columns={covariate_counts}")
    if covariate_names[:6] != ["sex", "gPC1", "gPC2", "gPC3", "gPC4", "gPC5"]:
        raise ValueError(f"Canonical design does not begin with sex plus GPC1-5: {covariate_names[:6]}")
    provenance = {
        "method_version": METHOD_VERSION,
        "config": str(args.config), "observed": observed,
        "n_outcome_batches": len(batches), "n_snp_tasks": len(snp_tasks),
        "n_analysis_tasks": len(analysis_tasks), "outcome_batch_size_max": int(batches.n_outcomes.max()),
        "panels_are_same_sample_sensitivity_callsets_not_replication": True,
        "pair_identity_sha256": stable_frame_hash(pairs, ["sv_id", "feature_id", "discovery_scope"]),
        "sample_order_fingerprint": str(design["sample_order_fingerprint"].item()),
        "design_fingerprint": str(design["design_fingerprint"].item()),
        "source_sha256": {str(path): sha256_file(path) for path in [hierarchy_path, global_path, job_path, design_path, args.config]},
        "implementation_sha256": {str(path): sha256_file(path) for path in [
            Path(__file__), Path("scripts/comprehensive_conditioning.py"),
            Path("scripts/prepare_comprehensive_conditioning_panels.py"),
            Path("scripts/extract_comprehensive_joint_panel.py"),
            Path("scripts/summarize_comprehensive_conditioning_panels.py"),
            Path("scripts/extract_comprehensive_conditioning_snps.py"),
            Path("scripts/run_comprehensive_conditioning_batch.py"),
            Path("scripts/run_comprehensive_susie_batch.R"),
            Path("scripts/summarize_comprehensive_conditioning.py"),
            Path("workflow/meqtl_comprehensive_conditioning/Snakefile"),
            Path("workflow/meqtl_comprehensive_conditioning/slurm_profile/config.v8+.yaml"),
            Path("scripts/submit_comprehensive_conditioning_snakemake.sh"),
        ]},
        "manifest_sha256": {str(path): sha256_file(path) for path in [pair_path, sv_path, batch_path, snp_path, analysis_path]},
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
