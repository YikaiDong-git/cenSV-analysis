#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from robust_meqtl import (
    METHOD_VERSION,
    array_fingerprint,
    canonical_sample_ids,
    design_fingerprint,
    multiplier_operators,
    sample_order_fingerprint,
)
from run_sv_methylation_assoc import MATRIX_META_COLUMNS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-manifest", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--genotype-covariates", required=True)
    parser.add_argument("--methylation-covariates", required=True)
    parser.add_argument("--n-multipliers", type=int, required=True)
    parser.add_argument("--multiplier-seed-base", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def header_samples(path: Path, metadata_columns: set[str]) -> list[str]:
    columns = pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).tolist()
    return [column for column in columns if column not in metadata_columns]


def split_names(value: str) -> list[str]:
    return [name for name in value.split(",") if name]


def main() -> int:
    args = parse_args()
    manifest = pd.read_csv(args.job_manifest, sep="\t")
    covariates = pd.read_csv(args.covariates, sep="\t")
    if "sample_id" not in covariates:
        raise ValueError("Covariate table lacks sample_id")
    genotype_names = split_names(args.genotype_covariates)
    methylation_names = split_names(args.methylation_covariates)
    covariate_names = [*genotype_names, *methylation_names]
    missing_covariates = sorted(set(covariate_names) - set(covariates.columns))
    if missing_covariates:
        raise ValueError(f"Missing configured covariates: {missing_covariates}")
    if len(covariate_names) != len(set(covariate_names)):
        raise ValueError("Configured covariate names are not unique")

    matrix_paths = [
        args.matrix_dir / f"{chunk_id}.methylation.tsv.gz"
        for chunk_id in sorted(manifest["chunk_id"].astype(str).unique())
    ]
    dosage_paths = [Path(path) for path in sorted(manifest["sv_dosages"].astype(str).unique())]
    collections = [
        header_samples(path, set(MATRIX_META_COLUMNS)) for path in matrix_paths
    ]
    collections.extend(header_samples(path, {"sv_index"}) for path in dosage_paths)
    samples = canonical_sample_ids(covariates["sample_id"].astype(str).tolist(), collections)
    covariates = covariates.assign(sample_id=covariates["sample_id"].astype(str)).set_index("sample_id")
    values = covariates.loc[samples, covariate_names].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError("Canonical covariate matrix contains non-finite values")

    seeds, mammen, rademacher = multiplier_operators(
        samples, args.multiplier_seed_base, args.n_multipliers
    )
    sample_fp = sample_order_fingerprint(samples)
    mammen_fp = array_fingerprint(mammen, labels=[*samples, "mammen"])
    rademacher_fp = array_fingerprint(rademacher, labels=[*samples, "rademacher"])
    design_fp = design_fingerprint(
        samples, covariate_names, values, seeds, mammen, rademacher
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        method_version=np.asarray(METHOD_VERSION),
        sample_ids=np.asarray(samples),
        sample_order_fingerprint=np.asarray(sample_fp),
        covariate_names=np.asarray(covariate_names),
        genotype_covariate_names=np.asarray(genotype_names),
        methylation_covariate_names=np.asarray(methylation_names),
        covariates=values,
        replicate_seeds=seeds,
        mammen=mammen,
        rademacher=rademacher,
        mammen_fingerprint=np.asarray(mammen_fp),
        rademacher_fingerprint=np.asarray(rademacher_fp),
        design_fingerprint=np.asarray(design_fp),
    )
    summary = {
        "method_version": METHOD_VERSION,
        "n_samples": len(samples),
        "sample_ids": samples,
        "sample_order_fingerprint": sample_fp,
        "n_genotype_covariates": len(genotype_names),
        "genotype_covariates": genotype_names,
        "n_methylation_covariates": len(methylation_names),
        "methylation_covariates": methylation_names,
        "n_multipliers": args.n_multipliers,
        "multiplier_seed_base": args.multiplier_seed_base,
        "mammen_fingerprint": mammen_fp,
        "rademacher_fingerprint": rademacher_fp,
        "design_fingerprint": design_fp,
        "n_matrix_headers_checked": len(matrix_paths),
        "n_dosage_headers_checked": len(dosage_paths),
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
