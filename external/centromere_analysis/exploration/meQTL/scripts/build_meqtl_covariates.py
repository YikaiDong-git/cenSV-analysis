#!/usr/bin/env python3
"""Build methylation association covariates for overlapping samples."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


MATRIX_META_COLUMNS = {"chrom", "start", "end", "feature_id", "n_observed"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methylation-manifest", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--upstream-covariates", type=Path, default=None)
    parser.add_argument("--include-sex", action="store_true")
    parser.add_argument("--n-genotype-pcs", type=int, default=5)
    parser.add_argument("--include-platform-covariates", action="store_true")
    parser.add_argument("--methylation-pcs", type=Path, default=None)
    parser.add_argument("--methylation-pc-matrix", type=Path, default=None)
    parser.add_argument("--n-methylation-pcs", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def read_samples(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def add_upstream_covariates(covariates: pd.DataFrame, upstream_path: Path, include_sex: bool, n_genotype_pcs: int) -> list[str]:
    added = []
    if upstream_path is None or str(upstream_path) in {"", "."} or not upstream_path.exists():
        return added
    upstream = pd.read_csv(upstream_path, sep="\t", index_col=0)
    samples = covariates.index.intersection(upstream.columns)
    if include_sex and "sex" in upstream.index:
        sex = upstream.loc["sex", samples].map({"XY": 1, "male": 1, "XX": 0, "female": 0})
        covariates.loc[samples, "sex"] = pd.to_numeric(sex, errors="coerce")
        added.append("sex")
    for pc_num in range(1, n_genotype_pcs + 1):
        row_name = f"PC{pc_num}"
        if row_name in upstream.index:
            col_name = f"gPC{pc_num}"
            covariates.loc[samples, col_name] = pd.to_numeric(upstream.loc[row_name, samples], errors="coerce")
            added.append(col_name)
    return added


def add_platform_covariates(covariates: pd.DataFrame, methyl_manifest: pd.DataFrame) -> list[str]:
    added = []
    manifest = methyl_manifest.set_index("sample_id")
    platform_cols = [col for col in ["pore", "kit", "caller", "model"] if col in manifest.columns]
    platform = manifest.reindex(covariates.index)[platform_cols].fillna("unknown").astype(str)
    max_levels = max(2, min(10, len(covariates) // 10))
    platform_cols = [col for col in platform_cols if platform[col].nunique() <= max_levels]
    if not platform_cols:
        return added
    dummies = pd.get_dummies(platform[platform_cols], prefix=platform_cols, drop_first=True, dtype=float)
    for col in dummies.columns:
        if dummies[col].nunique(dropna=True) > 1:
            covariates[col] = dummies[col]
            added.append(col)
    return added


def add_methylation_pcs(covariates: pd.DataFrame, matrix_path: Path, n_pcs: int) -> tuple[list[str], list[float]]:
    if n_pcs <= 0 or matrix_path is None or not matrix_path.exists():
        return [], []
    matrix = pd.read_csv(matrix_path, sep="\t")
    sample_cols = [col for col in matrix.columns if col not in MATRIX_META_COLUMNS]
    sample_cols = [col for col in sample_cols if col in covariates.index]
    if not sample_cols or matrix.empty:
        return [], []
    x = matrix[sample_cols].T.astype(float)
    x = x.apply(lambda col: col.fillna(col.mean()), axis=0)
    x = x.loc[:, x.var(axis=0) > 0]
    max_pcs = min(n_pcs, x.shape[0] - 1, x.shape[1])
    if max_pcs <= 0:
        return [], []
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=0).replace(0, 1)
    pca = PCA(n_components=max_pcs, svd_solver="randomized", random_state=13)
    pcs = pca.fit_transform(x.to_numpy(dtype=float))
    added = []
    for i in range(max_pcs):
        col = f"methylPC{i + 1}"
        covariates.loc[x.index, col] = pcs[:, i]
        added.append(col)
    return added, pca.explained_variance_ratio_.astype(float).tolist()


def add_precomputed_methylation_pcs(covariates: pd.DataFrame, pcs_path: Path, n_pcs: int) -> list[str]:
    if n_pcs <= 0 or pcs_path is None or not pcs_path.exists():
        return []
    pcs = pd.read_csv(pcs_path, sep="\t")
    if "sample_id" not in pcs.columns:
        raise ValueError(f"{pcs_path} is missing a sample_id column")
    pcs = pcs.set_index("sample_id")
    pc_cols = [col for col in pcs.columns if col.startswith("methylPC")]
    pc_cols = pc_cols[:n_pcs]
    added = []
    samples = covariates.index.intersection(pcs.index)
    for col in pc_cols:
        covariates.loc[samples, col] = pd.to_numeric(pcs.loc[samples, col], errors="coerce")
        added.append(col)
    return added


def main() -> int:
    args = parse_args()
    samples = read_samples(args.samples)
    covariates = pd.DataFrame(index=samples)
    methyl_manifest = pd.read_csv(args.methylation_manifest, sep="\t")

    added = []
    added.extend(add_upstream_covariates(covariates, args.upstream_covariates, args.include_sex, args.n_genotype_pcs))
    if args.include_platform_covariates:
        added.extend(add_platform_covariates(covariates, methyl_manifest))
    methyl_pc_var = []
    methyl_pc_cols = add_precomputed_methylation_pcs(covariates, args.methylation_pcs, args.n_methylation_pcs)
    if not methyl_pc_cols:
        methyl_pc_cols, methyl_pc_var = add_methylation_pcs(
            covariates,
            args.methylation_pc_matrix,
            args.n_methylation_pcs,
        )
    added.extend(methyl_pc_cols)

    covariates = covariates.apply(pd.to_numeric, errors="coerce")
    covariates = covariates.loc[:, covariates.nunique(dropna=True) > 1]
    covariates = covariates.fillna(covariates.mean(axis=0))
    covariates.insert(0, "sample_id", covariates.index)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    covariates.to_csv(args.output, sep="\t", index=False)
    summary = {
        "n_samples": len(samples),
        "n_covariates": covariates.shape[1] - 1,
        "covariates": [col for col in covariates.columns if col != "sample_id"],
        "requested_covariates": added,
        "methylation_pc_variance_explained": methyl_pc_var,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
