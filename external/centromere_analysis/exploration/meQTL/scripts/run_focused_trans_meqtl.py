#!/usr/bin/env python3
"""Scan prioritized SVs against one off-chromosome methylation chunk with permutations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from run_sv_methylation_assoc import MATRIX_META_COLUMNS, bh_rows_for_sv, residualize


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sv-dosages", type=Path, required=True)
    parser.add_argument("--sv-metadata", type=Path, required=True)
    parser.add_argument("--methylation-matrix", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--n-permutations", type=int, default=100)
    parser.add_argument("--permutation-seed-base", type=int, default=739391)
    parser.add_argument("--permutation-batch-size", type=int, default=10)
    parser.add_argument("--pvalue-threshold", type=float, default=1e-4)
    parser.add_argument("--top-per-sv", type=int, default=20)
    parser.add_argument("--min-residual-df", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--null-stats", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def association_statistics(
    x_resid: np.ndarray, y_resid: np.ndarray, y_ss: np.ndarray, df: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_ss = np.sum(x_resid * x_resid, axis=0)
    xy = x_resid.T @ y_resid
    beta = np.divide(xy, x_ss[:, None], out=np.zeros_like(xy), where=x_ss[:, None] > 0)
    rss = np.maximum(y_ss[None, :] - beta * beta * x_ss[:, None], 0.0)
    se = np.sqrt(
        np.divide(rss / df, x_ss[:, None], out=np.full_like(rss, np.inf), where=x_ss[:, None] > 0)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = beta / se
        r_value = xy / np.sqrt(x_ss[:, None] * y_ss[None, :])
    p_values = 2.0 * stats.t.sf(np.abs(t_stat), df=df)
    p_values = np.where(np.isfinite(p_values), p_values, 1.0)
    return beta, se, t_stat, r_value, p_values


def make_permutation_seeds(base: int, count: int) -> np.ndarray:
    return np.random.default_rng(base).integers(1, np.iinfo(np.int32).max, size=count, dtype=np.int64)


def main() -> int:
    args = parse_args()
    meth = pd.read_csv(args.methylation_matrix, sep="\t")
    cov = pd.read_csv(args.covariates, sep="\t").set_index("sample_id")
    dosage = pd.read_csv(args.sv_dosages, sep="\t").set_index("sv_index", drop=False)
    metadata = pd.read_csv(args.sv_metadata, sep="\t").set_index("sv_index", drop=False)

    methyl_chroms = meth["chrom"].astype(str).dropna().unique().tolist()
    if len(methyl_chroms) != 1:
        raise ValueError(f"Expected one methylation chromosome; found {methyl_chroms}")
    methyl_chrom = methyl_chroms[0]
    active_indices = [idx for idx in dosage.index if metadata.loc[idx, "chrom"] != methyl_chrom]
    dosage = dosage.loc[active_indices]
    metadata = metadata.loc[active_indices]

    sample_cols = [column for column in meth.columns if column not in MATRIX_META_COLUMNS]
    samples = [sample for sample in sample_cols if sample in dosage.columns and sample in cov.index]
    if not samples or meth.empty or dosage.empty:
        raise ValueError("Focused trans chunk has no samples, methylation features, or off-chromosome SVs")

    y_raw = meth[samples].to_numpy(dtype=float)
    n_observed = np.sum(np.isfinite(y_raw), axis=1)
    valid = n_observed > 0
    meth = meth.loc[valid].reset_index(drop=True)
    y_raw = y_raw[valid]
    n_observed = n_observed[valid]
    row_means = np.nansum(y_raw, axis=1) / n_observed
    missing_rows, missing_cols = np.where(~np.isfinite(y_raw))
    y_raw[missing_rows, missing_cols] = row_means[missing_rows]

    y_values = y_raw.T
    x_values = dosage[samples].to_numpy(dtype=float).T
    cov_values = cov.loc[samples].to_numpy(dtype=float)
    y_resid, cov_rank = residualize(y_values, cov_values)
    x_resid, _ = residualize(x_values, cov_values)
    y_ss = np.sum(y_resid * y_resid, axis=0)
    df = len(samples) - cov_rank - 1
    if df < args.min_residual_df:
        raise ValueError(f"Only {df} residual degrees of freedom; require {args.min_residual_df}")

    beta, se, t_stat, r_value, p_values = association_statistics(x_resid, y_resid, y_ss, df)
    feature_meta = meth[[column for column in MATRIX_META_COLUMNS if column in meth.columns]].copy()
    feature_meta["n_observed"] = n_observed
    rows: list[dict] = []
    for sv_position, sv_index in enumerate(dosage.index):
        rows.extend(
            bh_rows_for_sv(
                metadata.loc[sv_index],
                feature_meta,
                beta[sv_position],
                se[sv_position],
                t_stat[sv_position],
                r_value[sv_position],
                p_values[sv_position],
                args.pvalue_threshold,
                args.top_per_sv,
                False,
            )
        )
    observed = pd.DataFrame(rows)
    if not observed.empty:
        observed.insert(0, "n_samples", len(samples))
        observed = observed.sort_values("pvalue")

    seeds = make_permutation_seeds(args.permutation_seed_base, args.n_permutations)
    null_rows: list[dict] = []
    threshold_t = float(stats.t.isf(args.pvalue_threshold / 2.0, df=df))
    batch_size = max(1, args.permutation_batch_size)
    for batch_start in range(0, len(seeds), batch_size):
        batch_seeds = seeds[batch_start : batch_start + batch_size]
        residual_columns = []
        column_keys: list[tuple[int, int]] = []
        for seed in batch_seeds:
            permutation = np.random.default_rng(int(seed)).permutation(len(samples))
            permuted_resid, _ = residualize(x_values[permutation, :], cov_values)
            residual_columns.append(permuted_resid)
            column_keys.extend((int(seed), int(sv_index)) for sv_index in dosage.index)
        stacked = np.column_stack(residual_columns)
        x_ss = np.sum(stacked * stacked, axis=0)
        xy = stacked.T @ y_resid
        with np.errstate(divide="ignore", invalid="ignore"):
            r = xy / np.sqrt(x_ss[:, None] * y_ss[None, :])
            r2 = np.minimum(r * r, np.nextafter(1.0, 0.0))
            abs_t = np.sqrt(df * r2 / np.maximum(1.0 - r2, np.finfo(float).tiny))
        abs_t = np.where(np.isfinite(abs_t), abs_t, 0.0)
        maxima = np.max(abs_t, axis=1)
        counts = np.sum(abs_t >= threshold_t, axis=1)
        min_p = 2.0 * stats.t.sf(maxima, df=df)
        for position, (seed, sv_index) in enumerate(column_keys):
            null_rows.append(
                {
                    "permutation_seed": seed,
                    "sv_index": sv_index,
                    "sv_id": metadata.loc[sv_index, "id"],
                    "methyl_chrom": methyl_chrom,
                    "n_tests": int(len(meth)),
                    "max_abs_t": float(maxima[position]),
                    "min_pvalue": float(min_p[position]),
                    "n_p_le_threshold": int(counts[position]),
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.null_stats.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    observed.to_csv(args.output, sep="\t", index=False, compression="gzip")
    pd.DataFrame(null_rows).to_csv(args.null_stats, sep="\t", index=False, compression="gzip")
    summary = {
        "methylation_matrix": str(args.methylation_matrix),
        "methyl_chrom": methyl_chrom,
        "n_samples": len(samples),
        "n_svs": len(dosage),
        "n_features": len(meth),
        "n_tests": int(len(dosage) * len(meth)),
        "n_retained": len(observed),
        "min_pvalue": float(observed["pvalue"].min()) if not observed.empty else None,
        "df": int(df),
        "n_permutations": int(args.n_permutations),
        "permutation_seed_base": int(args.permutation_seed_base),
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
