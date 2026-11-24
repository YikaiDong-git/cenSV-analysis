#!/usr/bin/env python3
"""Run residualized SV dosage to methylation association for one job chunk."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


MATRIX_META_COLUMNS = ["chrom", "start", "end", "feature_id", "n_observed"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sv-dosages", type=Path, required=True)
    parser.add_argument("--sv-metadata", type=Path, required=True)
    parser.add_argument("--sv-indices", required=True)
    parser.add_argument("--methylation-matrix", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--pvalue-threshold", type=float, default=1e-4)
    parser.add_argument("--top-per-sv", type=int, default=10)
    parser.add_argument("--min-residual-df", type=int, default=10)
    parser.add_argument("--write-all", action="store_true")
    parser.add_argument("--permute-genotype-seed", type=int, default=None)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def residualize(matrix: np.ndarray, covariates: np.ndarray) -> tuple[np.ndarray, int]:
    intercept = np.ones((matrix.shape[0], 1), dtype=float)
    design = np.column_stack([intercept, covariates]) if covariates.size else intercept
    u, singular_values, _ = np.linalg.svd(design, full_matrices=False)
    if singular_values.size == 0:
        return matrix, 0
    tol = np.finfo(float).eps * max(design.shape) * float(singular_values[0])
    rank = int(np.sum(singular_values > tol))
    if rank == 0:
        return matrix, 0
    basis = u[:, :rank]
    return matrix - basis @ (basis.T @ matrix), rank


def bh_rows_for_sv(
    sv_row: pd.Series,
    feature_meta: pd.DataFrame,
    beta: np.ndarray,
    se: np.ndarray,
    t_stat: np.ndarray,
    r_value: np.ndarray,
    p_values: np.ndarray,
    pvalue_threshold: float,
    top_per_sv: int,
    write_all: bool,
) -> list[dict]:
    if write_all:
        keep = np.arange(len(p_values))
        reasons = np.full(len(p_values), "all", dtype=object)
    else:
        threshold_keep = set(np.where(p_values <= pvalue_threshold)[0].tolist())
        top_keep = set()
        if top_per_sv > 0 and len(p_values) > 0:
            top_n = min(top_per_sv, len(p_values))
            top_keep = set(np.argpartition(p_values, top_n - 1)[:top_n].tolist())
        keep = np.array(sorted(threshold_keep | top_keep), dtype=int)
        reasons = []
        for idx in keep:
            reason = []
            if idx in threshold_keep:
                reason.append("threshold")
            if idx in top_keep:
                reason.append("top")
            reasons.append(",".join(reason))

    rows = []
    for out_i, feature_idx in enumerate(keep):
        feature = feature_meta.iloc[feature_idx]
        rows.append(
            {
                "retained_reason": reasons[out_i] if not write_all else "all",
                "sv_index": sv_row["sv_index"],
                "sv_id": sv_row.get("id", ""),
                "sv_chrom": sv_row.get("chrom", ""),
                "sv_pos": sv_row.get("pos", np.nan),
                "sv_end": sv_row.get("end", np.nan),
                "svtype": sv_row.get("svtype", ""),
                "svlen": sv_row.get("svlen", np.nan),
                "sv_category": sv_row.get("sv_category", ""),
                "cv_r2": sv_row.get("cv_r2", np.nan),
                "svfreq": sv_row.get("svfreq", np.nan),
                "methyl_chrom": feature["chrom"],
                "methyl_start": feature["start"],
                "methyl_end": feature["end"],
                "feature_id": feature["feature_id"],
                "feature_n_observed": feature.get("n_observed", np.nan),
                "beta": beta[feature_idx],
                "se": se[feature_idx],
                "t_stat": t_stat[feature_idx],
                "r": r_value[feature_idx],
                "pvalue": p_values[feature_idx],
            }
        )
    return rows


def main() -> int:
    args = parse_args()
    sv_indices = [int(x) for x in args.sv_indices.split(",") if x]

    meth = pd.read_csv(args.methylation_matrix, sep="\t")
    cov = pd.read_csv(args.covariates, sep="\t").set_index("sample_id")
    dosage = pd.read_csv(args.sv_dosages, sep="\t")
    sv_meta = pd.read_csv(args.sv_metadata, sep="\t").set_index("sv_index", drop=False)

    dosage = dosage[dosage["sv_index"].isin(sv_indices)].set_index("sv_index", drop=False)
    dosage = dosage.loc[[idx for idx in sv_indices if idx in dosage.index]]
    sv_meta = sv_meta.loc[[idx for idx in dosage.index if idx in sv_meta.index]]

    sample_cols = [col for col in meth.columns if col not in MATRIX_META_COLUMNS]
    samples = [sample for sample in sample_cols if sample in dosage.columns and sample in cov.index]
    if not samples or meth.empty or dosage.empty:
        empty = pd.DataFrame()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        empty.to_csv(args.output, sep="\t", index=False, compression="gzip")
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(json.dumps({"n_tests": 0, "n_retained": 0}, indent=2) + "\n")
        return 0

    y = meth[samples].astype(float)
    n_observed = y.notna().sum(axis=1).to_numpy()
    y = y.apply(lambda col: col.fillna(col.mean()), axis=1)
    y_values = y.to_numpy(dtype=float).T
    x_values = dosage[samples].astype(float).to_numpy(dtype=float).T
    cov_values = cov.loc[samples].to_numpy(dtype=float)

    if args.permute_genotype_seed is not None:
        x_values = x_values[np.random.default_rng(args.permute_genotype_seed).permutation(len(samples)), :]

    x_resid, cov_rank = residualize(x_values, cov_values)
    y_resid, _ = residualize(y_values, cov_values)
    x_ss = np.sum(x_resid * x_resid, axis=0)
    y_ss = np.sum(y_resid * y_resid, axis=0)
    xy = x_resid.T @ y_resid
    df = len(samples) - cov_rank - 1
    if df < args.min_residual_df:
        raise ValueError(
            f"Covariate design leaves only {df} residual degrees of freedom "
            f"for {len(samples)} samples (rank={cov_rank}); require >= {args.min_residual_df}"
        )

    feature_meta = meth[[col for col in MATRIX_META_COLUMNS if col in meth.columns]].copy()
    rows = []
    for sv_pos, sv_index in enumerate(dosage.index):
        if x_ss[sv_pos] <= 0:
            continue
        beta = xy[sv_pos, :] / x_ss[sv_pos]
        rss = np.maximum(y_ss - beta * beta * x_ss[sv_pos], 0)
        se = np.sqrt((rss / df) / x_ss[sv_pos])
        with np.errstate(divide="ignore", invalid="ignore"):
            t_stat = beta / se
            r_value = xy[sv_pos, :] / np.sqrt(x_ss[sv_pos] * y_ss)
        p_values = 2.0 * stats.t.sf(np.abs(t_stat), df=df)
        p_values = np.where(np.isfinite(p_values), p_values, 1.0)
        sv_row = sv_meta.loc[sv_index] if sv_index in sv_meta.index else dosage.loc[sv_index]
        rows.extend(
            bh_rows_for_sv(
                sv_row,
                feature_meta,
                beta,
                se,
                t_stat,
                r_value,
                p_values,
                args.pvalue_threshold,
                args.top_per_sv,
                args.write_all,
            )
        )

    result = pd.DataFrame(rows)
    if not result.empty:
        result.insert(0, "n_samples", len(samples))
        result["feature_n_observed"] = result["feature_n_observed"].fillna(pd.Series(n_observed).median())
        result = result.sort_values("pvalue").reset_index(drop=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False, compression="gzip")

    summary = {
        "methylation_matrix": str(args.methylation_matrix),
        "n_samples": len(samples),
        "n_svs": int(len(dosage)),
        "n_features": int(len(meth)),
        "n_tests": int(len(dosage) * len(meth)),
        "n_retained": int(len(result)),
        "n_threshold": int(result["retained_reason"].str.contains("threshold").sum()) if not result.empty else 0,
        "min_pvalue": float(result["pvalue"].min()) if not result.empty else None,
        "df": int(df),
        "permute_genotype_seed": args.permute_genotype_seed,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
