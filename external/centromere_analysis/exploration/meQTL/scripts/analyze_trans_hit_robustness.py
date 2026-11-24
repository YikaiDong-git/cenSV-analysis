#!/usr/bin/env python3
"""Evaluate sample and carrier influence for focused trans-meQTL candidates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--sv-dosages", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def fit_effect(y: np.ndarray, x: np.ndarray, covariates: np.ndarray) -> dict[str, float]:
    design = np.column_stack([np.ones(len(y)), covariates, x])
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coefficients
    rank = np.linalg.matrix_rank(design)
    df = len(y) - rank
    inverse = np.linalg.pinv(design.T @ design)
    variance = np.sum(residual**2) / df
    se = np.sqrt(max(0.0, inverse[-1, -1] * variance))
    t_value = coefficients[-1] / se if se > 0 else np.nan

    leverage = np.sum((design @ inverse) * design, axis=1)
    adjusted = residual / np.maximum(1.0 - leverage, np.finfo(float).eps)
    meat = design.T @ ((adjusted * adjusted)[:, None] * design)
    hc3_covariance = inverse @ meat @ inverse
    hc3_se = np.sqrt(max(0.0, hc3_covariance[-1, -1]))
    hc3_t = coefficients[-1] / hc3_se if hc3_se > 0 else np.nan
    return {
        "beta": float(coefficients[-1]),
        "se": float(se),
        "pvalue": float(2 * stats.t.sf(abs(t_value), df=df)),
        "hc3_se": float(hc3_se),
        "hc3_pvalue": float(2 * stats.t.sf(abs(hc3_t), df=df)),
        "df": int(df),
        "max_leverage": float(np.max(leverage)),
    }


def main() -> int:
    args = parse_args()
    hits = pd.read_csv(args.hits, sep="\t")
    if hits.empty:
        raise ValueError("No trans hits supplied for robustness analysis")
    dosage = pd.read_csv(args.sv_dosages, sep="\t").set_index("sv_index")
    covariates = pd.read_csv(args.covariates, sep="\t").set_index("sample_id")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fit_rows: list[dict] = []
    sample_rows: list[dict] = []
    for hit in hits.itertuples(index=False):
        matrix = pd.read_csv(args.matrix_dir / f"{hit.chunk_id}.methylation.tsv.gz", sep="\t")
        feature = matrix[matrix["feature_id"].eq(hit.feature_id)]
        if len(feature) != 1:
            raise ValueError(f"Expected one row for {hit.feature_id}; found {len(feature)}")
        sample_columns = [sample for sample in covariates.index if sample in matrix.columns and sample in dosage.columns]
        y_series = pd.to_numeric(feature.iloc[0][sample_columns], errors="coerce")
        x = dosage.loc[int(hit.sv_index), sample_columns].to_numpy(dtype=float)
        cov = covariates.loc[sample_columns].to_numpy(dtype=float)
        observed = np.isfinite(y_series.to_numpy(dtype=float))

        for analysis, mask, impute in [
            ("mean_imputed", np.ones(len(sample_columns), dtype=bool), True),
            ("complete_case", observed, False),
        ]:
            y = y_series.to_numpy(dtype=float).copy()
            if impute:
                y[~observed] = np.nanmean(y)
            result = fit_effect(y[mask], x[mask], cov[mask])
            fit_rows.append(
                {
                    "sv_id": hit.sv_id,
                    "feature_id": hit.feature_id,
                    "analysis": analysis,
                    "n_samples": int(mask.sum()),
                    "n_carriers": int(np.sum(x[mask] > 0)),
                    **result,
                }
            )

        y = y_series.to_numpy(dtype=float).copy()
        y[~observed] = np.nanmean(y)
        full = fit_effect(y, x, cov)
        for sample_index, sample in enumerate(sample_columns):
            keep = np.ones(len(sample_columns), dtype=bool)
            keep[sample_index] = False
            leave_one_out = fit_effect(y[keep], x[keep], cov[keep])
            sample_rows.append(
                {
                    "sv_id": hit.sv_id,
                    "feature_id": hit.feature_id,
                    "sample_id": sample,
                    "dosage": x[sample_index],
                    "is_carrier": bool(x[sample_index] > 0),
                    "methylation": y[sample_index],
                    "methylation_observed": bool(observed[sample_index]),
                    "loo_beta": leave_one_out["beta"],
                    "loo_pvalue": leave_one_out["pvalue"],
                    "loo_hc3_pvalue": leave_one_out["hc3_pvalue"],
                    "beta_change_fraction": (leave_one_out["beta"] - full["beta"]) / abs(full["beta"]),
                }
            )

    fits = pd.DataFrame(fit_rows)
    samples = pd.DataFrame(sample_rows)
    fits.to_csv(args.output_dir / "model_sensitivity.tsv", sep="\t", index=False)
    samples.to_csv(args.output_dir / "leave_one_out.tsv", sep="\t", index=False)
    carrier_rows = samples[samples["is_carrier"]]
    mean_imputed = fits[fits["analysis"].eq("mean_imputed")].iloc[0]
    summary = {
        "n_hits": int(hits.shape[0]),
        "beta": float(mean_imputed["beta"]),
        "pvalue": float(mean_imputed["pvalue"]),
        "hc3_pvalue": float(mean_imputed["hc3_pvalue"]),
        "n_carriers": int(mean_imputed["n_carriers"]),
        "carrier_loo_direction_concordance": float(
            np.mean(np.sign(carrier_rows["loo_beta"]) == np.sign(mean_imputed["beta"]))
        ),
        "carrier_loo_min_beta": float(carrier_rows["loo_beta"].min()),
        "carrier_loo_max_beta": float(carrier_rows["loo_beta"].max()),
        "maximum_absolute_carrier_beta_change_fraction": float(
            carrier_rows["beta_change_fraction"].abs().max()
        ),
        "minimum_carrier_loo_pvalue": float(carrier_rows["loo_pvalue"].min()),
        "maximum_carrier_loo_pvalue": float(carrier_rows["loo_pvalue"].max()),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
