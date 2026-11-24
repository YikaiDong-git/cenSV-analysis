#!/usr/bin/env python3
"""Summarize MAC8 cis effect direction across prespecified scientific strata."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def summarize(frame: pd.DataFrame, dimension: str, value: str) -> dict:
    return {
        "dimension": dimension,
        "stratum": value,
        "n_pairs": int(len(frame)),
        "n_svs": int(frame.sv_id.nunique()),
        "n_cpgs": int(frame.feature_id.nunique()),
        "negative_raw_beta_fraction": float(frame.expected_raw_beta.lt(0).mean()),
        "median_raw_beta": float(frame.expected_raw_beta.median()),
        "negative_rint_beta_fraction": float(frame.baseline_beta.lt(0).mean()),
        "median_rint_beta": float(frame.baseline_beta.median()),
        "negative_depth_adjusted_beta_fraction": float(frame.depth_adjusted_beta.lt(0).mean()),
        "median_depth_adjusted_beta": float(frame.depth_adjusted_beta.median()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-results", type=Path, required=True)
    parser.add_argument("--censat", type=Path, required=True)
    parser.add_argument("--conditioning", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    frame = pd.read_csv(args.audit_results, sep="\t", low_memory=False)
    frame = frame.loc[frame.is_discovery.astype(bool)].copy()
    censat = pd.read_csv(args.censat, sep="\t", low_memory=False)[["feature_id", "satellite_class"]].drop_duplicates()
    frame = frame.merge(censat, on="feature_id", how="left")
    frame["distance_to_sv_bp"] = np.maximum(
        frame.target_start.to_numpy(float) - frame.source_end.to_numpy(float),
        frame.source_pos.to_numpy(float) - frame.target_end.to_numpy(float),
    ).clip(min=0)
    frame["distance_class"] = pd.cut(
        frame.distance_to_sv_bp,
        [-1, 100_000, 1_000_000, 5_000_000, 25_000_000, np.inf],
        labels=["<=100kb", "0.1-1Mb", "1-5Mb", "5-25Mb", ">25Mb"],
    ).astype(str)
    frame["sv_type"] = frame.sv_id.astype(str).str.extract(r"_(DEL|INS|DUP|INV|BND)_", expand=False).fillna("other")
    frame["coverage_robustness"] = np.where(
        frame.depth_adjusted_p.lt(0.05)
        & frame.depth_adjusted_effect_retention.abs().ge(0.5)
        & np.sign(frame.depth_adjusted_beta).eq(np.sign(frame.baseline_beta)),
        "depth_robust", "depth_sensitive",
    )
    conditioning = pd.read_csv(args.conditioning, sep="\t", low_memory=False)
    retained = set(
        zip(
            conditioning.loc[conditioning.upstream_aligned_supported_scope.fillna(False), "sv_id"].astype(str),
            conditioning.loc[conditioning.upstream_aligned_supported_scope.fillna(False), "feature_id"].astype(str),
        )
    )
    frame["primary_conditioning_status"] = [
        "primary_retained" if (str(sv), str(feature)) in retained else "not_primary_retained"
        for sv, feature in zip(frame.sv_id, frame.feature_id)
    ]
    rows = [summarize(frame, "overall", "all_cis_discoveries")]
    for dimension in ["distance_class", "satellite_class", "sv_type", "coverage_robustness", "primary_conditioning_status"]:
        for value, group in frame.groupby(dimension, observed=True, dropna=False):
            rows.append(summarize(group, dimension, str(value)))
    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    print(output.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
