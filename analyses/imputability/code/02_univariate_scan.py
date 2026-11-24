#!/usr/bin/env python3
"""Univariate predictor scan: for each feature, Spearman + AUC vs imputability."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"

INTRINSIC = [
    "svlen", "log_svlen", "maf", "log_maf",
    "n_carriers", "n_hom_alt",
    "svtype_INS", "svtype_DEL", "svtype_INV", "svtype_DUP",
    "cat_monophyletic", "cat_residual", "cat_non_decomposed", "cat_original_decomposed",
]
LOCUS = [
    "n_window_snps", "log_n_window_snps",
    "best_snp_dist_mean", "best_snp_dist_min", "best_snp_dist_max",
]
CHROM = [
    "active_HOR_span_bp", "log_active_HOR_span", "n_HOR_intervals",
    "pericen_len_100000", "pericen_len_200000", "pericen_len_500000",
    "pericen_len_1000000", "pericen_len_2000000",
    "n_train_haps", "acrocentric",
]
DIAGNOSTIC = [
    "best_l1_ratio", "best_alpha", "best_n_nonzero",
    "cv_r2_var", "cv_r2_range",
]


def scan(df, feat, y_col, kind):
    x = df[feat].astype(float).values
    y = df[y_col].astype(float).values
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 30:
        return None
    x, y = x[mask], y[mask]
    rho, p_rho = spearmanr(x, y)
    out = {"feature": feat, "kind": kind, "y": y_col,
           "n": int(mask.sum()),
           "spearman_rho": float(rho), "spearman_p": float(p_rho)}
    if set(np.unique(y).astype(int)).issubset({0, 1}) and len(np.unique(y)) == 2:
        try:
            auc = roc_auc_score(y, x)
            out["auc"] = float(auc)
            out["auc_distance"] = abs(auc - 0.5)
        except Exception:
            out["auc"] = np.nan
            out["auc_distance"] = np.nan
    return out


def main():
    feat_path = OUT / "per_sv_features.tsv"
    if not feat_path.exists():
        print(f"[FAIL] {feat_path} missing — run 01 first", flush=True)
        sys.exit(1)
    df = pd.read_csv(feat_path, sep="\t", low_memory=False)
    print(f"[INFO] loaded {len(df)} SVs × {df.shape[1]} cols", flush=True)

    rows = []
    feat_groups = [("intrinsic", INTRINSIC), ("locus", LOCUS),
                   ("chromosome", CHROM), ("diagnostic", DIAGNOSTIC)]
    for y_col in ["y_cv_r2", "y_pass_01", "y_pass_03"]:
        for kind, feats in feat_groups:
            for f in feats:
                if f not in df.columns:
                    continue
                rec = scan(df, f, y_col, kind)
                if rec is not None:
                    rows.append(rec)

    out = pd.DataFrame(rows).sort_values(
        ["y", "kind", "spearman_p"]).reset_index(drop=True)
    out_path = OUT / "univariate_ranking.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {out_path} ({len(out)} feature×outcome rows)", flush=True)

    for y_col in ["y_cv_r2", "y_pass_01", "y_pass_03"]:
        sub = out[out["y"] == y_col].copy()
        if "auc_distance" in sub.columns and sub["auc_distance"].notna().any():
            sub = sub.sort_values("auc_distance", ascending=False)
        else:
            sub = sub.sort_values("spearman_p")
        print(f"\nTop 12 features for {y_col}:", flush=True)
        print(sub.head(12).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
