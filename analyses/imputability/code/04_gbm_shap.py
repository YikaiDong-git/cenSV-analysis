#!/usr/bin/env python3
"""Gradient-boosted regression + permutation importance for CV R²."""
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor, HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

warnings.filterwarnings("ignore")

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"

FEATURES = [
    "log_svlen", "log_maf",
    "n_carriers", "n_hom_alt",
    "svtype_INS", "svtype_DEL", "svtype_INV", "svtype_DUP",
    "cat_monophyletic", "cat_residual",
    "cat_non_decomposed", "cat_original_decomposed",
    "log_n_window_snps", "best_snp_dist_mean",
    "log_active_HOR_span", "n_HOR_intervals",
    "n_train_haps", "acrocentric",
    "best_n_nonzero", "best_l1_ratio",
    "pericen_len_500000",
]


def run_one(X_tr, X_te, y_tr, y_te, names, kind, outcome_label):
    if kind == "regression":
        model = HistGradientBoostingRegressor(
            max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
            min_samples_leaf=50, l2_regularization=1.0, random_state=42)
    else:
        model = HistGradientBoostingClassifier(
            max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
            min_samples_leaf=50, l2_regularization=1.0, random_state=42)
    model.fit(X_tr, y_tr)
    train_score = model.score(X_tr, y_tr)
    test_score = model.score(X_te, y_te)
    print(f"[{outcome_label}] train={train_score:.3f} test={test_score:.3f}", flush=True)

    pi = permutation_importance(
        model, X_te, y_te, n_repeats=10, random_state=42, n_jobs=2,
        scoring="r2" if kind == "regression" else "roc_auc",
    )
    rows = []
    for f, m, s in zip(names, pi.importances_mean, pi.importances_std):
        rows.append({
            "feature": f, "outcome": outcome_label,
            "perm_imp_mean": float(m), "perm_imp_std": float(s),
            "train_score": float(train_score), "test_score": float(test_score),
        })
    return pd.DataFrame(rows)


def main():
    feat_path = OUT / "per_sv_features.tsv"
    df = pd.read_csv(feat_path, sep="\t", low_memory=False)
    print(f"[INFO] loaded {len(df)} SVs", flush=True)

    keep = [c for c in FEATURES if c in df.columns]
    sub = df[keep + ["y_cv_r2", "y_pass_03"]].dropna(subset=["y_cv_r2", "y_pass_03"])
    print(f"[INFO] using {len(sub)} SVs × {len(keep)} features (NaN allowed in X)",
          flush=True)
    nan_rate = sub[keep].isna().mean().sort_values(ascending=False)
    print("Per-feature NaN rate (top 5):", flush=True)
    print(nan_rate.head().to_string(), flush=True)
    X = sub[keep].astype(float).values

    Xtr, Xte, ytr, yte, ytr_c, yte_c = train_test_split(
        X, sub["y_cv_r2"].values, sub["y_pass_03"].astype(int).values,
        test_size=0.2, random_state=42,
    )

    print("\n=== HistGradientBoosting regressor on y_cv_r2 ===", flush=True)
    out_r = run_one(Xtr, Xte, ytr, yte, keep, "regression", "y_cv_r2")
    print("\n=== HistGradientBoosting classifier on y_pass_03 ===", flush=True)
    out_c = run_one(Xtr, Xte, ytr_c, yte_c, keep, "classification", "y_pass_03")

    out = pd.concat([out_r, out_c], ignore_index=True)
    out = out.sort_values(["outcome", "perm_imp_mean"], ascending=[True, False])
    out_path = OUT / "gbm_importances.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\n[OK] wrote {out_path}", flush=True)
    for outcome in out["outcome"].unique():
        sub_o = out[out["outcome"] == outcome].head(15)
        print(f"\nTop 15 features for {outcome}:", flush=True)
        print(sub_o.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
