#!/usr/bin/env python3
"""GBM feature importance INCLUDING LD features."""
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

BASE_FEATURES = [
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
LD_FEATURES = [
    "top_tag_r2", "n_high_ld_r2_05",
    "max_r2_within_50kb", "max_r2_within_500kb", "max_r2_within_5Mb",
    "mean_pair_r2_window", "n_indep_snps_window",
    "n_window_snps_actual",
]


def fit_one(X_tr, X_te, y_tr, y_te, names, kind):
    if kind == "regression":
        m = HistGradientBoostingRegressor(
            max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
            min_samples_leaf=50, l2_regularization=1.0, random_state=42)
    else:
        m = HistGradientBoostingClassifier(
            max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
            min_samples_leaf=50, l2_regularization=1.0, random_state=42)
    m.fit(X_tr, y_tr)
    train = m.score(X_tr, y_tr)
    test = m.score(X_te, y_te)
    pi = permutation_importance(m, X_te, y_te, n_repeats=10, random_state=42,
                                n_jobs=2,
                                scoring="r2" if kind == "regression" else "roc_auc")
    return train, test, pi.importances_mean, pi.importances_std


def main():
    feat_path = OUT / "per_sv_features_with_ld.tsv"
    if not feat_path.exists():
        print(f"[FAIL] {feat_path} missing — run 08 aggregator first", flush=True)
        sys.exit(1)
    df = pd.read_csv(feat_path, sep="\t", low_memory=False)
    print(f"[INFO] loaded {len(df)} SVs × {df.shape[1]} cols", flush=True)

    sub = df.dropna(subset=["y_cv_r2", "y_pass_03"])
    keep_base = [c for c in BASE_FEATURES if c in sub.columns]
    keep_ld = [c for c in LD_FEATURES if c in sub.columns]
    print(f"[INFO] base features: {len(keep_base)}, LD features: {len(keep_ld)}",
          flush=True)
    if not keep_ld:
        print("[FAIL] no LD features in table — re-check merge step", flush=True)
        sys.exit(1)

    X_base = sub[keep_base].astype(float).values
    X_full = sub[keep_base + keep_ld].astype(float).values
    Xtr_b, Xte_b, Xtr_f, Xte_f, ytr, yte, ytr_c, yte_c = train_test_split(
        X_base, X_full,
        sub["y_cv_r2"].values, sub["y_pass_03"].astype(int).values,
        test_size=0.2, random_state=42,
    )

    rows = []
    for kind, ytr_, yte_, label in [
        ("regression", ytr, yte, "y_cv_r2"),
        ("classification", ytr_c, yte_c, "y_pass_03"),
    ]:
        print(f"\n=== BASELINE (no LD) — {kind} on {label} ===", flush=True)
        tb, ts_b, _, _ = fit_one(Xtr_b, Xte_b, ytr_, yte_, keep_base, kind)
        print(f"  train={tb:.3f} test={ts_b:.3f}", flush=True)
        print(f"\n=== WITH-LD — {kind} on {label} ===", flush=True)
        tf, ts_f, imp_mean, imp_std = fit_one(Xtr_f, Xte_f, ytr_, yte_,
                                              keep_base + keep_ld, kind)
        print(f"  train={tf:.3f} test={ts_f:.3f}  Δtest={ts_f - ts_b:+.3f}", flush=True)
        for f, m, s in zip(keep_base + keep_ld, imp_mean, imp_std):
            rows.append({
                "feature": f, "outcome": label,
                "block": "LD" if f in keep_ld else "base",
                "perm_imp_mean": float(m), "perm_imp_std": float(s),
                "test_score_baseline": float(ts_b),
                "test_score_with_ld": float(ts_f),
                "delta_test": float(ts_f - ts_b),
            })

    out = pd.DataFrame(rows).sort_values(
        ["outcome", "perm_imp_mean"], ascending=[True, False])
    out.to_csv(OUT / "gbm_importances_with_ld.tsv", sep="\t", index=False)
    print(f"\n[OK] wrote {OUT / 'gbm_importances_with_ld.tsv'}", flush=True)

    summary = out.groupby("outcome").first()[
        ["test_score_baseline", "test_score_with_ld", "delta_test"]].reset_index()
    summary.to_csv(OUT / "gbm_with_vs_without_ld.tsv", sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)

    for outcome in out["outcome"].unique():
        print(f"\nTop 15 features for {outcome}:", flush=True)
        print(out[out["outcome"] == outcome].head(15).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
