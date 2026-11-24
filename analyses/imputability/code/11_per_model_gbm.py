#!/usr/bin/env python3
"""Phase A2 — per-model GBM x {WITH-CATEGORY, EX-CATEGORY}."""
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

warnings.filterwarnings("ignore")

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"

OUTCOMES = [
    "cv_r2_best6", "cv_r2_best8",
    "cv_r2_edge500", "cv_r2_edge1000", "cv_r2_edge2000",
    "cv_r2_fksize100k", "cv_r2_fksize200k", "cv_r2_fksize500k",
    "cv_r2_midpoint500", "cv_r2_midpoint1000",
]

FEATURES_FULL = [
    "log_svlen", "log_maf", "n_carriers", "n_hom_alt",
    "svtype_INS", "svtype_DEL", "svtype_INV", "svtype_DUP",
    "cat_monophyletic", "cat_residual",
    "cat_non_decomposed", "cat_original_decomposed",
    "log_n_window_snps", "best_snp_dist_mean",
    "log_active_HOR_span", "n_HOR_intervals",
    "n_train_haps", "acrocentric",
    "best_n_nonzero", "best_l1_ratio",
    "pericen_len_500000",
    "top_tag_r2", "n_high_ld_r2_05",
    "max_r2_within_50kb", "max_r2_within_500kb", "max_r2_within_5Mb",
    "mean_pair_r2_window", "n_indep_snps_window", "n_window_snps_actual",
]

CAT_COLS = ["cat_monophyletic", "cat_residual",
            "cat_non_decomposed", "cat_original_decomposed"]

POST_EN_LEAKS = ["best_n_nonzero", "best_l1_ratio", "best_snp_dist_mean"]


def fit_one(X_tr, X_te, y_tr, y_te, names, label):
    m = HistGradientBoostingRegressor(
        max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
        min_samples_leaf=50, l2_regularization=1.0, random_state=42,
    )
    m.fit(X_tr, y_tr)
    train = m.score(X_tr, y_tr)
    test = m.score(X_te, y_te)
    pi = permutation_importance(m, X_te, y_te, n_repeats=10, random_state=42,
                                n_jobs=2, scoring="r2")
    print(f"  [{label}] n_tr={len(X_tr)} train={train:.3f} test={test:.3f}",
          flush=True)
    return train, test, pi.importances_mean, pi.importances_std


def main():
    feat_path = OUT / "per_sv_features_v3.tsv"
    if not feat_path.exists():
        print(f"[FAIL] {feat_path} missing — run 10_build_features.py first",
              flush=True)
        sys.exit(1)
    df = pd.read_csv(feat_path, sep="\t", low_memory=False)
    print(f"[INFO] loaded {len(df)} SVs x {df.shape[1]} cols", flush=True)

    rows = []
    summary = []
    for outcome in OUTCOMES:
        if outcome not in df.columns:
            print(f"[SKIP] {outcome}: column missing", flush=True)
            continue
        sub = df.dropna(subset=[outcome]).copy()
        if len(sub) < 200:
            print(f"[SKIP] {outcome}: only {len(sub)} valid rows", flush=True)
            continue
        print(f"\n=== outcome: {outcome}  (n={len(sub)}) ===", flush=True)

        keep_full = [c for c in FEATURES_FULL if c in sub.columns]
        keep_excat = [c for c in keep_full
                      if c not in CAT_COLS and c not in POST_EN_LEAKS]

        X_full = sub[keep_full].astype(float).values
        X_excat = sub[keep_excat].astype(float).values
        y = sub[outcome].astype(float).values

        Xtr_f, Xte_f, Xtr_e, Xte_e, ytr, yte = train_test_split(
            X_full, X_excat, y, test_size=0.2, random_state=42)

        for mode, Xtr, Xte, names in [
            ("with_category", Xtr_f, Xte_f, keep_full),
            ("ex_category",   Xtr_e, Xte_e, keep_excat),
        ]:
            train_r2, test_r2, imp_mean, imp_std = fit_one(
                Xtr, Xte, ytr, yte, names, f"{outcome}/{mode}")
            for f, mn, sd in zip(names, imp_mean, imp_std):
                rows.append({
                    "outcome": outcome, "mode": mode, "feature": f,
                    "perm_imp_mean": float(mn), "perm_imp_std": float(sd),
                    "train_r2": float(train_r2), "test_r2": float(test_r2),
                    "n": len(sub),
                })

        with_test = next(r["test_r2"] for r in rows
                         if r["outcome"] == outcome and r["mode"] == "with_category")
        ex_test = next(r["test_r2"] for r in rows
                       if r["outcome"] == outcome and r["mode"] == "ex_category")
        ex_imp = sorted(
            [r for r in rows if r["outcome"] == outcome and r["mode"] == "ex_category"],
            key=lambda r: -r["perm_imp_mean"])[:5]
        summary.append({
            "outcome": outcome,
            "test_r2_with_category": with_test,
            "test_r2_ex_category": ex_test,
            "delta_test_r2": with_test - ex_test,
            "ex_top1": ex_imp[0]["feature"] if len(ex_imp) > 0 else "",
            "ex_top1_imp": ex_imp[0]["perm_imp_mean"] if len(ex_imp) > 0 else np.nan,
            "ex_top2": ex_imp[1]["feature"] if len(ex_imp) > 1 else "",
            "ex_top2_imp": ex_imp[1]["perm_imp_mean"] if len(ex_imp) > 1 else np.nan,
            "ex_top3": ex_imp[2]["feature"] if len(ex_imp) > 2 else "",
            "ex_top3_imp": ex_imp[2]["perm_imp_mean"] if len(ex_imp) > 2 else np.nan,
            "ex_top4": ex_imp[3]["feature"] if len(ex_imp) > 3 else "",
            "ex_top5": ex_imp[4]["feature"] if len(ex_imp) > 4 else "",
        })

    out_long = pd.DataFrame(rows).sort_values(
        ["outcome", "mode", "perm_imp_mean"], ascending=[True, True, False])
    out_long.to_csv(OUT / "v3_per_model_with_ex_cat.tsv", sep="\t", index=False)
    print(f"\n[OK] wrote {OUT/'v3_per_model_with_ex_cat.tsv'} ({len(out_long)} rows)",
          flush=True)

    out_summary = pd.DataFrame(summary)
    out_summary.to_csv(OUT / "v3_with_vs_excat_summary.tsv", sep="\t", index=False)
    print(f"[OK] wrote {OUT/'v3_with_vs_excat_summary.tsv'}", flush=True)
    print("\nPer-outcome with vs ex-category dtest R^2 + EX-CAT top-5:", flush=True)
    print(out_summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
