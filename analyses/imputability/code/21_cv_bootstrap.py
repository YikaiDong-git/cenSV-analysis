#!/usr/bin/env python3
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import KFold
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

warnings.filterwarnings("ignore")

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"

OUTCOMES = ["cv_r2_best6", "cv_r2_best8"]

BASE = [
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
CAT = ["cat_monophyletic", "cat_residual",
        "cat_non_decomposed", "cat_original_decomposed"]
LEAKS = ["best_n_nonzero", "best_l1_ratio", "best_snp_dist_mean"]
GENO = [
    "gc_10kb", "gc_100kb", "gc_1Mb",
    "sat_active", "sat_local_complexity",
    "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
    "af_skew_sv", "max_FST_pair_sv", "total_af_sv",
    "mean_AF_AFR", "mean_AF_AMR", "mean_AF_EAS", "mean_AF_EUR", "mean_AF_SAS",
    "AF_skew_mean_window", "max_FST_pair_window", "n_pop_snps_window",
]


def fit_one_fold(X_tr, X_te, y_tr, y_te, names):
    m = HistGradientBoostingRegressor(
        max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
        min_samples_leaf=50, l2_regularization=1.0, random_state=42,
    )
    m.fit(X_tr, y_tr)
    train = m.score(X_tr, y_tr)
    test = m.score(X_te, y_te)
    pi = permutation_importance(m, X_te, y_te, n_repeats=50,
                                random_state=42, n_jobs=2, scoring="r2")
    return train, test, pi.importances_mean, pi.importances_std


def main():
    df = pd.read_csv(OUT / "per_sv_features_with_genome.tsv",
                      sep="\t", low_memory=False)
    print(f"[INFO] loaded {len(df)} SVs", flush=True)

    sat_oh = [c for c in df.columns if c.startswith("sat_") and
              c not in {"sat_active", "sat_local_complexity", "sat_family"}]

    summary_rows = []
    perm_rows = []
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    for outcome in OUTCOMES:
        if outcome not in df.columns:
            continue
        sub = df.dropna(subset=[outcome]).copy()
        if len(sub) < 200:
            continue
        print(f"\n=== {outcome}  (n={len(sub)}) ===", flush=True)

        keep_full = [c for c in BASE if c in sub.columns]
        keep_excat = [c for c in keep_full if c not in CAT and c not in LEAKS]
        full_geno = keep_full + [c for c in GENO if c in sub.columns] + sat_oh
        excat_geno = keep_excat + [c for c in GENO if c in sub.columns] + sat_oh

        modes = {
            "with_category": keep_full,
            "with_category_geno": full_geno,
            "ex_category": keep_excat,
            "ex_category_geno": excat_geno,
        }
        y = sub[outcome].astype(float).values

        for mode, names in modes.items():
            X = sub[names].astype(float).values
            fold_train_r2 = []
            fold_test_r2 = []
            fold_perm_imps = []
            for fold_i, (tr_idx, te_idx) in enumerate(kf.split(X)):
                tr_r2, te_r2, imp_mean, imp_std = fit_one_fold(
                    X[tr_idx], X[te_idx], y[tr_idx], y[te_idx], names)
                fold_train_r2.append(tr_r2)
                fold_test_r2.append(te_r2)
                fold_perm_imps.append(imp_mean)
                print(f"  [{mode} fold{fold_i}] train={tr_r2:.3f} test={te_r2:.3f}",
                      flush=True)
            mean_imp = np.mean(np.vstack(fold_perm_imps), axis=0)
            std_imp = np.std(np.vstack(fold_perm_imps), axis=0)

            te_arr = np.array(fold_test_r2)
            tr_arr = np.array(fold_train_r2)
            summary_rows.append({
                "outcome": outcome, "mode": mode,
                "test_r2_mean_5fold": float(te_arr.mean()),
                "test_r2_std_5fold": float(te_arr.std()),
                "test_r2_min_5fold": float(te_arr.min()),
                "test_r2_max_5fold": float(te_arr.max()),
                "train_r2_mean_5fold": float(tr_arr.mean()),
            })
            for f, mn, sd in zip(names, mean_imp, std_imp):
                perm_rows.append({
                    "outcome": outcome, "mode": mode, "feature": f,
                    "perm_imp_mean_cv": float(mn),
                    "perm_imp_std_cv": float(sd),
                })

    summary = pd.DataFrame(summary_rows)
    perm = pd.DataFrame(perm_rows)
    summary.to_csv(OUT / "v3_cv_bootstrap_summary.tsv", sep="\t", index=False)
    perm.to_csv(OUT / "v3_cv_bootstrap_perm.tsv", sep="\t", index=False)
    print(f"\n[OK] wrote v3_cv_bootstrap_summary.tsv + perm.tsv", flush=True)

    print("\n=== Test R^2 mean ± std (5-fold CV) ===", flush=True)
    print(summary.to_string(index=False), flush=True)

    print("\n=== EX-CATEGORY top-10 (best8, perm-imp_mean ± std across CV folds) ===",
          flush=True)
    sub = perm[(perm["outcome"] == "cv_r2_best8") &
                (perm["mode"] == "ex_category_geno")].sort_values(
                "perm_imp_mean_cv", ascending=False).head(10)
    print(sub.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
