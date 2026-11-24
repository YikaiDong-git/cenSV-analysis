#!/usr/bin/env python3
"""Phase C — per-model GBM with genome features added."""
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

BASE_FEATURES = [
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

GENOME_FEATURES_NUM = [
    "gc_10kb", "gc_100kb", "gc_1Mb",
    "sat_active", "sat_local_complexity",
    "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
    "af_skew_sv", "max_FST_pair_sv", "total_af_sv",
    "mean_AF_AFR", "mean_AF_AMR", "mean_AF_EAS", "mean_AF_EUR", "mean_AF_SAS",
    "AF_skew_mean_window", "max_FST_pair_window", "n_pop_snps_window",
]


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
    print(f"  [{label}] train={train:.3f} test={test:.3f}", flush=True)
    return train, test, pi.importances_mean, pi.importances_std


def main():
    feat = pd.read_csv(OUT / "per_sv_features_with_genome.tsv",
                        sep="\t", low_memory=False)
    print(f"[INFO] loaded {len(feat)} SVs x {feat.shape[1]} cols", flush=True)

    sat_oh = [c for c in feat.columns if c.startswith("sat_") and
              c not in {"sat_active", "sat_local_complexity", "sat_family"}]
    GENOME_ALL = GENOME_FEATURES_NUM + sat_oh
    print(f"[INFO] genome features: {len(GENOME_ALL)} ({sat_oh})", flush=True)

    rows = []
    summary = []
    for outcome in OUTCOMES:
        if outcome not in feat.columns:
            continue
        sub = feat.dropna(subset=[outcome]).copy()
        if len(sub) < 200:
            continue
        print(f"\n=== outcome: {outcome}  (n={len(sub)}) ===", flush=True)

        keep_full = [c for c in BASE_FEATURES if c in sub.columns]
        keep_excat = [c for c in keep_full
                      if c not in CAT_COLS and c not in POST_EN_LEAKS]
        keep_full_g = keep_full + [c for c in GENOME_ALL if c in sub.columns]
        keep_excat_g = keep_excat + [c for c in GENOME_ALL if c in sub.columns]

        Xs = {
            "with_category": sub[keep_full].astype(float).values,
            "with_category_geno": sub[keep_full_g].astype(float).values,
            "ex_category": sub[keep_excat].astype(float).values,
            "ex_category_geno": sub[keep_excat_g].astype(float).values,
        }
        names_map = {
            "with_category": keep_full, "with_category_geno": keep_full_g,
            "ex_category": keep_excat, "ex_category_geno": keep_excat_g,
        }
        y = sub[outcome].astype(float).values

        Xtr_a, Xte_a, Xtr_b, Xte_b, Xtr_c, Xte_c, Xtr_d, Xte_d, ytr, yte = \
            train_test_split(
                Xs["with_category"], Xs["with_category_geno"],
                Xs["ex_category"], Xs["ex_category_geno"],
                y, test_size=0.2, random_state=42)

        Xs_tr = {"with_category": Xtr_a, "with_category_geno": Xtr_b,
                  "ex_category": Xtr_c, "ex_category_geno": Xtr_d}
        Xs_te = {"with_category": Xte_a, "with_category_geno": Xte_b,
                  "ex_category": Xte_c, "ex_category_geno": Xte_d}

        scores = {}
        for mode in ["with_category", "with_category_geno",
                      "ex_category", "ex_category_geno"]:
            train_r2, test_r2, imp_mean, imp_std = fit_one(
                Xs_tr[mode], Xs_te[mode], ytr, yte, names_map[mode],
                f"{outcome}/{mode}")
            scores[mode] = test_r2
            for f, mn, sd in zip(names_map[mode], imp_mean, imp_std):
                rows.append({
                    "outcome": outcome, "mode": mode, "feature": f,
                    "perm_imp_mean": float(mn), "perm_imp_std": float(sd),
                    "train_r2": float(train_r2), "test_r2": float(test_r2),
                    "n": len(sub),
                })
        summary.append({
            "outcome": outcome,
            "with_cat": scores["with_category"],
            "with_cat_geno": scores["with_category_geno"],
            "delta_with_cat": scores["with_category_geno"] - scores["with_category"],
            "ex_cat": scores["ex_category"],
            "ex_cat_geno": scores["ex_category_geno"],
            "delta_ex_cat": scores["ex_category_geno"] - scores["ex_category"],
        })

    out_long = pd.DataFrame(rows).sort_values(
        ["outcome", "mode", "perm_imp_mean"], ascending=[True, True, False])
    out_long.to_csv(OUT / "v3_per_model_with_genome.tsv", sep="\t", index=False)
    print(f"\n[OK] wrote {OUT/'v3_per_model_with_genome.tsv'} ({len(out_long)} rows)",
          flush=True)
    out_summary = pd.DataFrame(summary)
    out_summary.to_csv(OUT / "v3_genome_delta_summary.tsv", sep="\t", index=False)
    print(f"[OK] wrote {OUT/'v3_genome_delta_summary.tsv'}", flush=True)
    print("\nGenome-feature delta test R^2 summary:", flush=True)
    print(out_summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
