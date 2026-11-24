#!/usr/bin/env python3
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

SNP_WINDOW_FEATURES = [
    "mean_AF_AFR", "mean_AF_AMR", "mean_AF_EAS", "mean_AF_EUR", "mean_AF_SAS",
    "AF_skew_mean_window", "max_FST_pair_window", "n_pop_snps_window",
]

ALL_FEATURES_EXCAT_GENO = [
    "log_svlen", "log_maf", "n_carriers", "n_hom_alt",
    "svtype_INS", "svtype_DEL", "svtype_INV", "svtype_DUP",
    "log_n_window_snps",
    "log_active_HOR_span", "n_HOR_intervals",
    "n_train_haps", "acrocentric", "pericen_len_500000",
    "top_tag_r2", "n_high_ld_r2_05",
    "max_r2_within_50kb", "max_r2_within_500kb", "max_r2_within_5Mb",
    "mean_pair_r2_window", "n_indep_snps_window", "n_window_snps_actual",
    "gc_10kb", "gc_100kb", "gc_1Mb",
    "sat_active", "sat_local_complexity",
    "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
    "af_skew_sv", "max_FST_pair_sv", "total_af_sv",
] + SNP_WINDOW_FEATURES


def main():
    df = pd.read_csv(OUT / "per_sv_features_with_genome.tsv",
                      sep="\t", low_memory=False)
    sub = df.dropna(subset=["cv_r2_best8"]).copy()

    keep = [f for f in ALL_FEATURES_EXCAT_GENO if f in sub.columns]
    print(f"[INFO] {len(sub)} SVs x {len(keep)} features", flush=True)

    print(f"\nSNP-window feature coverage:", flush=True)
    for f in SNP_WINDOW_FEATURES:
        if f in sub.columns:
            print(f"  {f}: {sub[f].notna().mean()*100:.1f}% non-NA", flush=True)

    X1 = sub[keep].astype(float).values
    sub2 = sub.copy()
    for f in SNP_WINDOW_FEATURES:
        if f in sub2.columns:
            med = sub2[f].median()
            sub2[f] = sub2[f].fillna(med)
    X2 = sub2[keep].astype(float).values

    keep_nosnpwin = [f for f in keep if f not in SNP_WINDOW_FEATURES]
    X3 = sub[keep_nosnpwin].astype(float).values

    y = sub["cv_r2_best8"].astype(float).values

    Xtr_1, Xte_1, Xtr_2, Xte_2, Xtr_3, Xte_3, ytr, yte = train_test_split(
        X1, X2, X3, y, test_size=0.2, random_state=42)

    rows = []
    for label, Xtr, Xte, names in [
        ("variant_1_nan_handled", Xtr_1, Xte_1, keep),
        ("variant_2_median_imputed", Xtr_2, Xte_2, keep),
        ("variant_3_dropped_snpwin", Xtr_3, Xte_3, keep_nosnpwin),
    ]:
        m = HistGradientBoostingRegressor(
            max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
            min_samples_leaf=50, l2_regularization=1.0, random_state=42)
        m.fit(Xtr, ytr)
        train = m.score(Xtr, ytr); test = m.score(Xte, yte)
        pi = permutation_importance(m, Xte, yte, n_repeats=20, random_state=42,
                                     n_jobs=2, scoring="r2")
        print(f"\n=== {label} ===", flush=True)
        print(f"  test R^2: {test:.4f}  (train {train:.4f})", flush=True)
        for f, mn, sd in zip(names, pi.importances_mean, pi.importances_std):
            rows.append({
                "variant": label, "feature": f,
                "perm_imp_mean": float(mn),
                "perm_imp_std": float(sd),
                "test_r2": float(test),
            })
        print(f"  SNP-window feature importances:", flush=True)
        for f in SNP_WINDOW_FEATURES:
            if f in names:
                idx = names.index(f)
                print(f"    {f:30s}  {pi.importances_mean[idx]:+.4f} ± "
                      f"{pi.importances_std[idx]:.4f}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "v3_missingness_leakage.tsv", sep="\t", index=False)
    print(f"\n[OK] wrote v3_missingness_leakage.tsv", flush=True)


if __name__ == "__main__":
    main()
