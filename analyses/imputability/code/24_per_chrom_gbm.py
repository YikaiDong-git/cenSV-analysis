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

CHROMS = ["chr1", "chr17", "chr18", "chr19", "chr21"]

FEATURES = [
    "log_svlen", "log_maf", "n_carriers", "n_hom_alt",
    "svtype_INS", "svtype_DEL", "svtype_INV", "svtype_DUP",
    "log_n_window_snps",
    "log_active_HOR_span", "n_HOR_intervals",
    "n_train_haps", "acrocentric",
    "pericen_len_500000",
    "top_tag_r2", "n_high_ld_r2_05",
    "max_r2_within_50kb", "max_r2_within_500kb", "max_r2_within_5Mb",
    "mean_pair_r2_window", "n_indep_snps_window", "n_window_snps_actual",
    "gc_10kb", "gc_100kb", "gc_1Mb",
    "sat_local_complexity",
    "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
    "af_skew_sv", "max_FST_pair_sv", "total_af_sv",
]


def main():
    df = pd.read_csv(OUT / "per_sv_features_with_genome.tsv",
                      sep="\t", low_memory=False)
    print(f"[INFO] loaded {len(df)} SVs total", flush=True)

    rows = []
    for chrom in CHROMS:
        sub = df[(df["chrom"] == chrom) &
                  df["cv_r2_best8"].notna()].copy()
        if len(sub) < 100:
            print(f"[SKIP] {chrom}: only {len(sub)} SVs", flush=True)
            continue
        print(f"\n=== {chrom} ({len(sub)} SVs) ===", flush=True)

        keep = [f for f in FEATURES if f in sub.columns and
                 sub[f].nunique() > 1]
        dropped = [f for f in FEATURES if f in sub.columns and
                    sub[f].nunique() <= 1]
        print(f"  kept {len(keep)} features (dropped {len(dropped)} constants: "
              f"{dropped[:5]}...)", flush=True)
        X = sub[keep].astype(float).values
        y = sub["cv_r2_best8"].astype(float).values

        kf = KFold(n_splits=min(5, len(sub) // 100), shuffle=True, random_state=42)
        fold_imps = []
        fold_test = []
        for tr, te in kf.split(X):
            m = HistGradientBoostingRegressor(
                max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
                min_samples_leaf=20, l2_regularization=1.0, random_state=42)
            m.fit(X[tr], y[tr])
            r2 = m.score(X[te], y[te])
            fold_test.append(r2)
            pi = permutation_importance(m, X[te], y[te], n_repeats=20,
                                         random_state=42, n_jobs=2, scoring="r2")
            fold_imps.append(pi.importances_mean)
        mean_imp = np.mean(np.vstack(fold_imps), axis=0)
        std_imp = np.std(np.vstack(fold_imps), axis=0)
        order = np.argsort(mean_imp)[::-1]

        print(f"  test R^2 5-fold: {np.mean(fold_test):.3f} ± {np.std(fold_test):.3f}",
              flush=True)
        print(f"  top-5 features:", flush=True)
        for k in range(min(5, len(keep))):
            i = order[k]
            print(f"    {k+1}. {keep[i]:30s}  perm_imp = "
                  f"{mean_imp[i]:.3f} ± {std_imp[i]:.3f}", flush=True)
            rows.append({
                "chrom": chrom, "n_sv": len(sub), "rank": k + 1,
                "feature": keep[i],
                "perm_imp_mean": float(mean_imp[i]),
                "perm_imp_std": float(std_imp[i]),
                "test_r2_mean": float(np.mean(fold_test)),
            })

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "v3_per_chrom_top3.tsv", sep="\t", index=False)
    print(f"\n[OK] wrote v3_per_chrom_top3.tsv", flush=True)
    print("\n=== Top-3 per chrom (consistency check) ===", flush=True)
    pivot = out.pivot_table(index="chrom", columns="rank", values="feature",
                              aggfunc="first")
    print(pivot.to_string(), flush=True)


if __name__ == "__main__":
    main()
