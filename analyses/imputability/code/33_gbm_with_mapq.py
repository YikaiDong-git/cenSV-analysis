#!/usr/bin/env python3
from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJECT = Path(CENSV_ROOT + "")
OUT = PROJECT / "analyses/imputability/output"
MAPQ_DIR = OUT / "mapq_per_chrom"
V3_TABLE = OUT / "per_sv_features_with_genome.tsv"
MERGED_TABLE = OUT / "per_sv_features_v3_with_mapq.tsv"

RNG_SEED = 42
N_FOLDS = 5
N_PERM_REPEATS = 50

MAPQ_FEATURES = [
    "mean_MAPQ_at_locus", "median_MAPQ_at_locus", "n_contigs_at_locus",
    "mean_NM_at_locus", "frac_secondary_at_locus",
]


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.parent / f"{path.name}.tmp.{os.getpid()}"
    df.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)


def merge_mapq_into_v3() -> pd.DataFrame:
    if MERGED_TABLE.exists():
        print(f"[D2] {MERGED_TABLE} exists, loading...", flush=True)
        return pd.read_csv(MERGED_TABLE, sep="\t")

    files = sorted(MAPQ_DIR.glob("mapq_features_*.tsv"))
    if len(files) < 23:
        sys.exit(f"[ERR] only {len(files)}/23 per-chrom MAPQ TSVs present; "
                  f"wait for sbatch jobs to finish.")
    mapq = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    print(f"[D2] merged {len(mapq)} MAPQ rows from {len(files)} chroms", flush=True)

    v3 = pd.read_csv(V3_TABLE, sep="\t")
    print(f"[D2] v3 table: {len(v3)} rows", flush=True)

    if "sv_index" in v3.columns and "sv_index" in mapq.columns:
        merged = v3.merge(mapq[["chrom", "sv_index"] + MAPQ_FEATURES],
                          on=["chrom", "sv_index"], how="left")
    else:
        merged = v3.merge(mapq[["chrom", "pos"] + MAPQ_FEATURES],
                          on=["chrom", "pos"], how="left")
    print(f"[D2] merged table: {len(merged)} rows ({merged[MAPQ_FEATURES[0]].notna().sum()} with MAPQ)", flush=True)

    atomic_write_tsv(merged, MERGED_TABLE)
    return merged


def fit_layer2_cv(merged: pd.DataFrame, mode: str) -> tuple[dict, pd.DataFrame]:
    base = [
        "log_svlen", "log_maf", "n_carriers", "n_hom_alt",
        "svtype_INS", "svtype_DEL",
        "log_active_HOR_span", "n_HOR_intervals", "n_train_haps", "acrocentric",
        "top_tag_r2", "n_high_ld_r2_05",
        "max_r2_within_50kb", "max_r2_within_500kb", "max_r2_within_5Mb",
        "mean_pair_r2_window", "n_indep_snps_window", "n_window_snps_actual",
    ]
    cat = [
        "cat_monophyletic", "cat_residual", "cat_non_decomposed", "cat_original_decomposed",
        "best_n_nonzero", "best_l1_ratio", "best_snp_dist_mean",
    ]
    geno = [
        "gc_10kb", "gc_100kb", "gc_1Mb",
        "sat_local_complexity",
        "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
        "af_skew_sv", "max_FST_pair_sv", "total_af_sv",
        "mean_AF_AFR", "mean_AF_AMR", "mean_AF_EAS", "mean_AF_EUR", "mean_AF_SAS",
        "AF_skew_mean_window", "max_FST_pair_window", "n_pop_snps_window",
    ]
    mapq = [f for f in MAPQ_FEATURES if f in merged.columns]

    if mode == "ex_category_geno_mapq":
        features = base + geno + mapq
    elif mode == "ex_category_geno":
        features = base + geno
    elif mode == "with_category_geno_mapq":
        features = base + cat + geno + mapq
    else:
        raise ValueError(mode)

    features = [f for f in features if f in merged.columns]
    y = merged["cv_r2_best8"].values
    X = merged[features].values

    rng = np.random.default_rng(RNG_SEED)
    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=RNG_SEED)
    test_r2s, train_r2s = [], []
    perm_means = np.zeros(len(features))
    perm_stds = np.zeros(len(features))

    for tr, te in kf.split(X):
        gbm = HistGradientBoostingRegressor(
            max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
            min_samples_leaf=50, l2_regularization=1.0, random_state=RNG_SEED,
        )
        gbm.fit(X[tr], y[tr])
        train_r2s.append(r2_score(y[tr], gbm.predict(X[tr])))
        test_r2s.append(r2_score(y[te], gbm.predict(X[te])))
        pi = permutation_importance(gbm, X[te], y[te],
                                     n_repeats=N_PERM_REPEATS,
                                     random_state=RNG_SEED, n_jobs=2)
        perm_means += pi.importances_mean / N_FOLDS
        perm_stds += pi.importances_std / N_FOLDS

    summary = {
        "mode": mode,
        "n_features": len(features),
        "test_r2_mean": float(np.mean(test_r2s)),
        "test_r2_std": float(np.std(test_r2s)),
        "train_r2_mean": float(np.mean(train_r2s)),
    }
    perm_df = pd.DataFrame({
        "mode": mode, "feature": features,
        "perm_imp_mean": perm_means, "perm_imp_std": perm_stds,
    }).sort_values("perm_imp_mean", ascending=False).reset_index(drop=True)
    return summary, perm_df


def grouped_drop_column_deltaR2(merged: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(RNG_SEED)

    base = [
        "log_svlen", "log_maf", "n_carriers", "n_hom_alt",
        "svtype_INS", "svtype_DEL",
        "n_train_haps", "acrocentric",
        "top_tag_r2", "n_high_ld_r2_05",
        "max_r2_within_50kb", "max_r2_within_500kb", "max_r2_within_5Mb",
        "mean_pair_r2_window", "n_indep_snps_window", "n_window_snps_actual",
        "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
        "af_skew_sv", "max_FST_pair_sv", "total_af_sv",
    ]
    mapq_grp = [f for f in MAPQ_FEATURES if f in merged.columns]
    gc_grp = ["gc_10kb", "gc_100kb", "gc_1Mb"]
    hor_grp = ["log_active_HOR_span", "sat_local_complexity", "n_HOR_intervals"]

    groups = {
        "drop_MAPQ":      (base + gc_grp + hor_grp, base + gc_grp + hor_grp + mapq_grp),
        "drop_GC":        (base + mapq_grp + hor_grp, base + gc_grp + hor_grp + mapq_grp),
        "drop_HOR":       (base + gc_grp + mapq_grp, base + gc_grp + hor_grp + mapq_grp),
        "drop_MAPQ_GC_HOR": (base, base + gc_grp + hor_grp + mapq_grp),
    }

    y = merged["cv_r2_best8"].values
    rows = []
    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=RNG_SEED)
    for name, (no_grp, full) in groups.items():
        full_use = [f for f in full if f in merged.columns]
        no_use = [f for f in no_grp if f in merged.columns]
        X_full = merged[full_use].values
        X_no = merged[no_use].values
        full_r2s, no_r2s = [], []
        for tr, te in kf.split(X_full):
            for X, store in [(X_full, full_r2s), (X_no, no_r2s)]:
                gbm = HistGradientBoostingRegressor(
                    max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
                    min_samples_leaf=50, l2_regularization=1.0,
                    random_state=RNG_SEED,
                )
                gbm.fit(X[tr], y[tr])
                store.append(r2_score(y[te], gbm.predict(X[te])))
        rows.append({
            "group": name,
            "r2_full": float(np.mean(full_r2s)),
            "r2_full_std": float(np.std(full_r2s)),
            "r2_dropped": float(np.mean(no_r2s)),
            "r2_dropped_std": float(np.std(no_r2s)),
            "delta_r2": float(np.mean(full_r2s) - np.mean(no_r2s)),
        })
    return pd.DataFrame(rows)


def partial_correlation(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for f in ["gc_10kb", "gc_100kb", "gc_1Mb", "log_active_HOR_span",
              "sat_local_complexity", "n_hom_alt"]:
        if f not in merged.columns:
            continue
        for z in MAPQ_FEATURES:
            if z not in merged.columns:
                continue
            sub = merged[[f, "cv_r2_best8", z]].dropna()
            if len(sub) < 100:
                continue
            r_fy, _ = pearsonr(sub[f], sub["cv_r2_best8"])
            r_fz, _ = pearsonr(sub[f], sub[z])
            r_yz, _ = pearsonr(sub["cv_r2_best8"], sub[z])
            denom = np.sqrt(max(1e-12, (1 - r_fz ** 2) * (1 - r_yz ** 2)))
            r_partial = (r_fy - r_fz * r_yz) / denom
            rows.append({
                "feature_X": f, "partial_out_Z": z,
                "r_raw_XY": r_fy, "r_partial_XY|Z": r_partial,
                "r_XZ": r_fz, "r_YZ": r_yz,
                "drop_pct": 100 * (1 - abs(r_partial) / max(abs(r_fy), 1e-12)),
                "n": len(sub),
            })
    return pd.DataFrame(rows).sort_values(["feature_X", "drop_pct"], ascending=[True, False])


def main() -> int:
    merged = merge_mapq_into_v3()

    if merged[MAPQ_FEATURES[0]].notna().sum() < 1000:
        sys.exit(f"[ERR] only {merged[MAPQ_FEATURES[0]].notna().sum()} rows have MAPQ; "
                  "wait for sbatch jobs to finish.")

    print("[D3] fitting Layer-2 GBMs with MAPQ features...", flush=True)
    all_summaries, all_perms = [], []
    for mode in ["ex_category_geno", "ex_category_geno_mapq", "with_category_geno_mapq"]:
        print(f"  mode = {mode}", flush=True)
        summary, perm_df = fit_layer2_cv(merged, mode)
        all_summaries.append(summary)
        all_perms.append(perm_df)
        print(f"    test_R2 = {summary['test_r2_mean']:.4f} ± {summary['test_r2_std']:.4f}", flush=True)

    summary_df = pd.DataFrame(all_summaries)
    atomic_write_tsv(summary_df, OUT / "layer2_with_mapq_summary.tsv")
    perm_df = pd.concat(all_perms, ignore_index=True)
    atomic_write_tsv(perm_df, OUT / "layer2_with_mapq_perm.tsv")

    print("[D4a] partial correlation (weaker, single-variable)...", flush=True)
    pcor = partial_correlation(merged)
    atomic_write_tsv(pcor, OUT / "v3b_mapq_partial_cor.tsv")
    print(pcor.head(20).to_string(index=False), flush=True)

    print("[D4b] grouped drop-column ΔR² (Codex D-fix)...", flush=True)
    drop_df = grouped_drop_column_deltaR2(merged)
    atomic_write_tsv(drop_df, OUT / "grouped_dropcol_dR2.tsv")
    print(drop_df.to_string(index=False), flush=True)

    print("\n[OK] all v3b GBM outputs written.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
