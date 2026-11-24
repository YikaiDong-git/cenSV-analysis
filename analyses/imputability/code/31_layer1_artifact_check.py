#!/usr/bin/env python3
from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, fisher_exact
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJECT = Path(CENSV_ROOT + "")
OUT = PROJECT / "analyses/imputability/output"

INTERSECT = OUT / "hg002_truth_intersect.tsv"
V3_TABLE = OUT / "per_sv_features_with_genome.tsv"
SURVIVOR_TABLE = OUT / "v3_survivor_bias_perm.tsv"
LAYER1_FEATURES_TSV = OUT / "layer1_feature_matrix.tsv"

RNG_SEED = 42
POS_TOL = 100
SVLEN_RATIO_MIN = 0.5


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.parent / f"{path.name}.tmp.{os.getpid()}"
    df.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)


def wilson_ci(k: int, n: int, conf: float = 0.95) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    z = norm.ppf(0.5 * (1 + conf))
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (centre - half, centre + half)


def annotate_imputability(intersect: pd.DataFrame, v3: pd.DataFrame) -> pd.DataFrame:
    v3 = v3.copy()
    v3["pos"] = v3["pos"].astype(int)
    v3["svlen"] = v3["svlen"].astype(int) if "svlen" in v3.columns else 0
    if "svtype" not in v3.columns:
        v3["svtype"] = np.where(v3.get("svtype_INS", 0) == 1, "INS",
                          np.where(v3.get("svtype_DEL", 0) == 1, "DEL", "OTHER"))

    out_rows = []
    for _, row in intersect.iterrows():
        chrom_match = v3["chrom"] == row["chrom"]
        svtype_match = v3["svtype"] == row["svtype"]
        sub = v3[chrom_match & svtype_match]
        if sub.empty:
            out_rows.append({**row.to_dict(), "imputable": False})
            continue
        pos_diff = (sub["pos"] - int(row["pos"])).abs()
        sub = sub[pos_diff <= POS_TOL]
        if sub.empty:
            out_rows.append({**row.to_dict(), "imputable": False})
            continue
        prod_svlen = row.get("svlen", np.nan)
        if pd.isna(prod_svlen) or row.get("svtype") == "BND":
            out_rows.append({**row.to_dict(), "imputable": bool(not sub.empty)})
            continue
        prod_svlen = int(prod_svlen)
        sub = sub.dropna(subset=["svlen"])
        if sub.empty:
            out_rows.append({**row.to_dict(), "imputable": False})
            continue
        ratio = np.minimum(sub["svlen"], prod_svlen) / np.maximum(sub["svlen"], prod_svlen)
        sub = sub[ratio >= SVLEN_RATIO_MIN]
        out_rows.append({**row.to_dict(), "imputable": bool(not sub.empty)})

    return pd.DataFrame(out_rows)


def cross_tab(annotated: pd.DataFrame) -> pd.DataFrame:
    out = []
    for hor in [True, False]:
        for imp in [True, False]:
            sub = annotated[(annotated["in_active_HOR"] == hor) & (annotated["imputable"] == imp)]
            n = len(sub)
            k = int(sub["truth_match"].sum())
            lo, hi = wilson_ci(k, n)
            out.append({
                "in_active_HOR": hor,
                "imputable": imp,
                "n": n,
                "TP": k,
                "TP_rate": k / n if n else np.nan,
                "wilson_lo": lo,
                "wilson_hi": hi,
            })
    return pd.DataFrame(out)


def newcombe_wilson_diff(k1: int, n1: int, k2: int, n2: int, conf: float = 0.95) -> tuple[float, float, float]:
    if n1 == 0 or n2 == 0:
        return (np.nan, np.nan, np.nan)
    p1, p2 = k1 / n1, k2 / n2
    z = norm.ppf(0.5 * (1 + conf))
    def wilson_endpoints(k, n):
        denom = 1 + z * z / n
        centre = (k / n + z * z / (2 * n)) / denom
        half = z * np.sqrt(k * (n - k) / (n ** 3) + z * z / (4 * n * n)) / denom
        return (centre - half, centre + half)
    l1, u1 = wilson_endpoints(k1, n1)
    l2, u2 = wilson_endpoints(k2, n2)
    diff = p1 - p2
    ci_lo = diff - np.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    ci_hi = diff + np.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    return (diff, ci_lo, ci_hi)


def primary_contrast(annotated: pd.DataFrame) -> pd.DataFrame:
    hor = annotated[annotated["in_active_HOR"]].copy()
    imp = hor[hor["imputable"]]
    unimp = hor[~hor["imputable"]]
    k1, n1 = int(imp["truth_match"].sum()), len(imp)
    k2, n2 = int(unimp["truth_match"].sum()), len(unimp)

    table = np.array([[k1, n1 - k1], [k2, n2 - k2]])
    if n1 == 0 or n2 == 0:
        odds, p_fisher = np.nan, np.nan
    else:
        odds, p_fisher = fisher_exact(table, alternative="two-sided")

    diff_obs, ci_lo, ci_hi = newcombe_wilson_diff(k1, n1, k2, n2)

    rng = np.random.default_rng(RNG_SEED)
    chroms = hor["chrom"].unique()
    boot_diffs = []
    for _ in range(2000):
        sample_chroms = rng.choice(chroms, size=len(chroms), replace=True)
        sample = pd.concat([hor[hor["chrom"] == c] for c in sample_chroms],
                            ignore_index=True)
        s_imp = sample[sample["imputable"]]
        s_unimp = sample[~sample["imputable"]]
        if len(s_imp) == 0 or len(s_unimp) == 0:
            continue
        boot_diffs.append(s_imp["truth_match"].mean() - s_unimp["truth_match"].mean())
    boot_lo, boot_hi = (np.nan, np.nan)
    if boot_diffs:
        boot_lo = float(np.percentile(boot_diffs, 2.5))
        boot_hi = float(np.percentile(boot_diffs, 97.5))

    return pd.DataFrame([{
        "k_imputable": k1, "n_imputable": n1, "tp_rate_imputable": k1 / n1 if n1 else np.nan,
        "k_unimputable": k2, "n_unimputable": n2, "tp_rate_unimputable": k2 / n2 if n2 else np.nan,
        "risk_diff": diff_obs,
        "ncbm_wilson_lo": ci_lo, "ncbm_wilson_hi": ci_hi,
        "boot_lo_chrom_cluster": boot_lo, "boot_hi_chrom_cluster": boot_hi,
        "fisher_exact_odds": odds, "fisher_exact_p": p_fisher,
    }])


def refit_layer1_tp_only(annotated: pd.DataFrame, v3: pd.DataFrame, k_folds: int = 5) -> dict:
    layer1_features = [
        "gc_10kb", "gc_100kb", "gc_1Mb",
        "log_active_HOR_span", "sat_local_complexity",
        "n_HOR_intervals", "acrocentric",
        "log_svlen", "log_maf", "n_carriers", "n_hom_alt",
    ]
    layer1_features = [f for f in layer1_features if f in v3.columns]

    df = v3.copy()
    if "cv_r2_best8" not in df.columns:
        print("[WARN] cv_r2_best8 missing; using cv_r2 if present", flush=True)
        df["cv_r2_best8"] = df.get("cv_r2", np.nan)

    threshold = np.nanmedian(df["cv_r2_best8"])
    df["y"] = (df["cv_r2_best8"] >= threshold).astype(int)

    X = df[layer1_features].values
    y = df["y"].values
    keep = ~np.isnan(X).any(axis=1) & ~np.isnan(y)
    X, y = X[keep], y[keep]

    aucs = []
    kf = StratifiedKFold(n_splits=k_folds, shuffle=True, random_state=RNG_SEED)
    for tr, te in kf.split(X, y):
        clf = HistGradientBoostingClassifier(
            max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
            min_samples_leaf=50, l2_regularization=1.0,
            random_state=RNG_SEED,
        )
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], p))

    return {
        "n_used": int(keep.sum()),
        "auc_mean": float(np.mean(aucs)),
        "auc_std": float(np.std(aucs)),
        "auc_per_fold": [float(a) for a in aucs],
        "features": layer1_features,
        "threshold_cv_r2": float(threshold),
    }


def main() -> int:
    if not INTERSECT.exists():
        sys.exit(f"[ERR] {INTERSECT} not found — run 30_intersect_hg002_truth.py first")
    if not V3_TABLE.exists():
        sys.exit(f"[ERR] v3 table missing: {V3_TABLE}")

    intersect = pd.read_csv(INTERSECT, sep="\t")
    v3 = pd.read_csv(V3_TABLE, sep="\t")

    annotated = annotate_imputability(intersect, v3)
    out_annot = OUT / "hg002_intersect_with_imputability.tsv"
    atomic_write_tsv(annotated.sort_values(["chrom", "pos", "svtype"]).reset_index(drop=True), out_annot)
    print(f"[OK] annotated -> {out_annot}", flush=True)

    print("\n=== Headline cross-tab (TP rate by imputability × in_active_HOR) ===")
    ct = cross_tab(annotated)
    atomic_write_tsv(ct, OUT / "v3b_crosstab.tsv")
    print(ct.to_string(index=False))

    print("\n=== Primary contrast: TP_rate(imp) vs TP_rate(unimp) in active HOR ===")
    contrast = primary_contrast(annotated)
    atomic_write_tsv(contrast, OUT / "v3b_primary_contrast.tsv")
    print(contrast.to_string(index=False))

    hor = annotated[annotated["in_active_HOR"]]
    perchrom = hor.groupby("chrom").apply(
        lambda g: pd.Series({
            "n": len(g),
            "n_imputable": int(g["imputable"].sum()),
            "TP_imputable": int((g["imputable"] & g["truth_match"]).sum()),
            "TP_rate_imputable": float((g[g["imputable"]]["truth_match"].mean())) if g["imputable"].any() else np.nan,
            "n_unimputable": int((~g["imputable"]).sum()),
            "TP_unimputable": int((~g["imputable"] & g["truth_match"]).sum()),
            "TP_rate_unimputable": float(g[~g["imputable"]]["truth_match"].mean()) if (~g["imputable"]).any() else np.nan,
        }), include_groups=False
    ).reset_index()
    atomic_write_tsv(perchrom, OUT / "v3b_perchrom_tprate.tsv")
    print("\n=== Per-chrom TP rate (active HOR) ===")
    print(perchrom.to_string(index=False))

    print("\n=== Layer-1 refit (proxy on v3 table) ===")
    refit = refit_layer1_tp_only(annotated, v3)
    refit_df = pd.DataFrame([{
        "n_used": refit["n_used"],
        "auc_mean": refit["auc_mean"],
        "auc_std": refit["auc_std"],
        "auc_per_fold": ";".join(f"{x:.6f}" for x in refit["auc_per_fold"]),
        "features": ",".join(refit["features"]),
        "threshold_cv_r2": refit["threshold_cv_r2"],
    }])
    atomic_write_tsv(refit_df, OUT / "v3b_layer1_refit_proxy.tsv")
    print(refit_df.to_string(index=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
