#!/usr/bin/env python3
"""Per-strategy GBM feature importance."""
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

STRATEGIES = ["edge500", "edge1000", "edge2000",
              "fksize100k", "fksize200k", "fksize500k",
              "midpoint500", "midpoint1000"]

FEATURES = [
    "log_svlen", "log_maf",
    "n_carriers", "n_hom_alt",
    "svtype_INS", "svtype_DEL", "svtype_INV", "svtype_DUP",
    "cat_monophyletic", "cat_residual",
    "cat_non_decomposed", "cat_original_decomposed",
    "log_n_window_snps", "best_snp_dist_mean",
    "log_active_HOR_span", "n_HOR_intervals",
    "n_train_haps", "acrocentric",
    "pericen_len_500000",
]


def main():
    df = pd.read_csv(OUT / "per_sv_features.tsv", sep="\t", low_memory=False)
    print(f"[INFO] loaded {len(df)} SVs", flush=True)

    keep = [c for c in FEATURES if c in df.columns]
    rows = []
    summary_rows = []

    for strat in STRATEGIES:
        y_col = f"cv_r2_{strat}"
        if y_col not in df.columns:
            print(f"[SKIP] {strat}: no {y_col} column", flush=True)
            continue
        sub = df[keep + [y_col]].dropna(subset=[y_col]).copy()
        if len(sub) < 200:
            print(f"[SKIP] {strat}: only {len(sub)} valid rows", flush=True)
            continue
        X = sub[keep].astype(float).values
        y = sub[y_col].values
        Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=42)
        m = HistGradientBoostingRegressor(
            max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
            min_samples_leaf=50, l2_regularization=1.0, random_state=42,
        )
        m.fit(Xtr, ytr)
        train_r2 = m.score(Xtr, ytr)
        test_r2 = m.score(Xte, yte)
        print(f"[{strat}] n={len(sub)} train={train_r2:.3f} test={test_r2:.3f}",
              flush=True)
        pi = permutation_importance(m, Xte, yte, n_repeats=10, random_state=42,
                                    n_jobs=2, scoring="r2")
        for f, mn, sd in zip(keep, pi.importances_mean, pi.importances_std):
            rows.append({
                "strategy": strat, "feature": f,
                "perm_imp_mean": float(mn), "perm_imp_std": float(sd),
                "train_r2": float(train_r2), "test_r2": float(test_r2),
                "n": len(sub),
            })
        order = np.argsort(pi.importances_mean)[::-1][:5]
        top5 = [keep[i] for i in order]
        summary_rows.append({
            "strategy": strat, "n_svs": len(sub),
            "test_r2": float(test_r2),
            "top1": top5[0] if len(top5) > 0 else "",
            "top2": top5[1] if len(top5) > 1 else "",
            "top3": top5[2] if len(top5) > 2 else "",
            "top4": top5[3] if len(top5) > 3 else "",
            "top5": top5[4] if len(top5) > 4 else "",
        })

    long_df = pd.DataFrame(rows)
    long_df.to_csv(OUT / "per_strategy_gbm.tsv", sep="\t", index=False)
    print(f"\n[OK] wrote {OUT/'per_strategy_gbm.tsv'} ({len(long_df)} rows)",
          flush=True)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT / "per_strategy_top_features.tsv", sep="\t", index=False)
    print(f"[OK] wrote {OUT/'per_strategy_top_features.tsv'}", flush=True)
    print("\nPer-strategy top-5 features:", flush=True)
    print(summary_df.to_string(index=False), flush=True)

    pivot = long_df.pivot_table(index="feature", columns="strategy",
                                values="perm_imp_mean", aggfunc="mean")
    pivot["_avg"] = pivot.mean(axis=1)
    pivot = pivot.sort_values("_avg", ascending=False).drop(columns="_avg")
    pivot.to_csv(OUT / "per_strategy_gbm_pivot.tsv", sep="\t")
    print(f"[OK] wrote {OUT/'per_strategy_gbm_pivot.tsv'}", flush=True)


if __name__ == "__main__":
    main()
