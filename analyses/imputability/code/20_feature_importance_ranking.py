#!/usr/bin/env python3
"""Quantitative feature-importance ranking across 10 outcomes x 4 modes."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"

OUTCOMES = [
    "cv_r2_best6", "cv_r2_best8",
    "cv_r2_edge500", "cv_r2_edge1000", "cv_r2_edge2000",
    "cv_r2_fksize100k", "cv_r2_fksize200k", "cv_r2_fksize500k",
    "cv_r2_midpoint500", "cv_r2_midpoint1000",
]


def rank_block(df, mode):
    sub = df[df["mode"] == mode].copy()
    if len(sub) == 0:
        return None
    rows = []
    for feat, grp in sub.groupby("feature"):
        imps = grp.set_index("outcome")["perm_imp_mean"].reindex(OUTCOMES).dropna()
        if len(imps) == 0:
            continue
        rows.append({
            "mode": mode, "feature": feat, "n_outcomes": int(len(imps)),
            "mean_imp": float(imps.mean()),
            "median_imp": float(imps.median()),
            "std_imp": float(imps.std()),
            "min_imp": float(imps.min()),
            "max_imp": float(imps.max()),
            "n_negative": int((imps < 0).sum()),
        })
    top1 = top3 = top5 = {}
    rank_counts = {f: {"top1": 0, "top3": 0, "top5": 0} for f in sub["feature"].unique()}
    for outc in OUTCOMES:
        oc = sub[sub["outcome"] == outc].sort_values("perm_imp_mean", ascending=False)
        feats = oc["feature"].tolist()
        for k, ks in [(1, "top1"), (3, "top3"), (5, "top5")]:
            for f in feats[:k]:
                if f in rank_counts:
                    rank_counts[f][ks] += 1
    for r in rows:
        r["n_top1"] = rank_counts[r["feature"]]["top1"]
        r["n_top3"] = rank_counts[r["feature"]]["top3"]
        r["n_top5"] = rank_counts[r["feature"]]["top5"]
    return pd.DataFrame(rows).sort_values("mean_imp", ascending=False)


def main():
    pre = pd.read_csv(OUT / "v3_per_model_with_ex_cat.tsv", sep="\t",
                      low_memory=False)
    geno = pd.read_csv(OUT / "v3_per_model_with_genome.tsv", sep="\t",
                       low_memory=False)

    blocks = []
    for df, modes in [(pre, ["with_category", "ex_category"]),
                      (geno, ["with_category", "with_category_geno",
                              "ex_category", "ex_category_geno"])]:
        for mode in modes:
            r = rank_block(df, mode)
            if r is not None:
                blocks.append(r)

    full = pd.concat(blocks, ignore_index=True)
    full = full.drop_duplicates(subset=["mode", "feature"])
    full = full.sort_values(["mode", "mean_imp"], ascending=[True, False])
    out_path = OUT / "v3_feature_ranking.tsv"
    full.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {out_path}", flush=True)

    for mode in ["with_category", "with_category_geno", "ex_category",
                 "ex_category_geno"]:
        sub = full[full["mode"] == mode].head(30)
        if len(sub) == 0:
            continue
        print(f"\n=== mode: {mode}  (top 30 by mean perm-imp across 10 outcomes) ===",
              flush=True)
        print(sub[["feature", "mean_imp", "median_imp", "std_imp", "min_imp",
                   "max_imp", "n_negative", "n_top1", "n_top3", "n_top5"]
                ].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
