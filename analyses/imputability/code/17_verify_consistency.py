#!/usr/bin/env python3
"""Phase D verification — quantify EX-CATEGORY ranking consistency across 10 outcomes."""
import sys
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
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


def main():
    df = pd.read_csv(OUT / "v3_per_model_with_ex_cat.tsv", sep="\t",
                      low_memory=False)
    rows = []
    for mode in ["with_category", "ex_category"]:
        sub = df[df["mode"] == mode]
        piv = sub.pivot_table(index="feature", columns="outcome",
                                values="perm_imp_mean", aggfunc="mean")
        piv = piv.reindex(columns=OUTCOMES).dropna(how="all")
        rhos = []
        for o1, o2 in combinations(OUTCOMES, 2):
            if o1 in piv.columns and o2 in piv.columns:
                a = piv[o1].fillna(0); b = piv[o2].fillna(0)
                if a.std() > 0 and b.std() > 0:
                    r, _ = spearmanr(a, b)
                    rhos.append(r)
        rhos = np.array(rhos)
        rows.append({
            "mode": mode,
            "n_pairs": len(rhos),
            "mean_spearman_rho": float(rhos.mean()),
            "median_spearman_rho": float(np.median(rhos)),
            "min_spearman_rho": float(rhos.min()),
            "max_spearman_rho": float(rhos.max()),
            "frac_pairs_above_0.9": float(np.mean(rhos > 0.9)),
        })

        top3 = {}
        for outc in OUTCOMES:
            if outc not in piv.columns:
                continue
            ord_ = piv[outc].sort_values(ascending=False).head(3).index.tolist()
            top3[outc] = ord_
        feat_count = {}
        for outc, t in top3.items():
            for f in t:
                feat_count[f] = feat_count.get(f, 0) + 1
        feat_count_sorted = sorted(feat_count.items(), key=lambda x: -x[1])
        rows[-1]["top3_features_count_top5"] = " | ".join(
            f"{f}:{n}/10" for f, n in feat_count_sorted[:5])

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "v3_consistency_metrics.tsv", sep="\t", index=False)
    print(out.to_string(index=False), flush=True)
    print(f"\n[OK] wrote {OUT/'v3_consistency_metrics.tsv'}", flush=True)


if __name__ == "__main__":
    main()
