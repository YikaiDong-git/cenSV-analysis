#!/usr/bin/env python3
"""Summarize + empirically calibrate the per-SV FinnGen burden grid."""
import sys
import math
import numpy as np
import pandas as pd

df = pd.read_csv(sys.argv[1], sep="\t")
out_rows = []
GWS = 5e-8

for sv in sorted(df["sv"].unique()):
    d = df[df["sv"] == sv]
    null = d[d["cell"] == "NEG_CTRL_unrelated"]["z_SV"].to_numpy()
    null = null[np.isfinite(null)]
    mu, sd = float(np.mean(null)), float(np.std(null, ddof=1))
    med = float(np.median(null))
    mad = float(np.median(np.abs(null - med)))
    rsd = 1.4826 * mad if mad > 0 else sd
    cv_r2 = float(d["cv_r2"].iloc[0])
    n_gws_null = int((d[d["cell"] == "NEG_CTRL_unrelated"]["region_n_gws"] > 0).sum())
    print(f"\n================  {sv}  (cv_r2={cv_r2:.3f})  ================")
    print(f"EMPIRICAL NULL over {len(null)} unrelated FinnGen traits: "
          f"mean={mu:+.3f} sd={sd:.3f} | robust median={med:+.3f} robust_sd={rsd:.3f} "
          f"| min={null.min():+.2f} max={null.max():+.2f} | |z|>=2: {int((np.abs(null)>=2).sum())} "
          f"|z|>=3: {int((np.abs(null)>=3).sum())}")
    print(f"  (SD>1 => burden z inflated vs N(0,1) under FinnGen LD; use robust_sd to recalibrate)")

    for cellname in ("TARGET", "NEG_CTRL_crosstrait"):
        sub = d[d["cell"] == cellname].copy()
        if len(sub) == 0:
            continue
        sub = sub.sort_values("p_SV")
        print(f"\n  --- {cellname} ---")
        hdr = f"  {'endpoint':<18}{'role':<9}{'n_tag':>6}{'z_SV':>8}{'p_raw':>10}{'z_cal':>8}{'p_cal':>10}{'emp_p':>9}{'reg_minp':>11}{'reg_GWS':>8}"
        print(hdr)
        for _, r in sub.iterrows():
            z = r["z_SV"]
            z_cal = (z - med) / rsd if np.isfinite(z) else float("nan")
            p_cal = math.erfc(abs(z_cal) / math.sqrt(2)) if np.isfinite(z_cal) else float("nan")
            emp_p = (1 + int((np.abs(null - med) >= abs(z - med)).sum())) / (1 + len(null)) if np.isfinite(z) else float("nan")
            print(f"  {r['endpoint']:<18}{str(r['role']):<9}{int(r['n_tag_match']):>6}{z:>8.3f}"
                  f"{r['p_SV']:>10.4g}{z_cal:>8.3f}{p_cal:>10.4g}{emp_p:>9.4g}"
                  f"{r['region_min_p']:>11.3g}{int(r['region_n_gws']):>8}")
            out_rows.append(dict(
                sv=sv, cv_r2=cv_r2, endpoint=r["endpoint"], panel=r["panel"], role=r["role"], cell=cellname,
                n_tag_match=int(r["n_tag_match"]), z_SV=round(float(z), 4), p_raw=float(r["p_SV"]),
                null_median=round(med, 4), null_robust_sd=round(rsd, 4),
                z_cal=round(float(z_cal), 4), p_cal=float(p_cal), emp_p=float(emp_p),
                region_min_p=float(r["region_min_p"]), region_n_gws=int(r["region_n_gws"]),
                region_top_snp=r["region_top_snp"], region_top_p=float(r["region_top_p"])))

pd.DataFrame(out_rows).to_csv(sys.argv[2], sep="\t", index=False, float_format="%.4g")
print(f"\n[wrote] {sys.argv[2]}")
