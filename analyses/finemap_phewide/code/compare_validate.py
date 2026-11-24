#!/usr/bin/env python3
import glob
import pandas as pd

KI = "output/keyinterest"
ok = True
for arm in ["impute", "matched"]:
    fu = pd.read_csv(f"output/validate_chr16_{arm}.tsv", sep="\t")
    fu = fu[fu["mode"] == "hitbundle"]
    kis = glob.glob(f"{KI}/chr16*__{arm}__diag.tsv")
    ki = pd.concat([pd.read_csv(f, sep="\t") for f in kis], ignore_index=True)
    ki = ki[ki["mode"] == "hitbundle"]
    m = fu.merge(ki, on=["sv_locus", "trait"], suffixes=("_f", "_k"))
    if len(m) == 0:
        print(f"{arm}: NO OVERLAP (fused SVs={sorted(fu['sv_locus'].unique())[:3]})"); ok = False; continue
    dp = (m["sv_pip_f"] - m["sv_pip_k"]).abs()
    dr = (m["svSNP_r_f"] - m["svSNP_r_k"]).abs()
    dz = (m["z_SV_f"] - m["z_SV_k"]).abs()
    dcs = (m["in_cs_f"] != m["in_cs_k"]).sum()
    print(f"{arm}: n={len(m)} cells, SVs={sorted(m['sv_locus'].str.split('_').str[:2].str.join('_').unique())} "
          f"max|dPIP|={dp.max():.2e} max|dsvSNPr|={dr.max():.2e} max|dz|={dz.max():.2e} in_cs_mismatch={dcs}")
    if dp.max() > 1e-3 or dr.max() > 1e-3 or dz.max() > 1e-3:
        ok = False
        w = m.loc[dp.idxmax()]
        print(f"   WORST PIP: {w['sv_locus']} {w['trait']} fused={w['sv_pip_f']:.5f} two-step={w['sv_pip_k']:.5f}")
        w = m.loc[dr.idxmax()]
        print(f"   WORST svSNP_r: {w['sv_locus']} {w['trait']} fused={w['svSNP_r_f']:.5f} two-step={w['svSNP_r_k']:.5f}")
print("VALIDATION", "PASS (fused == two-step)" if ok else "FAIL (divergence > 1e-3)")
