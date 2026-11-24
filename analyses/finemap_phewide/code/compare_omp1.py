#!/usr/bin/env python3
import pandas as pd
ok = True
for arm in ["impute", "matched"]:
    fu = pd.read_csv(f"output/validate_chr16_{arm}.tsv", sep="\t"); fu = fu[fu["mode"] == "hitbundle"]
    ts = pd.read_csv(f"output/twostep_omp1_{arm}.tsv", sep="\t"); ts = ts[ts["mode"] == "hitbundle"]
    m = fu.merge(ts, on=["sv_locus", "trait"], suffixes=("_f", "_t"))
    if len(m) == 0:
        print(f"{arm}: NO OVERLAP"); ok = False; continue
    dp = (m["sv_pip_f"] - m["sv_pip_t"]).abs()
    dr = (m["svSNP_r_f"] - m["svSNP_r_t"]).abs()
    dz = (m["z_SV_f"] - m["z_SV_t"]).abs()
    dcs = (m["in_cs_f"] != m["in_cs_t"]).sum()
    print(f"{arm}: n={len(m)} max|dPIP|={dp.max():.2e} mean|dPIP|={dp.mean():.2e} "
          f"max|dsvSNPr|={dr.max():.2e} max|dz|={dz.max():.2e} in_cs_mismatch={dcs}")
    if dp.max() > 1e-4:
        ok = False
        w = m.loc[dp.idxmax()]
        print(f"   worst: {w['sv_locus']} {w['trait']} fused={w['sv_pip_f']:.6f} two-step={w['sv_pip_t']:.6f}")
print("OMP=1 MATCHED THREADING:", "BIT-FAITHFUL (fused == two-step)" if ok else "residual > 1e-4 (inspect)")
