#!/usr/bin/env python
"""What each ADMIT FAIL label actually measures, with values."""
import pandas as pd, numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
ADM = f"{ROOT}/analyses/sv_catalog/output/ref_ukbb/admit/admit_decompose_cv01.tsv"
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"

a = pd.read_csv(ADM, sep="\t")
print("=== columns available ===")
print([c for c in a.columns])
cols = [c for c in ["n_weight_snps", "n_matched", "ld_match_frac", "n_in_block", "cv_r2"] if c in a.columns]

print("\n=== 1. the LD-panel quantities that DEFINE the labels (median [IQR]) ===")
for st, g in a.groupby("status"):
    bits = []
    for c in cols:
        v = pd.to_numeric(g[c], errors="coerce")
        bits.append(f"{c}={v.median():.3g} [{v.quantile(.25):.3g}-{v.quantile(.75):.3g}]")
    print(f"  {st:24s} n={len(g):5d}  " + "  ".join(bits))

print("\n=== 2. how coarse is ld_match_frac? distribution of n_weight_snps (its denominator) ===")
if "n_weight_snps" in a.columns:
    w = pd.to_numeric(a["n_weight_snps"], errors="coerce")
    print(f"  n_weight_snps over all 11,129: median={w.median():.0f}  "
          f"min={w.min():.0f}  max={w.max():.0f}  frac with <=5 SNPs = {(w<=5).mean():.3f}")
    for st, g in a.groupby("status"):
        v = pd.to_numeric(g["n_weight_snps"], errors="coerce")
        m = pd.to_numeric(g["n_matched"], errors="coerce")
        print(f"    {st:24s} n_weight median={v.median():5.0f}   n_matched median={m.median():5.0f}"
              f"   frac n_matched==0: {(m==0).mean():.3f}")

print("\n=== 3. THE DECISIVE COMPARISON: assembly-side quality, ADMIT vs excluded ===")
print("    (these are the only quantities the CDR-QTL actually consumes)")
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
vt["excluded"] = np.where(vt.status == "ADMIT", "ADMIT (kept)", "excluded by the gate")
for grp, g in vt.groupby("status"):
    print(f"  {grp:24s} n={len(g):5d}  n_miss/locus median={g.n_miss.median():.0f} "
          f"(mean {g.n_miss.mean():.2f})   maf median={g.maf.median():.4f}   "
          f"exact-reconstruction={100*(g.fidelity=='exact').mean():.1f}%")
print()
for grp, g in vt.groupby("excluded"):
    print(f"  {grp:22s} n={len(g):5d}  n_miss mean={g.n_miss.mean():.3f}  "
          f"maf median={g.maf.median():.4f}  af<0.05: {100*(g.maf<0.05).mean():.1f}%  "
          f"exact={100*(g.fidelity=='exact').mean():.1f}%")
adm_m, exc_m = vt[vt.status == "ADMIT"].n_miss, vt[vt.status != "ADMIT"].n_miss
from scipy import stats as st_
print(f"\n  missingness ADMIT vs excluded: Mann-Whitney p = {st_.mannwhitneyu(adm_m, exc_m).pvalue:.3g}")
print(f"  total missing cells: ADMIT {int(adm_m.sum())}, excluded {int(exc_m.sum())} "
      f"(out of {len(vt)*188:,} genotype cells)")
