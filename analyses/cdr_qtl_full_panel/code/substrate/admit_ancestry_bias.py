#!/usr/bin/env python
"""Is the ADMIT gate ancestry-biased?"""
import pandas as pd, numpy as np, os
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
PED = f"{ROOT}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
samp = np.array(open(f"{WU}/geno/full_samples.txt").read().split())
ped = pd.read_csv(PED, sep=r"\s+")
sc = [c for c in ped.columns if c.lower() in ("sampleid", "sample")][0]
spc = [c for c in ped.columns if "super" in c.lower()][0]
spmap = ped.set_index(sc)[spc].to_dict()
sp = np.array([spmap.get(s, "NA") for s in samp])
print("samples per superpopulation:", {p: int((sp == p).sum()) for p in POPS},
      " unassigned:", int((sp == "NA").sum()))

maf = {}
for p in POPS:
    m = sp == p
    sub = G[:, m]
    af = np.nanmean(sub, axis=1) / 2
    maf[p] = np.minimum(af, 1 - af)
    vt[f"maf_{p}"] = maf[p]

vt["grp"] = np.where(vt.status == "ADMIT", "ADMIT (kept)", "excluded")
print("\n=== mean within-superpopulation MAF, kept vs excluded ===")
print(f"{'group':14s} {'n':>6s} " + " ".join(f"{p:>8s}" for p in POPS))
for g, d in vt.groupby("grp"):
    print(f"{g:14s} {len(d):6d} " + " ".join(f"{np.nanmean(d[f'maf_{p}']):8.4f}" for p in POPS))

print("\n=== which superpopulation carries each SV's HIGHEST MAF? (composition, %) ===")
arg = vt[[f"maf_{p}" for p in POPS]].values
best = np.array(POPS)[np.nanargmax(np.nan_to_num(arg, nan=-1), axis=1)]
vt["top_pop"] = best
ct = pd.crosstab(vt.grp, vt.top_pop, normalize="index") * 100
print(ct.round(1).to_string())

print("\n=== EUR-invisible variants: MAF_EUR ~ 0 but common (>=0.05) elsewhere ===")
eur0 = (vt.maf_EUR.fillna(0) < 0.01)
elsewhere = vt[[f"maf_{p}" for p in POPS if p != "EUR"]].max(axis=1) >= 0.05
vt["eur_invisible"] = eur0 & elsewhere
for g, d in vt.groupby("grp"):
    print(f"  {g:14s} {100*d.eur_invisible.mean():5.1f}%  ({int(d.eur_invisible.sum()):4d} / {len(d)})")
tab = pd.crosstab(vt.grp, vt.eur_invisible)
odds, p = stats.fisher_exact(tab.values) if tab.shape == (2, 2) else (np.nan, np.nan)
print(f"  Fisher exact: OR={odds:.2f}  p={p:.3g}   (OR>1 => excluded set enriched for EUR-invisible SVs)")

print("\n=== mirror check: AFR-invisible but common elsewhere (control direction) ===")
afr0 = (vt.maf_AFR.fillna(0) < 0.01)
elsew2 = vt[[f"maf_{p}" for p in POPS if p != "AFR"]].max(axis=1) >= 0.05
vt["afr_invisible"] = afr0 & elsew2
for g, d in vt.groupby("grp"):
    print(f"  {g:14s} {100*d.afr_invisible.mean():5.1f}%  ({int(d.afr_invisible.sum()):4d} / {len(d)})")
tab2 = pd.crosstab(vt.grp, vt.afr_invisible)
odds2, p2 = stats.fisher_exact(tab2.values) if tab2.shape == (2, 2) else (np.nan, np.nan)
print(f"  Fisher exact: OR={odds2:.2f}  p={p2:.3g}")

out = f"{WU}/output/substrate/admit_ancestry_bias.tsv"
vt[["sv_locus", "chrom", "pos", "svtype", "svlen", "af", "maf", "status", "grp", "top_pop",
    "eur_invisible", "afr_invisible"] + [f"maf_{p}" for p in POPS]].to_csv(
    out, sep="\t", index=False, float_format="%.6g")
print(f"\nwrote {out}")
