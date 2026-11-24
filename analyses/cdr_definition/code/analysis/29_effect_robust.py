#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Is the raw-scale beta behind the kb effect size carried by a few extreme haplotypes?"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
U = ROOT / "analyses/cdr_definition"
WU = ROOT / "analyses/cdr_qtl_full_panel"
CANON = ROOT / "analyses/cdr_qtl"
sys.path.insert(0, str(CANON / "code/scan"))
from b6_cdr_pos_qtl_perchrom import read_gcta_grm_bin  # noqa: E402

CHROM, PHENO, LEAD_POS = "chr8", "cdr_pos_offset", 44300221
GRM = str(ROOT / "subproject_methylation/results/greml/pericen_1m_all")

SV_DIR = WU / "geno/perchrom" / CHROM
_npz = np.load(SV_DIR / ("%s_sv_genotypes.npz" % CHROM))
sv_dosages = _npz["arr_0"] if "arr_0" in _npz else _npz[list(_npz.keys())[0]]
sv_samples = open(SV_DIR / ("%s_sv_samples.txt" % CHROM)).read().split()
sv_variants = pd.read_csv(SV_DIR / ("%s_sv_variants.tsv" % CHROM), sep="\t")
if sv_dosages.shape[0] != len(sv_samples):
    sv_dosages = sv_dosages.T

pc = pd.read_csv(U / "output/compat/compat_dip_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
pc_chrom = pc[pc["chm13_chrom"] == CHROM]
pheno_ind = pc_chrom.groupby("sample")[PHENO].mean().reset_index()
bc_map = pc_chrom.groupby("sample")["basecaller"].first().to_dict()

K, grm_samples = read_gcta_grm_bin(GRM)
ev, evec = np.linalg.eigh(K)
PCs = evec[:, np.argsort(ev)[::-1][:10]]
pc_df = pd.DataFrame(PCs, columns=["PC%d" % (i + 1) for i in range(10)])
pc_df["sample"] = grm_samples

keep = sorted(set(bc_map) & set(pheno_ind["sample"]) & set(grm_samples) & set(sv_samples))
pmap = pheno_ind.set_index("sample")[PHENO].to_dict()
cmap = pc_df.set_index("sample")[["PC%d" % (i + 1) for i in range(10)]].to_dict("index")
sidx = {s: i for i, s in enumerate(sv_samples)}

y = np.array([pmap[s] for s in keep], float)
PC_mat = np.array([[cmap[s]["PC%d" % (i + 1)] for i in range(10)] for s in keep], float)
bc = np.array([1.0 if bc_map.get(s) == "R1041" else 0.0 for s in keep], float)
g = sv_dosages[[sidx[s] for s in keep], :].astype(float)[
    :, int(np.where(sv_variants["pos"].values == LEAD_POS)[0][0])]

m = ~np.isnan(y)
y, PC_mat, bc, g = y[m], PC_mat[m], bc[m], g[m]
keep = [k for i, k in enumerate(keep) if m[i]]
n = len(y)
X = np.column_stack([np.ones(n), PC_mat, bc]) if len(set(bc)) > 1 else np.column_stack([np.ones(n), PC_mat])


def beta_of(yv, gv, Xv):
    ry = yv - Xv @ np.linalg.lstsq(Xv, yv, rcond=None)[0]
    rg = gv - Xv @ np.linalg.lstsq(Xv, gv, rcond=None)[0]
    return float(np.dot(rg, ry) / np.dot(rg, rg)), ry, rg


B0, ry, rg = beta_of(y, g, X)
print("OLS raw-scale beta = %+.5f SPLICE/allele" % B0)

loo = np.array([beta_of(np.delete(y, i), np.delete(g, i), np.delete(X, i, 0))[0] for i in range(n)])
w = int(np.argmax(np.abs(loo - B0)))
print("\n1. LEAVE-ONE-OUT   min %+.5f  max %+.5f  worst shift %+.5f (%.1f%%) on sample %s"
      % (loo.min(), loo.max(), loo[w] - B0, 100 * abs(loo[w] - B0) / abs(B0), keep[w]))

lo, hi = np.percentile(ry, [2.5, 97.5])
t_mask = (ry >= lo) & (ry <= hi)
Bt = beta_of(y[t_mask], g[t_mask], X[t_mask])[0]
print("2. TRIMMED 5%%      beta %+.5f  (%.1f%% from full)  n %d -> %d"
      % (Bt, 100 * abs(Bt - B0) / abs(B0), n, t_mask.sum()))

d = rg.copy()
b, c = B0, 1.345 * np.median(np.abs(ry - B0 * rg)) / 0.6745
for _ in range(100):
    r = ry - b * d
    wgt = np.clip(c / np.maximum(np.abs(r), 1e-12), None, 1.0)
    bn = float(np.sum(wgt * d * ry) / np.sum(wgt * d * d))
    if abs(bn - b) < 1e-12:
        b = bn
        break
    b = bn
print("3. HUBER ROBUST    beta %+.5f  (%.1f%% from OLS)" % (b, 100 * abs(b - B0) / abs(B0)))

print("\nPER-HAPLOTYPE INTERPRETATION")
print("  y is the mean of an individual's two haplotype pos_splice values;")
print("  dosage counts alleles 0/1/2. A per-ALLELE slope on a two-haplotype MEAN is already")
print("  the per-haplotype effect scaled by 1/2 only if each allele acts on its own haplotype.")
gr = np.round(g)
for v in sorted(set(gr)):
    s = gr == v
    print("  dosage %d  n=%3d  mean %+.4f  sd %.4f" % (v, s.sum(), y[s].mean(), y[s].std()))
print("  0->1 step %+.4f   1->2 step %+.4f   (linear model assumes these are equal)"
      % (y[gr == 1].mean() - y[gr == 0].mean(), y[gr == 2].mean() - y[gr == 1].mean()))

AL = pc_chrom.groupby("sample")["active_len"].mean()
half = np.array([AL.get(s, np.nan) for s in keep], float) / 2.0
print("\nKB AT EACH ESTIMATE (median chr8 array, half = %.0f kb)" % (np.median(half) / 1e3))
for nm, bb in (("OLS", B0), ("trimmed", Bt), ("Huber", b)):
    print("  %-8s %+.1f kb per allele" % (nm, bb * np.median(half) / 1e3))
