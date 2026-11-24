#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Is the per-allele kb effect the per-HAPLOTYPE effect, or half of it?"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
U = ROOT / "analyses/cdr_definition"
WU = ROOT / "analyses/cdr_qtl_full_panel"
CHROM, PHENO, LEAD_POS = "chr8", "cdr_pos_offset", 44300221

pc = pd.read_csv(U / "output/compat/compat_dip_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
pcc = pc[pc["chm13_chrom"] == CHROM].copy()

print("PHENOTYPE TABLE STRUCTURE")
print("  rows for %s: %d over %d samples" % (CHROM, len(pcc), pcc["sample"].nunique()))
hapcols = [c for c in pcc.columns if "hap" in c.lower()]
print("  haplotype-identifying columns: %s" % hapcols)

SV_DIR = WU / "geno/perchrom" / CHROM
_npz = np.load(SV_DIR / ("%s_sv_genotypes.npz" % CHROM))
keys = list(_npz.keys())
arr = _npz["arr_0"] if "arr_0" in _npz else _npz[keys[0]]
sv_samples = open(SV_DIR / ("%s_sv_samples.txt" % CHROM)).read().split()
print("\nSV GENOTYPE SOURCE")
print("  npz keys: %s   matrix %s   n_sample_ids %d" % (keys, arr.shape, len(sv_samples)))
print("  sample ids unique: %d" % len(set(sv_samples)))
if arr.shape[0] != len(sv_samples):
    arr = arr.T
j = int(np.where(pd.read_csv(SV_DIR / ("%s_sv_variants.tsv" % CHROM), sep="\t")["pos"].values == LEAD_POS)[0][0])
lead = arr[:, j]
vals = np.unique(lead[~np.isnan(lead)])
print("  lead dosage distinct values: %s" % np.round(vals, 4)[:10])
per_hap_geno = len(sv_samples) != len(set(sv_samples))
print("  => genotypes are %s" % ("PER-HAPLOTYPE" if per_hap_geno else "PER-INDIVIDUAL (0/1/2)"))

two = pcc.groupby("sample").filter(lambda d: len(d) == 2)
g = two.groupby("sample")["active_len"]
L1, L2 = g.nth(0).values, g.nth(1).values
d = np.abs(L1 - L2)
print("\nPER-HAPLOTYPE ARRAY LENGTH (the kb denominator)")
print("  two-hap samples: %d   median |L1-L2| %.0f kb (%.1f%%)   max %.0f kb (%.1f%%)"
      % (len(d), np.median(d) / 1e3, 100 * np.median(d / ((L1 + L2) / 2)),
         d.max() / 1e3, 100 * (d / ((L1 + L2) / 2)).max()))

BETA = -0.26241
half_mean = np.median(pcc.groupby("sample")["active_len"].mean()) / 2.0
half_all = np.median(pcc["active_len"]) / 2.0
print("\nKB UNDER TWO DEFENSIBLE DENOMINATORS")
print("  median of per-individual MEAN length : half = %7.0f kb -> %+.1f kb per allele"
      % (half_mean / 1e3, BETA * half_mean / 1e3))
print("  median over ALL haplotypes           : half = %7.0f kb -> %+.1f kb per allele"
      % (half_all / 1e3, BETA * half_all / 1e3))
print("  difference between the two: %.1f kb (%.1f%%)"
      % (abs(BETA * half_mean / 1e3 - BETA * half_all / 1e3),
         100 * abs(half_mean - half_all) / half_mean))

print("\nCIS/TRANS ESTIMAND")
if per_hap_geno:
    print("  Per-haplotype genotypes EXIST -> the decisive haplotype-level refit is runnable.")
else:
    print("  Genotypes are PER-INDIVIDUAL only, so a haplotype-level refit is NOT runnable.")
    print("  The design therefore CANNOT distinguish cis from trans, and the estimand is fixed:")
    print("  beta is the per-allele effect on the individual's HAPLOTYPE-AVERAGED CDR position.")
    print("  Under a cis mechanism the per-haplotype displacement would be ~2x (%+.0f kb);"
          % (2 * BETA * half_mean / 1e3))
    print("  that value is a model-dependent implication, NOT an estimate this design supports.")
