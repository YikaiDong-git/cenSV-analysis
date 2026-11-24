#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Is r = 0.9952 a finding or a tautology?"""
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
U = ROOT + "/analyses/cdr_definition"
WU = ROOT + "/analyses/cdr_qtl_full_panel"
LEAD = 44300221
KEY = ["sample", "hap_local", "chm13_chrom"]

a = pd.read_csv(U + "/output/compat/compat_dip_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
b = pd.read_csv(U + "/output/compat/compat_gao_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
m = a.merge(b, on=KEY, suffixes=("_d", "_g"))
m = m[m.cdr_pos_offset_d.notna() & m.cdr_pos_offset_g.notna()].copy()
m["delta"] = m.cdr_pos_offset_g - m.cdr_pos_offset_d
m["adelta"] = m.delta.abs()

print("haplotype-chromosomes with a position in both units: %d" % len(m))
print()
print("|Gao - dip| on the [-1, +1] position scale:")
for q in (0.5, 0.75, 0.9, 0.95, 0.99, 1.0):
    print("   %5.0f%% of rows differ by <= %.4f" % (100 * q, m.adelta.quantile(q)))
for thr in (0.01, 0.05, 0.10, 0.20, 0.50):
    print("   %6.2f%% of rows differ by more than %.2f"
          % (100.0 * (m.adelta > thr).mean(), thr))
print()
print("identical on %d rows (%.1f%%) -- these are the single-dip arrays, where merging is a no-op"
      % (int((m.adelta == 0).sum()), 100.0 * (m.adelta == 0).mean()))
one = m[m.n_cdrs_d == 1]
print("   of the %d rows with exactly one dip, %d (%.1f%%) are identical"
      % (len(one), int((one.adelta == 0).sum()),
         100.0 * (one.adelta == 0).mean() if len(one) else float("nan")))

print()
print("disagreement grows with the number of dips merged, as it must:")
print("%-12s %8s %12s %12s" % ("n dips", "rows", "median |d|", "p90 |d|"))
for lo, hi, lbl in [(1, 1, "1"), (2, 3, "2-3"), (4, 6, "4-6"), (7, 10, "7-10"), (11, 999, "11+")]:
    s = m[(m.n_cdrs_d >= lo) & (m.n_cdrs_d <= hi)]
    if not len(s):
        continue
    print("%-12s %8d %12.4f %12.4f" % (lbl, len(s), s.adelta.median(), s.adelta.quantile(0.9)))

vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
loc = vt[(vt.chrom == "chr8") & (vt.pos == LEAD)].sv_locus.iloc[0]
G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
sm = open(WU + "/geno/full_samples.txt").read().split()
k = {l: j for j, l in enumerate(vt.sv_locus)}[loc]
dose = pd.Series(G[k], index=sm, name="dose")

c8 = m[m.chm13_chrom == "chr8"].join(dose, on="sample")
c8 = c8[c8.dose.notna()]
print()
print("chr8 only, %d haplotypes with a lead genotype:" % len(c8))
print("%-10s %8s %14s %14s" % ("dosage", "n", "median |d|", "mean |d|"))
for d in sorted(c8.dose.unique()):
    s = c8[c8.dose == d]
    print("%-10.1f %8d %14.4f %14.4f" % (d, len(s), s.adelta.median(), s.adelta.mean()))
from scipy.stats import kruskal
grp = [c8[c8.dose == d].adelta for d in sorted(c8.dose.unique())]
h, p = kruskal(*grp)
print("Kruskal-Wallis on |Gao - dip| across dosage classes: H = %.2f, P = %.3g" % (h, p))
if p > 0.05:
    print("-> the disagreement between units is NOT concentrated in any dosage class, so the")
    print("   lead surviving is not an accident of where the two units happen to agree.")
else:
    print("-> the disagreement DOES differ by dosage class; the survival of the lead needs a")
    print("   closer look before it can be called unit-independent.")
