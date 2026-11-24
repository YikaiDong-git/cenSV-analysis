#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Is the number of dips per array roughly constant, independent of array length?"""
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
UNIT = ROOT + "/analyses/cdr_phenotypes"
ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

ph = pd.read_csv(UNIT + "/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv", sep="\t")
ph = ph[ph.chm13_chrom.isin(ORDER)].copy()
ph["mb"] = ph.active_len_sum / 1e6
d = ph[["mb", "n_cdrs", "chm13_chrom", "cdrs_per_mb_sum", "cdr_fraction_sum",
        "total_cdr_bp"]].dropna(subset=["mb", "n_cdrs"])

print("haplotype-chromosomes: %d" % len(d))
print("")
print("=" * 80)
print("ACROSS ALL HAPLOTYPES")
print("=" * 80)
for col, name in (("n_cdrs", "dip COUNT"), ("cdrs_per_mb_sum", "dips per Mb"),
                  ("total_cdr_bp", "total dip bp"), ("cdr_fraction_sum", "dip fraction")):
    if col not in d.columns:
        continue
    s = d[[col, "mb"]].dropna()
    r, p = stats.spearmanr(s[col], s.mb)
    print("  rho(%-14s, array length) = %+.3f   P = %.2e   n = %d" % (name, r, p, len(s)))

print("")
print("=" * 80)
print("PER-CHROMOSOME MEDIANS  (is count flat while density tracks 1/L?)")
print("=" * 80)
g = d.groupby("chm13_chrom").agg(mb=("mb", "median"), n=("n_cdrs", "median"),
                                 dens=("cdrs_per_mb_sum", "median")).reindex(ORDER)
print("%-7s %9s %9s %10s %12s" % ("chrom", "med Mb", "med dips", "dips/Mb", "pred 1/L*k"))
k = float(np.median(g.n))
for c, r in g.iterrows():
    print("%-7s %9.2f %9.1f %10.2f %12.2f" % (c, r.mb, r.n, r.dens, k / r.mb))
print("")
print("  median dip count across chromosomes: %.1f   (range %.1f - %.1f)"
      % (k, g.n.min(), g.n.max()))
rc, pc = stats.spearmanr(g.mb, g.n)
rd, pd_ = stats.spearmanr(g.mb, g.dens)
print("  rho(median length, median dip COUNT)   = %+.3f  P = %.3f" % (rc, pc))
print("  rho(median length, median dips per Mb) = %+.3f  P = %.3g" % (rd, pd_))

pred = k / g.mb
r2 = 1 - np.sum((g.dens - pred) ** 2) / np.sum((g.dens - g.dens.mean()) ** 2)
print("  a constant-count model (density = %.1f / L) explains R2 = %.3f of the "
      "between-chromosome density variation" % (k, r2))

print("")
print("=" * 80)
print("READING")
if abs(rc) < 0.35 and rd < -0.7:
    print("  Dip COUNT is essentially flat across chromosomes while dips-per-Mb tracks 1/L.")
    print("  Dip density is therefore close to a deterministic function of array length, and")
    print("  carries little information the length does not already carry.")
elif rc > 0.6:
    print("  Dip count SCALES with length -- density is the informative quantity and the")
    print("  hyperbolae in N5 are a plotting artefact, not a finding.")
else:
    print("  Intermediate: count rises with length but sub-proportionally. Report both.")
