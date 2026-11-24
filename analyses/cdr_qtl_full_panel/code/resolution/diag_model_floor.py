#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Which condition decided whether a monophyletic component reached an imputation model."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from resfig import all_components, metrics, resolved_svs, thousands  # noqa: E402

ORD = ["non_decomposed", "residual", "original_decomposed", "monophyletic"]
ac = all_components()
met = metrics()
svs = resolved_svs()

print("=" * 72)
print("does 'each component retaining at least ten carriers' describe the modelled set?")
print("=" * 72)
mod, un = ac[ac.modelled], ac[~ac.modelled]
print("  components %d   modelled %d   not modelled %d" % (len(ac), len(mod), len(un)))
print("  modelled     carriers  min %d  max %d" % (mod.n_carriers.min(), mod.n_carriers.max()))
print("  not modelled carriers  min %d  max %d" % (un.n_carriers.min(), un.n_carriers.max()))
for t in (5, 8, 10, 12, 15):
    print("    floor %2d: modelled below it %4d   unmodelled at or above it %4d"
          % (t, int((mod.n_carriers < t).sum()), int((un.n_carriers >= t).sum())))

both = ac[(ac.n_carriers >= 10) & (ac.n_carriers <= 25)]
print("\n  where the two rules disagree, by carrier count")
print(both.groupby([pd.cut(both.n_carriers, [9, 12, 15, 18, 21, 25]), "modelled"],
                   observed=True).size().to_string())
print("\n  unmodelled components above twelve carriers")
print(ac[(~ac.modelled) & (ac.n_carriers > 12)][
    ["chrom", "branch_id", "n_carriers", "branch_r2", "clade_tips", "iteration"]
].to_string(index=False))

print("\n  smallest carrier frequency among modelled records, by class")
for k in ORD:
    s = met[met.sv_category == k]
    print("    %-22s n %7s  min %.4f  median %.4f"
          % (k, thousands(len(s)), s.svfreq.min(), s.svfreq.median()))
print("  0.0266 is ten carrier haplotypes in 376, which is the floor expressed as a frequency")

print("\n" + "=" * 72)
print("the same, restricted to the 138 SVs with more than one origin")
print("=" * 72)
multi = svs[svs.n_origins >= 2][["chrom", "original_sv_id"]]
mc = ac.merge(multi, on=["chrom", "original_sv_id"], how="inner")
print("  components of those SVs %d   modelled %d" % (len(mc), int(mc.modelled.sum())))
first = mc[mc.origin == "first"]
later = mc[mc.origin == "later"]
print("  first origins  %3d, modelled %3d" % (len(first), int(first.modelled.sum())))
print("  later origins  %3d, modelled %3d" % (len(later), int(later.modelled.sum())))
print("  SVs with EVERY origin modelled                    %d of %d"
      % (int(mc.groupby(["chrom", "original_sv_id"]).modelled.all().sum()), len(multi)))
ge10 = mc[mc.n_carriers >= 10]
print("  SVs whose 10-or-more-carrier origins were modelled %d of %d"
      % (int(ge10.groupby(["chrom", "original_sv_id"]).modelled.all().sum()), len(multi)))

print("\n" + "=" * 72)
print("what chromosomes 15 and 18 contribute")
print("=" * 72)
tab = pd.crosstab(met.chrom, met.sv_category).reindex(columns=ORD, fill_value=0)
tot, two = tab.sum(), tab.loc[["chr15", "chr18"]].sum()
for k in ORD:
    print("  %-22s all %7s   both %5d   %5.1f%%"
          % (k, thousands(tot[k]), two[k], 100.0 * two[k] / tot[k]))
keep = met[~met.chrom.isin(["chr15", "chr18"])]
for k in ORD:
    print("    %-22s published mean %.3f   without the two %.3f"
          % (k, met[met.sv_category == k].cv_r2.mean(), keep[keep.sv_category == k].cv_r2.mean()))
g = svs.assign(m=svs.n_origins >= 2).groupby("chrom").agg(svs=("m", "size"), multi=("m", "sum"))
print("  resolved SVs %d of %d, multi-origin %d of %d"
      % (int(g.loc[["chr15", "chr18"], "svs"].sum()), int(g.svs.sum()),
         int(g.loc[["chr15", "chr18"], "multi"].sum()), int(g.multi.sum())))
print("  records to re-impute if they were regenerated: %s, times eight strategies %s fits"
      % (thousands(int(two.sum())), thousands(8 * int(two.sum()))))
