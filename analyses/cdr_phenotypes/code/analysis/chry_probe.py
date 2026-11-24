#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""What Extended Data Fig."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
UNIT = ROOT + "/analyses/cdr_phenotypes"
SHOWN = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
ALPHA = 0.05 / 23

dp = pd.read_csv(UNIT + "/output/phenotypes/per_dip_positions.tsv", sep="\t")
print("per_dip_positions.tsv: %d dip rows" % len(dp))
print("columns: %s" % ", ".join(dp.columns[:12]))

present = sorted(dp.chm13_chrom.unique(), key=lambda c: (len(c), c))
print("\nchromosomes present in the table: %s" % " ".join(present))
print("chrY present: %s" % ("chrY" in present))
if "chrY" not in present:
    raise SystemExit("chrY carries no dip rows at all; the exclusion is a data limit, not a "
                     "choice, and nothing below applies")

dp = dp.copy()
if "rel" not in dp.columns:
    dp["rel"] = (dp.pos_splice + 1.0) / 2.0
if "len_mb" not in dp.columns:
    dp["len_mb"] = dp.active_len_sum / 1e6

hap = (dp.groupby(["sample", "hap_local", "chm13_chrom"])
         .agg(rel=("rel", "median"), len_mb=("len_mb", "first"))
         .reset_index())

print("\n%-6s %6s %6s %9s %9s %9s %8s %10s %s"
      % ("chrom", "haps", "dips", "len med", "len min", "len max", "rho", "P", ""))
print("-" * 82)
rows = []
for c in SHOWN + ["chrY"]:
    h = hap[hap.chm13_chrom == c]
    d = dp[dp.chm13_chrom == c]
    if len(h) < 3:
        print("%-6s %6d %6d   (too few haplotypes to correlate)" % (c, len(h), len(d)))
        rows.append((c, len(h), len(d), np.nan, np.nan))
        continue
    rho, p = spearmanr(h.len_mb, h.rel)
    star = "*" if p < ALPHA else ""
    mark = "   <-- NOT DRAWN" if c == "chrY" else ""
    print("%-6s %6d %6d %9.2f %9.2f %9.2f %+8.3f %10.2e %s%s"
          % (c, len(h), len(d), h.len_mb.median(), h.len_mb.min(), h.len_mb.max(),
             rho, p, star, mark))
    rows.append((c, len(h), len(d), rho, p))

y = hap[hap.chm13_chrom == "chrY"]
drawn = hap[hap.chm13_chrom.isin(SHOWN)]
print("\n--- chrY in context ---")
print("haplotypes: %d, against a median of %d over the 23 drawn (min %d on %s)"
      % (len(y),
         int(np.median([n for _c, n, _d, _r, _p in rows[:-1]])),
         min(n for _c, n, _d, _r, _p in rows[:-1]),
         min(((n, c) for c, n, _d, _r, _p in rows[:-1]))[1]))
print("dips per haplotype: chrY %.2f, drawn chromosomes %.2f"
      % (len(dp[dp.chm13_chrom == "chrY"]) / max(len(y), 1),
         len(dp[dp.chm13_chrom.isin(SHOWN)]) / len(drawn)))
print("active array length, chrY: median %.3f Mb, range %.3f-%.3f Mb"
      % (y.len_mb.median(), y.len_mb.min(), y.len_mb.max()))
print("active array length, the 23 drawn: median %.3f Mb, range %.3f-%.3f Mb"
      % (drawn.len_mb.median(), drawn.len_mb.min(), drawn.len_mb.max()))
print("smallest drawn array (chr21 in the caption): %.3f Mb median"
      % hap[hap.chm13_chrom == "chr21"].len_mb.median())
print("distinct samples contributing chrY: %d" % y["sample"].nunique())
print("distinct samples contributing chrX: %d"
      % hap[hap.chm13_chrom == "chrX"]["sample"].nunique())
print("distinct samples in the whole table: %d" % hap["sample"].nunique())

print("\nrelative dip position on chrY: median %.3f, IQR %.3f-%.3f"
      % (y.rel.median(), y.rel.quantile(.25), y.rel.quantile(.75)))
print("relative dip position, drawn: median %.3f" % drawn.rel.median())
