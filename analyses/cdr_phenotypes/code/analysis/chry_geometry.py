#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Is chromosome Y's rho = -0.77 a finding, or arithmetic?"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
UNIT = ROOT + "/analyses/cdr_phenotypes"
ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]
EDGE_WINDOW_MB = 0.5

dp = pd.read_csv(UNIT + "/output/phenotypes/per_dip_positions.tsv", sep="\t")
dp = dp[dp.chm13_chrom.isin(ORDER)].copy()
dp["rel"] = (dp.pos_splice + 1.0) / 2.0
dp["len_mb"] = dp.active_len_sum / 1e6
dp["abs_mb"] = dp.rel * dp.len_mb

hap = (dp.groupby(["sample", "hap_local", "chm13_chrom"])
       .agg(rel=("rel", "mean"), abs_mb=("abs_mb", "mean"), len_mb=("len_mb", "first"))
       .reset_index())

print("%-6s %5s %8s %8s %9s %9s %9s %9s  %s"
      % ("chrom", "haps", "len med", "CV len", "CV rel", "CV abs", "rho rel", "rho abs", "read"))
print("-" * 96)
rows = []
for c in ORDER:
    h = hap[hap.chm13_chrom == c]
    if len(h) < 5:
        continue
    cv_len = h.len_mb.std() / h.len_mb.mean()
    cv_rel = h.rel.std() / h.rel.mean()
    cv_abs = h.abs_mb.std() / h.abs_mb.mean()
    rho_rel, _p1 = spearmanr(h.len_mb, h.rel)
    rho_abs, p2 = spearmanr(h.len_mb, h.abs_mb)
    read = "geometric" if (cv_abs < cv_rel and abs(rho_abs) < 0.35) else \
           ("tracks length" if rho_abs > 0.5 else "mixed")
    rows.append((c, read))
    print("%-6s %5d %8.2f %8.2f %9.2f %9.2f %+9.2f %+9.2f  %s"
          % (c, len(h), h.len_mb.median(), cv_len, cv_rel, cv_abs, rho_rel, rho_abs, read))

print("\nHaplotypes whose ACTIVE ARRAY is shorter than CDR-Finder's %.1f Mb edge window"
      % EDGE_WINDOW_MB)
tot = 0
for c in ORDER:
    h = hap[hap.chm13_chrom == c]
    if not len(h):
        continue
    k = int((h.len_mb < EDGE_WINDOW_MB).sum())
    tot += k
    if k:
        print("  %-6s %4d of %4d  (%.0f%%)" % (c, k, len(h), 100.0 * k / len(h)))
print("  total %d of %d haplotype-chromosomes" % (tot, len(hap)))

y = hap[hap.chm13_chrom == "chrY"]
print("\nchrY in detail")
print("  absolute offset from array start: median %.3f Mb, IQR %.3f-%.3f, CV %.2f"
      % (y.abs_mb.median(), y.abs_mb.quantile(.25), y.abs_mb.quantile(.75),
         y.abs_mb.std() / y.abs_mb.mean()))
print("  relative position:                median %.3f, IQR %.3f-%.3f, CV %.2f"
      % (y.rel.median(), y.rel.quantile(.25), y.rel.quantile(.75), y.rel.std() / y.rel.mean()))
print("  array length:                     median %.3f Mb, IQR %.3f-%.3f, CV %.2f"
      % (y.len_mb.median(), y.len_mb.quantile(.25), y.len_mb.quantile(.75),
         y.len_mb.std() / y.len_mb.mean()))
r, p = spearmanr(y.len_mb, y.abs_mb)
print("  absolute offset against length:   rho = %+.3f, P = %.2g" % (r, p))
print("  dips per haplotype:               %.2f" % (len(dp[dp.chm13_chrom == "chrY"]) / len(y)))
print("  haplotypes with a single dip:     %d of %d"
      % (int((dp[dp.chm13_chrom == "chrY"].groupby(["sample", "hap_local"]).size() == 1).sum()),
         len(y)))
