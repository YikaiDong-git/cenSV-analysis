#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The tertile-KS covariation test at the HAPLOTYPE level, beside the per-CDR version."""
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, spearmanr
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + ""
        "/analyses/cdr_phenotypes")
ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
ALPHA = 0.05 / len(ORDER)

dp = pd.read_csv(UNIT + "/output/phenotypes/per_dip_positions.tsv", sep="\t")
dp = dp[dp.chm13_chrom.isin(ORDER)].copy()
dp["rel"] = (dp.pos_splice + 1.0) / 2.0
dp["len_mb"] = dp.active_len_sum / 1e6

hap = (dp.groupby(["sample", "hap_local", "chm13_chrom"])
       .agg(rel=("rel", "mean"), len_mb=("len_mb", "first")).reset_index())

print("units: %d CDRs -> %d haplotype-chromosomes (%.2f CDRs each)"
      % (len(dp), len(hap), len(dp) / len(hap)))
print()
print("%-7s %6s %6s %6s   %10s %8s   %10s %8s"
      % ("chrom", "n_hap", "n_sh", "n_lo", "KS P", "KS sig", "rho", "rho sig"))
print("-" * 78)
ks_hit, sp_hit, testable = [], [], []
for c in ORDER:
    h = hap[hap.chm13_chrom == c]
    q1, q2 = h.len_mb.quantile([1 / 3, 2 / 3])
    sh = h[h.len_mb <= q1].rel
    lo = h[h.len_mb >= q2].rel
    rho, prho = spearmanr(h.len_mb, h.rel)
    if prho < ALPHA:
        sp_hit.append(c)
    if len(sh) >= 20 and len(lo) >= 20:
        testable.append(c)
        _, pks = ks_2samp(sh, lo)
        if pks < ALPHA:
            ks_hit.append(c)
        pks_s = "%10.2e" % pks
        sig = "*" if pks < ALPHA else ""
    else:
        pks_s, sig = "%10s" % "n too low", ""
    print("%-7s %6d %6d %6d   %s %8s   %+10.3f %8s"
          % (c, len(h), len(sh), len(lo), pks_s, sig, rho, "*" if prho < ALPHA else ""))
print("-" * 78)
print("tertile-KS at HAPLOTYPE level : %d of %d testable chromosomes  (%s)"
      % (len(ks_hit), len(testable), ", ".join(x.replace("chr", "") for x in ks_hit)))
print("Spearman at HAPLOTYPE level   : %d of %d chromosomes            (%s)"
      % (len(sp_hit), len(ORDER), ", ".join(x.replace("chr", "") for x in sp_hit)))
print()
print("for comparison, the per-CDR tertile-KS count in the current draft is 19 of 23.")
print("the accepted supplementary figure marks SPEARMAN significance, so quoting %d of %d in"
      % (len(sp_hit), len(ORDER)))
print("Results lets a reader count the asterisks in that figure and get the same number.")
