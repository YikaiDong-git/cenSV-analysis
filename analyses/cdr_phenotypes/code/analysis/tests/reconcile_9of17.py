#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""STATUS records the haplotype-level tertile-KS count as 9 of 17."""
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
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

print("%8s  %9s  %5s  %s" % ("min n_hap", "testable", "sig", "significant chromosomes"))
print("-" * 78)
for thr in (0, 60, 100, 120, 140, 150, 155, 158, 160, 200, 240):
    testable, hit = 0, []
    for c in ORDER:
        h = hap[hap.chm13_chrom == c]
        if len(h) < thr:
            continue
        q1, q2 = h.len_mb.quantile([1 / 3, 2 / 3])
        sh, lo = h[h.len_mb <= q1].rel, h[h.len_mb >= q2].rel
        if len(sh) < 20 or len(lo) < 20:
            continue
        testable += 1
        if ks_2samp(sh, lo)[1] < ALPHA:
            hit.append(c.replace("chr", ""))
    mark = ""
    if (len(hit), testable) == (9, 17):
        mark = "   <-- reproduces the 9 of 17 in STATUS"
    if (len(hit), testable) == (12, 23):
        mark = "   <-- the permissive re-run"
    print("%8d  %9d  %5d  %s%s" % (thr, testable, len(hit), ", ".join(hit), mark))
print("-" * 78)
print("Both counts are the same test. They differ only in how many haplotypes a chromosome must")
print("have to be tested at all -- an arbitrary cut that changes the headline. Spearman needs no")
print("such cut and is what the accepted supplementary figure marks.")
