#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""What actually enters the regression, and how similar is it between units?"""
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition/output")
PH = ["cdr_pos_offset", "cdr_pos_offset_deepest", "cdr_fraction", "cdrs_per_mb"]


def rint(a):
    r = stats.rankdata(a)
    return stats.norm.ppf((r - 0.5) / len(r))


arms = {}
for arm in ("dip", "gao"):
    a = pd.read_csv("%s/compat/compat_%s_splice/per_hap_chrom_phenotypes.tsv" % (U, arm),
                    sep="\t")
    b = pd.read_csv("%s/compat/compat_%s_splice/per_hap_chrom_phenotypes_deepest.tsv" % (U, arm),
                    sep="\t")
    a["cdr_pos_offset_deepest"] = b["cdr_pos_offset_deepest"]
    arms[arm] = a

print("%-24s %8s %12s %12s %12s %12s" %
      ("phenotype", "chroms", "raw Pearson", "raw Spearman", "RINT Pearson", "min RINT r"))
print("-" * 84)
for ph in PH:
    rp, rs, ri, per = [], [], [], []
    for chrom in sorted(set(arms["dip"].chm13_chrom)):
        v = {}
        for arm in ("dip", "gao"):
            s = arms[arm]
            s = s[s.chm13_chrom == chrom]
            v[arm] = s.groupby("sample")[ph].mean()
        j = pd.concat([v["dip"].rename("d"), v["gao"].rename("g")], axis=1).dropna()
        if len(j) < 20:
            continue
        rp.append(j.d.corr(j.g))
        rs.append(j.d.corr(j.g, method="spearman"))
        x, y = rint(j.d.to_numpy()), rint(j.g.to_numpy())
        r = float(np.corrcoef(x, y)[0, 1])
        ri.append(r)
        per.append((chrom, len(j), r))
    if not ri:
        continue
    print("%-24s %8d %12.4f %12.4f %12.4f %12.4f"
          % (ph, len(ri), float(np.median(rp)), float(np.median(rs)),
             float(np.median(ri)), min(ri)))

print()
print("Per chromosome, cdr_pos_offset after averaging haplotypes and RINT -- the exact vector the")
print("model regresses on:")
ph = "cdr_pos_offset"
rows = []
for chrom in sorted(set(arms["dip"].chm13_chrom)):
    v = {}
    for arm in ("dip", "gao"):
        s = arms[arm]
        s = s[s.chm13_chrom == chrom]
        v[arm] = s.groupby("sample")[ph].mean()
    j = pd.concat([v["dip"].rename("d"), v["gao"].rename("g")], axis=1).dropna()
    if len(j) < 20:
        continue
    x, y = rint(j.d.to_numpy()), rint(j.g.to_numpy())
    disc = int((stats.rankdata(j.d) != stats.rankdata(j.g)).sum())
    rows.append((chrom, len(j), float(np.corrcoef(x, y)[0, 1]), disc))
rows.sort(key=lambda t: t[2])
print("%-8s %8s %14s %20s" % ("chrom", "n", "RINT r", "individuals reranked"))
for c, n, r, d in rows:
    print("%-8s %8d %14.4f %20d" % (c, n, r, d))

print()
print("Reading: if RINT r is far above raw r, the units disagree on MAGNITUDE but rank individuals")
print("the same way, and the association is unit-free because the model only uses ranks. If the")
print("two are close, the units genuinely measure the same quantity. Both are legitimate reasons")
print("for the result to hold, but they are different claims and only one of them is true here.")
