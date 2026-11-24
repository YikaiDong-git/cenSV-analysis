#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Which phenotypes remain testable under the Gao CDR unit?"""
import pandas as pd
import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition/output")
PH = ["cdr_pos_offset", "cdr_pos_offset_deepest", "cdr_pos_iqr", "cdr_pos_var",
      "n_cdrs", "cdrs_per_mb", "cdr_fraction", "median_cdr_width", "mean_5mC", "mean_5mC_covwt"]

arms = {}
for arm in ("dip", "gao"):
    a = pd.read_csv("%s/compat/compat_%s_splice/per_hap_chrom_phenotypes.tsv" % (U, arm),
                    sep="\t")
    b = pd.read_csv("%s/compat/compat_%s_splice/per_hap_chrom_phenotypes_deepest.tsv" % (U, arm),
                    sep="\t")
    a["cdr_pos_offset_deepest"] = b["cdr_pos_offset_deepest"]
    arms[arm] = a

print("%-24s %10s %10s   %10s %10s   %8s" %
      ("phenotype", "dip_sd", "gao_sd", "dip_%zero", "gao_%zero", "r(d,g)"))
print("-" * 82)
verdict = {}
for p in PH:
    d, g = arms["dip"][p], arms["gao"][p]
    ok = d.notna() & g.notna()
    zd = 100.0 * (d == 0).sum() / d.notna().sum()
    zg = 100.0 * (g == 0).sum() / g.notna().sum()
    r = d[ok].corr(g[ok])
    print("%-24s %10.4g %10.4g   %9.1f%% %9.1f%%   %8.4f"
          % (p, d.std(), g.std(), zd, zg, r))
    per_chrom_sd = arms["gao"].groupby("chm13_chrom")[p].std()
    frac_chrom_flat = float((per_chrom_sd.fillna(0) == 0).mean())
    verdict[p] = dict(gao_sd=g.std(), gao_pct_zero=zg, chrom_flat=frac_chrom_flat)

print()
print("Per-chromosome variability under the Gao unit (a scan is run WITHIN a chromosome, so a")
print("phenotype that is constant within chromosomes cannot associate with anything there):")
print("%-24s %22s %14s" % ("phenotype", "chroms with sd = 0", "verdict"))
print("-" * 64)
for p in PH:
    v = verdict[p]
    testable = v["gao_sd"] > 0 and v["chrom_flat"] < 0.5 and v["gao_pct_zero"] < 90
    print("%-24s %21.0f%% %14s" % (p, 100 * v["chrom_flat"],
                                   "TESTABLE" if testable else "DEGENERATE"))
