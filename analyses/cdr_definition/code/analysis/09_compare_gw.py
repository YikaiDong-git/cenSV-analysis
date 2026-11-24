#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Genome-wide comparison of the two CDR-unit arms, and the escalation list."""
import os
import sys
import glob
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition")
S = sys.argv[1] if len(sys.argv) > 1 else U + "/output/scan/scan_gw_1e5"
NPERM = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000
ALPHA = 0.05 / (11129 * 8)
FLOOR = 1.5 / (NPERM + 1)

rows, cells, missing = [], [], []
for p in sorted(glob.glob(S + "/dip/perchrom/*_perm.tsv")):
    base = os.path.basename(p)
    q = S + "/gao/perchrom/" + base
    if not os.path.exists(q):
        missing.append(base)
        continue
    chrom = base.split("_")[0]
    pheno = base[len(chrom) + 1:-len("_perm.tsv")]
    a = pd.read_csv(p, sep="\t")
    b = pd.read_csv(q, sep="\t")
    m = a.merge(b, on=["chrom", "pos", "sv_locus"], suffixes=("_dip", "_gao"))
    if m.empty:
        continue
    ok = m.t_obs_dip.notna() & m.t_obs_gao.notna()
    sp = m.loc[ok, "t_obs_dip"].corr(m.loc[ok, "t_obs_gao"], method="spearman")
    k = min(20, len(m))
    ov = len(set(m.nlargest(k, "abs_t_obs_dip").sv_locus)
             & set(m.nlargest(k, "abs_t_obs_gao").sv_locus))
    fl_d = int((m.p_emp_perm_dip <= FLOOR).sum())
    fl_g = int((m.p_emp_perm_gao <= FLOOR).sum())
    rows.append(dict(chrom=chrom, phenotype=pheno, n_loci=len(m),
                     max_abs_t_dip=m.abs_t_obs_dip.max(), max_abs_t_gao=m.abs_t_obs_gao.max(),
                     spearman_t=sp, top20_overlap=ov,
                     n_floor_dip=fl_d, n_floor_gao=fl_g))
    if fl_d or fl_g:
        cells.append((chrom, pheno, fl_d, fl_g))

if missing:
    print("MISSING gao counterpart for %d files: %s" % (len(missing), missing[:6]))
d = pd.DataFrame(rows).sort_values(["phenotype", "chrom"])
os.makedirs(U + "/output/gw", exist_ok=True)
d.to_csv(U + "/output/gw/gw_arm_comparison.tsv", sep="\t", index=False, float_format="%.6g")

print("cells compared: %d   (%d chromosomes x %d phenotypes)"
      % (len(d), d.chrom.nunique(), d.phenotype.nunique()))
print()
print("Per phenotype, across chromosomes:")
print("%-24s %8s %10s %10s %12s %10s" %
      ("phenotype", "cells", "min rho", "median rho", "max|t| dip", "max|t| gao"))
print("-" * 80)
for ph, g in d.groupby("phenotype"):
    print("%-24s %8d %10.4f %10.4f %12.3f %10.3f"
          % (ph, len(g), g.spearman_t.min(), g.spearman_t.median(),
             g.max_abs_t_dip.max(), g.max_abs_t_gao.max()))

print()
print("Cells where the two arms disagree most on ranking (lowest Spearman rho):")
print("%-8s %-24s %8s %10s %10s %10s" %
      ("chrom", "phenotype", "n_loci", "rho", "max|t|dip", "max|t|gao"))
print("-" * 76)
for _, r in d.nsmallest(10, "spearman_t").iterrows():
    print("%-8s %-24s %8d %10.4f %10.3f %10.3f"
          % (r.chrom, r.phenotype, r.n_loci, r.spearman_t,
             r.max_abs_t_dip, r.max_abs_t_gao))

print()
print("Cells holding a locus at the %g permutation floor (p <= %.3g):" % (NPERM, FLOOR))
print("%-8s %-24s %10s %10s" % ("chrom", "phenotype", "n_floor_dip", "n_floor_gao"))
print("-" * 56)
for c, ph, fd, fg in sorted(cells):
    flag = "" if (fd > 0) == (fg > 0) else "   <-- ARMS DISAGREE"
    print("%-8s %-24s %10d %10d%s" % (c, ph, fd, fg, flag))

esc = U + "/output/tasks/escalate_cells.txt"
with open(esc, "w") as f:
    for c, ph, _, _ in sorted(cells):
        f.write("%s\t%s\n" % (c, ph))
print()
print("[written] %s : %d cells to escalate to 1e7 (union over both arms)" % (esc, len(cells)))
print("[written] %s/output/gw/gw_arm_comparison.tsv" % U)
print()
print("alpha (study-wide, manuscript denominator 11,129 x 8) = %.3g" % ALPHA)
print("NOTHING is significant at %g permutations; the floor here is %.3g." % (NPERM, FLOOR))
