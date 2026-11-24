#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The 1e7 escalation, compared arm against arm."""
import os
import glob
import sys
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition")
S = sys.argv[1] if len(sys.argv) > 1 else U + "/output/scan/scan_gw_1e7"
ALPHA = 0.05 / (11129 * 8)

rows, missing = [], []
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
    m["phenotype"] = pheno
    rows.append(m)

if missing:
    print("MISSING gao counterpart for %d files: %s" % (len(missing), missing[:8]))
if not rows:
    raise SystemExit("no escalated cells found under %s" % S)
M = pd.concat(rows, ignore_index=True)
gcol = "geno_group_dip" if "geno_group_dip" in M.columns else None
print("escalated rows: %d over %d cells (%d chromosomes)"
      % (len(M), M.groupby(["chrom", "phenotype"]).ngroups, M.chrom.nunique()))
print("alpha = %.3g" % ALPHA)

sd = M.p_emp_perm_dip < ALPHA
sg = M.p_emp_perm_gao < ALPHA


def ngrp(mask):
    if gcol:
        return M.loc[mask].groupby(["chrom", "phenotype", gcol]).ngroups
    return int(mask.sum())


print()
print("study-wide significant:")
print("   dip arm : %5d rows / %4d genotype groups" % (int(sd.sum()), ngrp(sd)))
print("   gao arm : %5d rows / %4d genotype groups" % (int(sg.sum()), ngrp(sg)))
print("   both    : %5d rows / %4d genotype groups" % (int((sd & sg).sum()), ngrp(sd & sg)))
print("   dip only: %5d rows / %4d genotype groups" % (int((sd & ~sg).sum()), ngrp(sd & ~sg)))
print("   gao only: %5d rows / %4d genotype groups" % (int((~sd & sg).sum()), ngrp(~sd & sg)))

print()
print("by (chromosome, phenotype):")
print("%-8s %-24s %8s %8s %8s %10s %10s" %
      ("chrom", "phenotype", "n_dip", "n_gao", "n_both", "max|t|dip", "max|t|gao"))
print("-" * 82)
for (c, ph), g in M.groupby(["chrom", "phenotype"]):
    d = g.p_emp_perm_dip < ALPHA
    o = g.p_emp_perm_gao < ALPHA
    if not (d.any() or o.any()):
        continue
    flag = "" if d.any() == o.any() else "   <-- UNIT-DEPENDENT"
    print("%-8s %-24s %8d %8d %8d %10.3f %10.3f%s"
          % (c, ph, int(d.sum()), int(o.sum()), int((d & o).sum()),
             g.abs_t_obs_dip.max(), g.abs_t_obs_gao.max(), flag))

print()
print("loci significant in ONE arm only (the symmetric difference):")
X = M[sd ^ sg].copy()
if X.empty:
    print("   none -- the two units return the same study-wide set")
else:
    X["arm"] = ["dip only" if r else "gao only" for r in (X.p_emp_perm_dip < ALPHA)]
    cols = ["chrom", "pos", "phenotype", "sv_locus", "arm",
            "t_obs_dip", "p_emp_perm_dip", "t_obs_gao", "p_emp_perm_gao"]
    if gcol:
        cols.insert(4, gcol)
    X = X.sort_values(["phenotype", "chrom", "pos"])
    print(X[cols].to_string(index=False, max_rows=60, float_format=lambda v: "%.4g" % v))

os.makedirs(U + "/output/gw", exist_ok=True)
M.to_csv(U + "/output/gw/escalated_1e7_both_arms.tsv", sep="\t", index=False,
         float_format="%.6g")
print()
print("[written] %s/output/gw/escalated_1e7_both_arms.tsv" % U)
print()
print("NOTE: significance here is the RAW scan. The manuscript's robustness screens "
      "(array-length\nstability, endpoint conditioning, ancestry stratification) are a "
      "separate subtractive pass and\nare not applied above.")
