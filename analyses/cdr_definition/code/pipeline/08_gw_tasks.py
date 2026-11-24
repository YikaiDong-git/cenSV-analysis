#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Task list for the genome-wide CDR-unit head-to-head."""
import os
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition")
PHENO = ["cdr_pos_offset", "cdr_pos_offset_deepest",
         "n_cdrs", "cdrs_per_mb", "cdr_fraction", "median_cdr_width"]

p = pd.read_csv(U + "/output/compat/compat_gao_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
GENO = (CENSV_ROOT + "/analyses/"
        "cdr_qtl_full_panel/geno/perchrom")
have_geno = set()
for f in os.listdir(GENO):
    for tok in f.replace(".", "_").split("_"):
        if tok.startswith("chr"):
            have_geno.add(tok)
            break

chroms = sorted(set(p.chm13_chrom) & have_geno,
                key=lambda c: (len(c), c))
missing = sorted(set(p.chm13_chrom) - have_geno)
print("chromosomes with both phenotypes and genotypes: %d  %s" % (len(chroms), chroms))
if missing:
    print("EXCLUDED (no genotype file): %s" % missing)

os.makedirs(U + "/output/tasks", exist_ok=True)
out = U + "/output/tasks/gw_tasks.txt"
n = 0
with open(out, "w") as f:
    for arm in ("dip", "gao"):
        for c in chroms:
            for ph in PHENO:
                f.write("%s\t%s\t%s\n" % (arm, c, ph))
                n += 1
print("[written] %s : %d tasks (%d chroms x %d phenotypes x 2 arms)"
      % (out, n, len(chroms), len(PHENO)))
assert n == len(chroms) * len(PHENO) * 2
