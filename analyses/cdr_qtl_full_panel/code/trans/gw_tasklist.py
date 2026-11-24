#!/usr/bin/env python
"""Write the (phenotype, target chromosome) task list for the genome-wide trans scan."""
import sys
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")
CANON = (CENSV_ROOT + "/"
         "analyses/cdr_qtl")
sys.path.insert(0, f"{WU}/code")
from phenotype_family import PRIMARY

PC = pd.read_csv(f"{CANON}/output/per_hap_chrom_phenotypes.tsv", sep="\t",
                 engine="python")
chroms = sorted(set(PC.chm13_chrom.dropna()))
lines = [f"{p} {c}" for p in PRIMARY for c in chroms]
with open(f"{WU}/output/trans/gw_tasks.txt", "w") as fh:
    fh.write("\n".join(lines) + "\n")
print(f"{len(PRIMARY)} phenotypes x {len(chroms)} chroms = {len(lines)} tasks")
print("chroms:", " ".join(chroms))
