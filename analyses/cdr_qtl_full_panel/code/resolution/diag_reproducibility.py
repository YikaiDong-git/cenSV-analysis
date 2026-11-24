#!/usr/bin/env python
"""Which of the 23 centromere trees can be rebuilt today from the VCF on disk."""
import importlib.util
import os
import re

import numpy as np
import pandas as pd
from skbio.tree import TreeNode
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_"
       "qgap200_rgap200")
DEC = os.path.join(ROOT, "sv_extraction/decompose", SUF)
VCF = os.path.join(ROOT, "variant_calling/svim_asm", "cohort_svim_asm_%s.vcf.gz" % SUF)
LOGS = os.path.join(ROOT, "logs/decompose")
HC = os.path.join(ROOT, "workflow/scripts/haplogroup_clustering.py")
WU = os.path.join(ROOT, "analyses/cdr_qtl_full_panel")
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

spec = importlib.util.spec_from_file_location("hc", HC)
hc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hc)

rows = []
for C in CHROMS:
    logged = np.nan
    for f in sorted(os.listdir(LOGS)):
        if f.startswith("dec_%s_" % C):
            for line in open(os.path.join(LOGS, f)):
                m = re.search(r"Dedup: \d+ -> (\d+) SVs", line)
                if m:
                    logged = int(m.group(1))
    matrix, variants, samples = hc.parse_vcf_region(VCF, C, 2, 1.0, None)
    hg = pd.read_csv("%s/%s/%s_haplogroups.tsv" % (DEC, C, C), sep="\t")
    TAB = hg["sample"].astype(str).tolist()
    at = {s: i for i, s in enumerate(samples)}
    order = TAB + [s for s in samples if s not in set(TAB)]
    matrix = matrix[[at[s] for s in order]]
    res = hc.perform_clustering(matrix, order, (2, 30), 5, "cosine_nj")
    bv = hc.extract_branch_vectors(res["nj_tree"], res["samples"], min_tips=5)
    mine = {frozenset(s for s, v in zip(res["samples"], b["vector"]) if v) for b in bv}

    saved = TreeNode.read("%s/%s/%s_nj.nwk" % (DEC, C, C)).root_at_midpoint()
    S = set(res["samples"])
    ship = set()
    for n in saved.traverse(include_self=False):
        if n.is_tip():
            continue
        u = n.subset() & S
        if 5 <= len(u) <= len(S) - 5:
            ship.add(frozenset(u))

    var = pd.read_csv("%s/%s/%s_sv_variants.tsv" % (DEC, C, C), sep="\t", low_memory=False,
                      usecols=["original_sv_id"])
    shipped_ids = set(var.original_sv_id.astype(str))
    mine_ids = set(v["id"] for v in variants)
    lb = sum(float(n.length or 0) for n in res["nj_tree"].traverse(include_self=False))
    ls = sum(float(n.length or 0) for n in TreeNode.read(
        "%s/%s/%s_nj.nwk" % (DEC, C, C)).traverse(include_self=False))
    rows.append(dict(chrom=C, svs_now=len(variants), svs_logged=logged,
                     haplotypes=len(res["samples"]), table_haplotypes=len(TAB),
                     ids_missing=len(shipped_ids - mine_ids),
                     ids_extra=len(mine_ids - shipped_ids),
                     clades_now=len(mine), clades_shipped=len(ship),
                     clades_shared=len(mine & ship),
                     length_now=lb, length_shipped=ls,
                     reproduces=bool(mine == ship)))
    r = rows[-1]
    print("  %-6s SVs %6d/%-6s  clades %3d/%3d shared %3d  %s"
          % (C, r["svs_now"], r["svs_logged"], r["clades_now"], r["clades_shipped"],
             r["clades_shared"], "OK" if r["reproduces"] else "DIFFERS"))

out = pd.DataFrame(rows)
out.to_csv(WU + "/output/resolution/tree_reproducibility.tsv", sep="\t", index=False)
print("\n%d of %d chromosomes rebuild to the same tree" % (int(out.reproduces.sum()), len(out)))
bad = out[~out.reproduces]
if len(bad):
    print("\nthe ones that differ:")
    print(bad[["chrom", "svs_now", "svs_logged", "ids_missing", "ids_extra", "clades_now",
               "clades_shipped", "clades_shared", "length_now", "length_shipped"]].to_string(
        index=False))
print("\nwrote %s/output/resolution/tree_reproducibility.tsv" % WU)
