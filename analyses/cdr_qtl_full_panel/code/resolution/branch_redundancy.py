#!/usr/bin/env python
"""How many distinct genotypes the 1,962 monophyletic components actually carry."""
import os

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_"
       "qgap200_rgap200")
DEC = os.path.join(ROOT, "sv_extraction/decompose", SUF)
WU = os.path.join(ROOT, "analyses/cdr_qtl_full_panel")
A = os.path.join(WU, "output/resolution")
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

met = pd.read_csv(A + "/metrics_all.tsv", sep="\t", low_memory=False)
mono = met[met.sv_category == "monophyletic"][
    ["chrom", "original_sv_id", "mapped_haplogroup", "cv_r2", "svlen", "svtype", "svpos"]]
modelled = set(zip(mono.chrom, mono.original_sv_id, mono.mapped_haplogroup.astype(str)))

rows, pairs = [], []
for C in CHROMS:
    hg = pd.read_csv("%s/%s/%s_haplogroups.tsv" % (DEC, C, C), sep="\t")
    haps = hg["sample"].astype(str).tolist()
    d = pd.read_csv("%s/%s/%s_decomposition.tsv" % (DEC, C, C), sep="\t")
    d = d[d.type == "monophyletic"].copy()
    d["carriers"] = d.carrier_indices.fillna("").apply(
        lambda s: frozenset(int(i) for i in s.split(",") if s))
    d["is_modelled"] = [(C, s, str(b)) in modelled for s, b in zip(d.sv_id, d.branch_id)]

    for r in d.itertuples():
        rows.append(dict(chrom=C, sv_id=r.sv_id, branch_id=str(r.branch_id),
                         n_carriers=len(r.carriers), modelled=r.is_modelled,
                         key=hash(r.carriers)))

    m = d[d.is_modelled]
    for b, sub in m.groupby(m.branch_id.astype(str)):
        cs = list(sub.carriers)
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                pairs.append(dict(chrom=C, branch_id=b, same_branch=True,
                                  jaccard=len(cs[i] & cs[j]) / float(len(cs[i] | cs[j]))))
    rng = np.random.default_rng(20260910)
    cs = list(m.carriers)
    bs = list(m.branch_id.astype(str))
    for _ in range(min(4000, len(cs) * 8)):
        if len(cs) < 2:
            break
        i, j = rng.integers(len(cs), size=2)
        if i == j or bs[i] == bs[j]:
            continue
        pairs.append(dict(chrom=C, branch_id="", same_branch=False,
                          jaccard=len(cs[i] & cs[j]) / float(len(cs[i] | cs[j]))))
    print("  %-6s %4d components, %4d modelled, %3d distinct branches used"
          % (C, len(d), int(d.is_modelled.sum()), m.branch_id.nunique()))

g = pd.DataFrame(rows)
pr = pd.DataFrame(pairs)
g.to_csv(A + "/component_redundancy.tsv", sep="\t", index=False)
pr.to_csv(A + "/component_pairs.tsv", sep="\t", index=False)

m = g[g.modelled]
assert len(m) == 1962, "%d modelled components, expected 1,962" % len(m)
print("\nmodelled components                      %5d" % len(m))
print("distinct (chromosome, branch) pairs used %5d" % m.groupby(
    ["chrom", "branch_id"]).ngroups)
print("distinct carrier sets among them         %5d" % m.groupby(["chrom", "key"]).ngroups)
dup = m.groupby(["chrom", "key"]).size()
print("carrier sets carried by more than one component: %d, covering %d components"
      % (int((dup > 1).sum()), int(dup[dup > 1].sum())))
print("largest group of components sharing one carrier set: %d" % int(dup.max()))
sb = pr[pr.same_branch]
xb = pr[~pr.same_branch]
print("\npairwise carrier overlap of modelled components")
print("  sharing a branch     n=%6d  median Jaccard %.3f  identical %.1f%%"
      % (len(sb), sb.jaccard.median(), 100 * (sb.jaccard == 1).mean()))
print("  on different branches n=%6d  median Jaccard %.3f  identical %.1f%%"
      % (len(xb), xb.jaccard.median(), 100 * (xb.jaccard == 1).mean()))
print("\nwrote %s/component_redundancy.tsv" % A)
