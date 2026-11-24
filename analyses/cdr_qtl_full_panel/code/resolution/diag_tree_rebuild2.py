#!/usr/bin/env python
"""Where the rebuilt tree and the shipped tree part company, and whether it reaches the result."""
import importlib.util
import os

import numpy as np
import pandas as pd
import skbio
from skbio.tree import TreeNode
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_"
       "qgap200_rgap200")
DEC = os.path.join(ROOT, "sv_extraction/decompose", SUF)
VCF = os.path.join(ROOT, "variant_calling/svim_asm", "cohort_svim_asm_%s.vcf.gz" % SUF)
HC = os.path.join(ROOT, "workflow/scripts/haplogroup_clustering.py")
CHROM = "chr9"

spec = importlib.util.spec_from_file_location("hc", HC)
hc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hc)
print("scikit-bio %s" % skbio.__version__)

matrix, variants, samples = hc.parse_vcf_region(VCF, CHROM, 2, 1.0, None)
res = hc.perform_clustering(matrix, samples, (2, 30), 5, "cosine_nj")
fresh_raw = res["nj_tree"]
saved_raw = TreeNode.read(os.path.join(DEC, CHROM, "%s_nj.nwk" % CHROM))
S = set(res["samples"])


def splits(tree):
    tips = sorted(t.name for t in tree.tips())
    anchor = tips[0]
    out = set()
    for n in tree.traverse(include_self=False):
        if n.is_tip():
            continue
        u = n.subset() & S
        if 2 <= len(u) <= len(S) - 2:
            out.add(frozenset(u) if anchor not in u else frozenset(S - u))
    return out


a, b = splits(fresh_raw), splits(saved_raw)
print("\nUNROOTED splits: rebuilt %d, shipped %d, shared %d (%.1f%% of shipped)"
      % (len(a), len(b), len(a & b), 100 * len(a & b) / max(len(b), 1)))

la = sum(float(n.length or 0) for n in fresh_raw.traverse(include_self=False))
lb = sum(float(n.length or 0) for n in saved_raw.traverse(include_self=False))
print("total branch length: rebuilt %.6f, shipped %.6f" % (la, lb))
neg_a = sum(1 for n in fresh_raw.traverse(include_self=False) if (n.length or 0) < 0)
neg_b = sum(1 for n in saved_raw.traverse(include_self=False) if (n.length or 0) < 0)
print("negative branch lengths: rebuilt %d, shipped %d" % (neg_a, neg_b))

fr, sr = fresh_raw.root_at_midpoint(), saved_raw.root_at_midpoint()
fk = [frozenset(t.name for t in c.tips()) for c in fr.children]
sk = [frozenset(t.name for t in c.tips()) for c in sr.children]
print("\nmidpoint root splits the tips %s in the rebuild and %s in the shipped tree"
      % ("/".join(str(len(x)) for x in fk), "/".join(str(len(x)) for x in sk)))
print("the two root partitions agree: %s" % (set(fk) == set(sk)))

bv_fresh = hc.extract_branch_vectors(fresh_raw, res["samples"], min_tips=5)
recs, _, _, stats = hc.decompose_svs(res["matrix"], bv_fresh, variants, res["samples"],
                                     r2_threshold=0.3, min_carriers=5)
new = pd.DataFrame([r for r in recs if r["type"] == "monophyletic"])
old = pd.read_csv(os.path.join(DEC, CHROM, "%s_decomposition.tsv" % CHROM), sep="\t")
old = old[old.type == "monophyletic"].copy()
old["carriers"] = old.carrier_indices.fillna("").apply(
    lambda s: frozenset(int(i) for i in s.split(",") if s))
new["carriers"] = new.carrier_indices.apply(frozenset)

print("\nmonophyletic components: rebuilt %d, shipped %d" % (len(new), len(old)))
print("SVs resolved: rebuilt %d, shipped %d"
      % (new.sv_id.nunique(), old.sv_id.nunique()))
print("SVs with more than one origin: rebuilt %d, shipped %d"
      % (int((new.sv_id.value_counts() >= 2).sum()),
         int((old.sv_id.value_counts() >= 2).sum())))

byid = {}
for r in new.itertuples():
    byid.setdefault(r.sv_id, []).append(r.carriers)
j, exact = [], 0
for r in old.itertuples():
    cands = byid.get(r.sv_id, [])
    best = max((len(r.carriers & c) / float(len(r.carriers | c)) for c in cands), default=0.0)
    j.append(best)
    exact += best == 1.0
j = np.asarray(j)
print("\ncarrier sets of the shipped components, recovered by the rebuild:")
print("  identical            %4d of %d (%.1f%%)" % (exact, len(j), 100 * exact / len(j)))
print("  Jaccard >= 0.9       %4d (%.1f%%)" % (int((j >= .9).sum()), 100 * (j >= .9).mean()))
print("  Jaccard >= 0.5       %4d (%.1f%%)" % (int((j >= .5).sum()), 100 * (j >= .5).mean()))
print("  median Jaccard       %.3f" % float(np.median(j)))
print("  no counterpart       %4d" % int((j == 0).sum()))
