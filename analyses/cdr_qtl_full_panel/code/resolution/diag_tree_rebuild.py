#!/usr/bin/env python
"""Does rebuilding the tree from the matrix give back the tree that was shipped?"""
import importlib.util
import os

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

matrix, variants, samples = hc.parse_vcf_region(VCF, CHROM, 2, 1.0, None)
res = hc.perform_clustering(matrix, samples, (2, 30), 5, "cosine_nj")
fresh = hc.extract_branch_vectors(res["nj_tree"], res["samples"], min_tips=5)

saved = TreeNode.read(os.path.join(DEC, CHROM, "%s_nj.nwk" % CHROM)).root_at_midpoint()
sset = set(samples)
old = {}
for i, n in enumerate(saved.traverse(include_self=False)):
    if n.is_tip():
        continue
    u = n.subset() & sset
    if 5 <= len(u) <= len(samples) - 5:
        old["B%d" % i] = frozenset(u)

new = {b["branch_id"]: frozenset(s for s, v in zip(samples, b["vector"]) if v) for b in fresh}

print("shipped tree: %d candidate branches   rebuilt tree: %d" % (len(old), len(new)))
print("ids in common: %d" % len(set(old) & set(new)))
same_id_same_clade = sum(1 for b in set(old) & set(new) if old[b] == new[b])
print("ids in common AND the same clade: %d" % same_id_same_clade)

oc, nc = set(old.values()), set(new.values())
print("\nclades (ignoring ids): shipped %d, rebuilt %d, shared %d"
      % (len(oc), len(nc), len(oc & nc)))
print("shipped clades recovered by the rebuild: %.1f%%" % (100 * len(oc & nc) / len(oc)))

miss = sorted(oc - nc, key=len)
print("\nshipped clades the rebuild does not contain: %d" % len(miss))
for c in miss[:8]:
    best = max((len(c & d) / float(len(c | d)), len(d)) for d in nc)
    print("   %3d tips -> best rebuilt match Jaccard %.3f (%d tips)" % (len(c), best[0], best[1]))

used = set()
for c in CHROM,:
    pass
import pandas as pd                                                        # noqa: E402
v = pd.read_csv(os.path.join(DEC, CHROM, "%s_sv_variants.tsv" % CHROM), sep="\t",
                low_memory=False, usecols=["sv_category", "mapped_haplogroup"])
used = set(v[v.sv_category == "monophyletic"].mapped_haplogroup.dropna().astype(str))
print("\nbranches that actually carry a component: %d" % len(used))
print("  of these, the rebuild reproduces the same clade under the same id: %d"
      % sum(1 for b in used if b in new and new[b] == old[b]))
print("  of these, the clade exists in the rebuild under ANY id: %d"
      % sum(1 for b in used if old[b] in nc))
