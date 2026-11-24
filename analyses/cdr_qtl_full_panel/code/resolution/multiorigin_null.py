#!/usr/bin/env python
"""Is a second origin farther from the first than two clades picked at chance would be?"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from skbio.tree import TreeNode
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_"
       "qgap200_rgap200")
DEC = os.path.join(ROOT, "sv_extraction/decompose", SUF)
WU = os.path.join(ROOT, "analyses/cdr_qtl_full_panel")
BT = os.path.join(ROOT, "discovery/2026-08-30_censv_individual_views/code/build/build_trees.py")
DRAWS = 2000
TOL = 0.25
MIN_OUTSIDE = 5
SEED = 20260910

spec = importlib.util.spec_from_file_location("build_trees", BT)
build_trees = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build_trees)

geo = pd.read_csv(WU + "/output/multiorigin_tree_geometry.tsv", sep="\t")
rng = np.random.default_rng(SEED)


def tips_of(node):
    return frozenset(t.name for t in node.tips())


def to_root(node):
    out = []
    while node is not None:
        out.append(node)
        node = node.parent
    return out


def pair_stats(a, b):
    ta, tb = tips_of(a), tips_of(b)
    if ta <= tb or tb <= ta:
        kind = "nested"
    elif ta & tb:
        kind = "overlapping"
    elif a.parent is b.parent:
        kind = "sister"
    else:
        kind = "disjoint"
    up = {id(n): d for d, n in enumerate(to_root(a))}
    mrca = next(n for n in to_root(b) if id(n) in up)

    def blen(node):
        s = 0.0
        while node is not mrca and node is not None:
            s += float(node.length or 0.0)
            node = node.parent
        return s

    return kind, blen(a) + blen(b)


rows = []
for chrom, sub in geo.groupby("chrom"):
    rec, why = build_trees.build(chrom)
    if why:
        sys.exit("tree reconstruction failed: %s" % why)
    tree = TreeNode.read(os.path.join(DEC, chrom, "%s_nj.nwk" % chrom)).root_at_midpoint()
    nodes = list(tree.traverse(include_self=False))
    node_of = {b: nodes[int(b[1:])] for b in rec["branch_n"]}
    tipset = {b: tips_of(n) for b, n in node_of.items()}
    size = {b: len(v) for b, v in tipset.items()}
    ids = list(node_of)

    for r in sub.itertuples():
        bs = r.branches.split(",")
        for j in range(1, len(bs)):
            b1, b2 = bs[0], bs[j]
            kind, sep = pair_stats(node_of[b1], node_of[b2])
            n1, n2 = size[b1], size[b2]
            pool1 = [b for b in ids if abs(size[b] - n1) <= TOL * n1]
            pool2 = [b for b in ids if abs(size[b] - n2) <= TOL * n2]
            seps, kinds = [], []
            tries = 0
            while len(seps) < DRAWS and tries < DRAWS * 40:
                tries += 1
                c1 = pool1[rng.integers(len(pool1))]
                c2 = pool2[rng.integers(len(pool2))]
                if c1 == c2 or len(tipset[c2] - tipset[c1]) < MIN_OUTSIDE:
                    continue
                k, s = pair_stats(node_of[c1], node_of[c2])
                kinds.append(k)
                seps.append(s)
            seps = np.asarray(seps)
            rows.append(dict(
                chrom=chrom, original_sv_id=r.original_sv_id, pos=r.pos, svtype=r.svtype,
                branch_a=b1, branch_b=b2, n_tips_a=n1, n_tips_b=n2,
                relation=kind, separation=sep, n_draws=len(seps),
                null_median=float(np.median(seps)), null_p10=float(np.quantile(seps, .1)),
                null_p90=float(np.quantile(seps, .9)),
                percentile=float((seps < sep).mean()),
                null_frac_nested=float(np.mean([k == "nested" for k in kinds])),
                null_frac_sister=float(np.mean([k == "sister" for k in kinds])),
                null_frac_disjoint=float(np.mean([k == "disjoint" for k in kinds]))))
    print("  %-6s %3d pairs" % (chrom, sum(len(x.split(",")) - 1 for x in sub.branches)))

out = pd.DataFrame(rows)
out.to_csv(WU + "/output/multiorigin_null.tsv", sep="\t", index=False)
assert len(out) == 144, "%d origin pairs here against the 144 already reported" % len(out)

print("\n%d origin pairs" % len(out))
print("observed relation: %s" % out.relation.value_counts().to_dict())
print("under the size-matched null a pair is nested %.1f%%, sister %.1f%%, disjoint %.1f%% "
      "of the time (mean over pairs)"
      % (100 * out.null_frac_nested.mean(), 100 * out.null_frac_sister.mean(),
         100 * out.null_frac_disjoint.mean()))
exp_ns = float((out.null_frac_nested + out.null_frac_sister).sum())
print("expected nested-or-sister pairs: %.1f of %d; observed 0" % (exp_ns, len(out)))
print("\nseparation percentile against the null:")
print(out.percentile.describe().to_string())
print("pairs above the null median: %d of %d" % (int((out.percentile > .5).sum()), len(out)))
print("\nwrote %s" % (WU + "/output/multiorigin_null.tsv"))
