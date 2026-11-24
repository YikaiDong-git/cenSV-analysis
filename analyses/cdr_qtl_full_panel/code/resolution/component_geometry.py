#!/usr/bin/env python
"""Where on its own tree each monophyletic component sits, for all 2,247 of them."""
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
OUT = os.path.join(WU, "output/resolution")
BT = os.path.join(ROOT, "discovery/2026-08-30_censv_individual_views/code/build/build_trees.py")
TRACKS = os.path.join(ROOT, "cenSV_web/etl")
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
HAND = {"HG002": "EUR"}
SUP = ["AFR", "AMR", "EAS", "EUR", "SAS"]

spec = importlib.util.spec_from_file_location("build_trees", BT)
build_trees = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build_trees)
sys.path.insert(0, TRACKS)
import tracks                                                              # noqa: E402

rows, tsum = [], []
for C in CHROMS:
    rec, why = build_trees.build(C)
    if why:
        sys.exit("tree reconstruction failed: %s" % why)
    tree = TreeNode.read("%s/%s/%s_nj.nwk" % (DEC, C, C)).root_at_midpoint()
    nodes = list(tree.traverse(include_self=False))
    node_of = {b: nodes[int(b[1:])] for b in rec["branch_n"]}

    depth, steps = {}, {}
    for n in tree.preorder(include_self=True):
        p = n.parent
        depth[id(n)] = 0.0 if p is None else depth[id(p)] + float(n.length or 0.0)
        steps[id(n)] = 0 if p is None else steps[id(p)] + 1
    height = max(depth[id(t)] for t in tree.tips())

    hg = pd.read_csv("%s/%s/%s_haplogroups.tsv" % (DEC, C, C), sep="\t")
    haps = hg["sample"].astype(str).tolist()
    hgroup = dict(zip(haps, hg.haplogroup_id.astype(str)))

    d = pd.read_csv("%s/%s/%s_decomposition.tsv" % (DEC, C, C), sep="\t")
    d["carriers"] = d.carrier_indices.fillna("").apply(
        lambda s: [haps[int(i)] for i in s.split(",") if s])
    mono = d[d.type == "monophyletic"]
    for r in mono.itertuples():
        clade = {t.name for t in node_of[str(r.branch_id)].tips()}
        assert set(r.carriers) <= clade, (
            "%s %s: %d carriers fall outside their own branch's clade"
            % (C, r.branch_id, len(set(r.carriers) - clade)))

    tsum.append(dict(chrom=C, n_tips=len(haps), n_branch=rec["n_branch"], height=height,
                     n_haplogroup=hg.haplogroup_id.nunique(),
                     largest_haplogroup=int(hg.haplogroup_id.value_counts().iloc[0])))

    for r in mono.itertuples():
        b = str(r.branch_id)
        nd = node_of[b]
        clade = sorted(t.name for t in nd.tips())
        car = sorted(r.carriers)
        rows.append(dict(
            chrom=C, original_sv_id=r.sv_id, branch_id=b, iteration=int(r.iteration),
            branch_r2=float(r.r2), n_carriers=len(car), clade_tips=len(clade),
            purity=len(car) / float(len(clade)),
            clade_depth=depth[id(nd)], clade_depth_rel=depth[id(nd)] / height,
            node_steps=steps[id(nd)], branch_length=float(nd.length or 0.0),
            clade_haplogroups=len({hgroup[t] for t in clade}),
            carrier_haplogroups=len({hgroup[t] for t in car}),
            top_haplogroup_share=max(
                pd.Series([hgroup[t] for t in clade]).value_counts(normalize=True)),
            clade=";".join(clade), carriers=";".join(car)))
    print("  %-6s %4d components" % (C, len(mono)))

g = pd.DataFrame(rows)
samples = sorted({t.split("_hap")[0] for s in g.clade for t in s.split(";")})
sup, basis, byhand, panel_pop, why = tracks.resolve_ancestry(
    [s for s in samples if s not in HAND])
assert why is None, "ancestry join failed: %s" % why
sup.update(HAND)


def anc_counts(names):
    v = pd.Series([sup[t.split("_hap")[0]] for t in names.split(";")]).value_counts()
    return [int(v.get(k, 0)) for k in SUP]


for frame, pre in (("clade", "clade"), ("carriers", "carrier")):
    cc = np.array([anc_counts(s) for s in g[frame]])
    for j, k in enumerate(SUP):
        g["%s_%s" % (pre, k)] = cc[:, j]
    g["%s_n_ancestry" % pre] = (cc > 0).sum(axis=1)
    frac = cc / cc.sum(axis=1, keepdims=True)
    g["%s_top_ancestry_share" % pre] = frac.max(axis=1)
    g["%s_top_ancestry" % pre] = [SUP[i] for i in cc.argmax(axis=1)]

g.drop(columns=["clade", "carriers"]).to_csv(OUT + "/component_geometry.tsv", sep="\t",
                                             index=False)
pd.DataFrame(tsum).to_csv(OUT + "/tree_summary.tsv", sep="\t", index=False)

assert len(g) == 2247, "%d components here against the 2,247 already reported" % len(g)
print("\n%d components across %d chromosomes" % (len(g), g.chrom.nunique()))
print("\n--- clade purity: how much of the branch its carriers fill ---")
print(g.groupby(g.iteration.clip(upper=1)).purity.describe()[
    ["count", "25%", "50%", "75%"]].to_string())
print("\n--- ancestries a component's carriers span, by origin rank ---")
print(g.groupby(g.iteration.clip(upper=1)).carrier_n_ancestry.agg(
    ["size", "mean", "median"]).to_string())
print("\n--- how concentrated a clade is in one ancestry, against the cohort ---")
print(g.clade_top_ancestry_share.describe().to_string())
print("\nwrote %s/component_geometry.tsv" % OUT)
