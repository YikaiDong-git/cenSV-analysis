#!/usr/bin/env python
"""Every multi-origin SV, described by where its origins sit on the centromere tree."""
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
TRACKS = os.path.join(ROOT, "cenSV_web/etl")
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


build_trees = _load(BT, "build_trees")
sys.path.insert(0, TRACKS)
import tracks                                                              # noqa: E402


def clade_tips(node):
    return frozenset(t.name for t in node.tips())


def path_to_root(node):
    out = []
    while node is not None:
        out.append(node)
        node = node.parent
    return out


def relate(a, b):
    ta, tb = clade_tips(a), clade_tips(b)
    if ta <= tb:
        kind = "nested"
    elif tb <= ta:
        kind = "nested"
    elif ta & tb:
        kind = "overlapping"
    elif a.parent is b.parent:
        kind = "sister"
    else:
        kind = "disjoint"
    up_a = {id(n): d for d, n in enumerate(path_to_root(a))}
    mrca, depth_b = None, None
    for d, n in enumerate(path_to_root(b)):
        if id(n) in up_a:
            mrca, depth_b = n, d
            break
    depth_a = up_a[id(mrca)]

    def blen(node, stop):
        s = 0.0
        while node is not stop and node is not None:
            s += float(node.length or 0.0)
            node = node.parent
        return s

    return kind, depth_a, depth_b, blen(a, mrca) + blen(b, mrca), len(ta & tb)


rows = []
comp_rows = []
for chrom in CHROMS:
    rec, why = build_trees.build(chrom)
    if why:
        sys.exit("tree reconstruction failed: %s" % why)

    tree = TreeNode.read(os.path.join(DEC, chrom, "%s_nj.nwk" % chrom)).root_at_midpoint()
    nodes = list(tree.traverse(include_self=False))
    for b, n in rec["branch_n"].items():
        got = len(clade_tips(nodes[int(b[1:])]))
        assert got == n, ("%s %s spans %d tips here and %d in the validated reconstruction; "
                          "the two traversals disagree" % (chrom, b, got, n))
    node_of = {b: nodes[int(b[1:])] for b in rec["branch_n"]}

    hg = pd.read_csv(os.path.join(DEC, chrom, "%s_haplogroups.tsv" % chrom), sep="\t")
    haps = hg["sample"].astype(str).tolist()
    hgroup = dict(zip(haps, hg.haplogroup_id.astype(str)))

    d = pd.read_csv(os.path.join(DEC, chrom, "%s_decomposition.tsv" % chrom), sep="\t")
    d["carriers"] = d.carrier_indices.fillna("").apply(
        lambda s: [haps[int(i)] for i in s.split(",") if s])

    mono = d[d.type == "monophyletic"]
    outside = 0
    for r in mono.itertuples():
        outside += len(set(r.carriers) - clade_tips(node_of[str(r.branch_id)]))
    assert outside == 0, (
        "%s: %d carriers of monophyletic components fall outside their own branch's clade, so "
        "carrier_indices does not index the haplogroup table's row order" % (chrom, outside))

    v = pd.read_csv(os.path.join(DEC, chrom, "%s_sv_variants.tsv" % chrom), sep="\t",
                    low_memory=False,
                    usecols=["chrom", "pos", "end", "id", "svtype", "svlen", "ac", "af",
                             "original_sv_id", "mapped_haplogroup", "sv_category"])
    v = v[v.sv_category == "monophyletic"]

    n_or = v.groupby("original_sv_id").size()
    multi = set(n_or[n_or >= 2].index)
    if not multi:
        continue

    mono_by_sv = {}
    for r in mono.itertuples():
        mono_by_sv.setdefault(r.sv_id, []).append(r)
    resid = {r.sv_id: r for r in d[d.type == "residual"].itertuples()}

    for sv in sorted(multi):
        parts = v[v.original_sv_id == sv].sort_values("ac", ascending=False)
        logs = {str(r.branch_id): r for r in mono_by_sv.get(sv, [])}
        branches = [str(x) for x in parts.mapped_haplogroup]
        lead = branches[0]
        allc = set()
        for b in branches:
            if b in logs:
                allc |= set(logs[b].carriers)
        res = resid.get(sv)
        n_res = len(res.carriers) if res is not None else 0

        for i, r in enumerate(parts.itertuples()):
            b = str(r.mapped_haplogroup)
            lg = logs.get(b)
            car = set(lg.carriers) if lg is not None else set()
            anc = [c.split("_hap")[0] for c in car]
            comp_rows.append(dict(
                chrom=chrom, original_sv_id=sv, id=r.id, branch=b, rank=i,
                pos=int(r.pos), svtype=r.svtype, svlen=int(r.svlen), ac=int(r.ac),
                iteration=int(lg.iteration) if lg is not None else -1,
                branch_r2=float(lg.r2) if lg is not None else np.nan,
                n_carriers=len(car), clade_tips=len(clade_tips(node_of[b])),
                purity=len(car) / float(len(clade_tips(node_of[b]))),
                haplogroups=",".join(sorted({hgroup[c] for c in car})),
                individuals=len(set(anc)),
                homozygous=len(anc) - len(set(anc)),
                carriers=";".join(sorted(car)),
            ))

        kinds, seps, depths = [], [], []
        for j in range(1, len(branches)):
            k, da, db, bl, ov = relate(node_of[lead], node_of[branches[j]])
            kinds.append(k)
            seps.append(bl)
            depths.append(da + db)
        rows.append(dict(
            chrom=chrom, original_sv_id=sv, pos=int(parts.pos.iloc[0]),
            svtype=parts.svtype.iloc[0], svlen=int(parts.svlen.iloc[0]),
            n_origins=len(parts), branches=",".join(branches),
            carriers_total=len(allc), residual_carriers=n_res,
            ac_by_origin=",".join(str(int(a)) for a in parts.ac),
            r2_by_origin=",".join("%.3f" % (logs[b].r2 if b in logs else np.nan)
                                  for b in branches),
            relation=",".join(kinds), node_steps=",".join(str(x) for x in depths),
            branch_sep=",".join("%.4f" % s for s in seps),
            min_relation=("nested" if "nested" in kinds else
                          ("sister" if "sister" in kinds else "disjoint")),
            max_sep=max(seps) if seps else np.nan,
        ))
    print("  %-6s %4d multi-origin SVs" % (chrom, len(multi)))

svs = pd.DataFrame(rows)
comps = pd.DataFrame(comp_rows)

HAND = {"HG002": "EUR"}

samples = sorted({c.split("_hap")[0] for row in comps.carriers for c in row.split(";") if c})
sup, basis, hand, panel_pop, why = tracks.resolve_ancestry(
    [s for s in samples if s not in HAND])
assert why is None, "ancestry join failed: %s" % why
used = set(samples) & set(HAND)
assert used == set(HAND), ("the hand ancestry entries %s are not carriers of any multi-origin "
                           "component; delete them rather than leaving them unread"
                           % sorted(set(HAND) - used))
sup.update(HAND)
print("\nancestry: %d samples from the staged tables, %d by hand (%s)"
      % (len(sup) - len(HAND), len(HAND), ", ".join("%s=%s" % kv for kv in HAND.items())))

comps["ancestry"] = comps.carriers.apply(
    lambda s: ",".join("%s:%d" % (k, v) for k, v in sorted(
        pd.Series([sup[c.split("_hap")[0]] for c in s.split(";") if c]).value_counts().items(),
        key=lambda kv: -kv[1])))
comps["n_ancestries"] = comps.ancestry.apply(lambda s: len([x for x in s.split(",") if x]))
print("components carrying no decomposition-log row (no carrier list on disk): %d"
      % int((comps.iteration < 0).sum()))

os.makedirs(WU + "/output", exist_ok=True)
svs.to_csv(WU + "/output/multiorigin_tree_geometry.tsv", sep="\t", index=False)
comps.drop(columns=["carriers"]).to_csv(
    WU + "/output/multiorigin_components_detail.tsv", sep="\t", index=False)
comps.to_csv(WU + "/output/multiorigin_components_carriers.tsv", sep="\t", index=False)

assert len(svs) == 135, ("%d multi-origin SVs here against the 135 already reported; the two "
                         "counts read the same files and must agree" % len(svs))
print("\n%d multi-origin SVs, %d components" % (len(svs), len(comps)))
print("\n--- how the origins of one SV sit relative to each other ---")
print(svs.min_relation.value_counts().to_string())
print("\n--- pairwise relation over all origin pairs ---")
print(pd.Series([k for s in svs.relation for k in s.split(",")]).value_counts().to_string())
print("\n--- ancestries spanned by a component ---")
print(comps.groupby("iteration").n_ancestries.agg(["size", "mean", "median"]).to_string())
print("\nwrote %s" % (WU + "/output/multiorigin_tree_geometry.tsv"))
