#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import importlib.util
import os

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
HC = os.path.join(ROOT, "workflow/scripts/haplogroup_clustering.py")
KEY = ["pos", "end", "svtype", "svlen", "ac", "an", "n_miss"]

spec = importlib.util.spec_from_file_location("hc", HC)
hc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hc)


def clades(tree, samples):
    S = set(samples)
    out = set()
    for n in tree.traverse(include_self=False):
        if n.is_tip():
            continue
        u = n.subset() & S
        if 5 <= len(u) <= len(S) - 5:
            out.add(frozenset(u))
    return out


BT = os.path.join(ROOT, "discovery/2026-08-30_censv_individual_views/code/build/build_trees.py")
_bt = importlib.util.spec_from_file_location("build_trees", BT)
build_trees = importlib.util.module_from_spec(_bt)
_bt.loader.exec_module(build_trees)


def carrying_clades(C):
    rec, why = build_trees.build(C)
    if why:
        raise SystemExit(why)
    tips, L = rec["tips"], rec["layout"]
    var = pd.read_csv("%s/%s/%s_sv_variants.tsv" % (DEC, C, C), sep="\t", low_memory=False,
                      usecols=["sv_category", "mapped_haplogroup"])
    used = set(var[var.sv_category == "monophyletic"].mapped_haplogroup.dropna().astype(str))
    out = set()
    for b in used:
        node = rec["branch_node"][b]
        lo, hi = L["tmin"][node], L["tmax"][node]
        out.add(frozenset(tips[lo:hi + 1]))
    return out


def run(C):
    print("=" * 70)
    print(C)
    matrix, variants, samples = hc.parse_vcf_region(VCF, C, 2, 1.0, None)
    now = pd.DataFrame([{k: v[k] for k in KEY} for v in variants])
    now["col"] = np.arange(len(now))
    then = pd.read_csv("%s/%s/%s_variants.tsv" % (DEC, C, C), sep="\t", usecols=KEY,
                       low_memory=False)
    print("  parsed today %d SVs; the run logged %d" % (len(now), len(then)))

    for df, lab in ((now, "today"), (then, "the run")):
        d = df.duplicated(subset=KEY).sum()
        print("  %-8s rows sharing a key with another row: %d" % (lab, d))

    j = then.merge(now, on=KEY, how="left")
    lost = int(j.col.isna().sum())
    print("  run-era SVs with no row in today's parse: %d" % lost)
    extra = len(now) - len(j.dropna(subset=["col"]).col.unique())
    print("  today's SVs the run did not keep:          %d" % extra)
    if lost:
        print("  NOT recoverable by restriction: %d of the run's SVs are absent from the "
              "current VCF" % lost)
        miss = j[j.col.isna()][KEY]
        gone = now.merge(then.assign(seen=1), on=KEY, how="left")
        gone = gone[gone.seen.isna()][KEY]
        print("\n  the run kept these, today's parse does not:")
        print(miss.to_string(index=False))
        print("\n  today's parse keeps these, the run did not:")
        print(gone.to_string(index=False))
        same = set(map(tuple, miss[["pos", "end", "svtype"]].values)) & \
            set(map(tuple, gone[["pos", "end", "svtype"]].values))
        print("\n  positions appearing on both sides: %d of %d" % (len(same), len(miss)))
        pm = set(map(tuple, miss[["pos", "ac"]].values))
        pg = set(map(tuple, gone[["pos", "ac"]].values))
        print("  (position, allele count) pairs appearing on both sides: %d of %d"
              % (len(pm & pg), len(pm)))

        keep = j.dropna(subset=["col"]).col.astype(int).values
        hg = pd.read_csv("%s/%s/%s_haplogroups.tsv" % (DEC, C, C), sep="\t")
        TAB = hg["sample"].astype(str).tolist()
        at = {s: i for i, s in enumerate(samples)}
        order = TAB + [s for s in samples if s not in set(TAB)]
        saved = TreeNode.read("%s/%s/%s_nj.nwk" % (DEC, C, C)).root_at_midpoint()
        used = carrying_clades(C)
        for lab, cols in (("all of today's", np.arange(len(now))),
                          ("the agreed subset", keep)):
            m = matrix[[at[s] for s in order]][:, cols]
            res = hc.perform_clustering(m, order, (2, 30), 5, "cosine_nj")
            mine = clades(res["nj_tree"], res["samples"])
            ship = clades(saved, res["samples"])
            print("  %-18s %5d SVs -> %3d clades, %3d shared with the shipped %d"
                  % (lab, m.shape[1], len(mine), len(mine & ship), len(ship)))
            print("  %-18s of the %d clades carrying a component, %d survive"
                  % ("", len(used), len(used & mine)))
        return

    cols = j.col.astype(int).values
    hg = pd.read_csv("%s/%s/%s_haplogroups.tsv" % (DEC, C, C), sep="\t")
    TAB = hg["sample"].astype(str).tolist()
    at = {s: i for i, s in enumerate(samples)}
    order = TAB + [s for s in samples if s not in set(TAB)]
    m = matrix[[at[s] for s in order]][:, cols]
    print("  restricted matrix: %d haplotypes x %d SVs" % m.shape)

    res = hc.perform_clustering(m, order, (2, 30), 5, "cosine_nj")
    mine = clades(res["nj_tree"], res["samples"])
    saved = TreeNode.read("%s/%s/%s_nj.nwk" % (DEC, C, C)).root_at_midpoint()
    ship = clades(saved, res["samples"])
    lb = sum(float(n.length or 0) for n in res["nj_tree"].traverse(include_self=False))
    ls = sum(float(n.length or 0) for n in TreeNode.read(
        "%s/%s/%s_nj.nwk" % (DEC, C, C)).traverse(include_self=False))
    print("  clades rebuilt %d, shipped %d, shared %d" % (len(mine), len(ship),
                                                          len(mine & ship)))
    print("  total branch length rebuilt %.12g, shipped %.12g" % (lb, ls))
    print("  REPRODUCES" if mine == ship else "  STILL DIFFERS")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default="chr15,chr18")
    a = ap.parse_args()
    for C in a.chroms.split(","):
        run(C)
