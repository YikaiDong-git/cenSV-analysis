#!/usr/bin/env python
"""Rebuild each centromere tree WITHOUT the SV it was used to explain, and look again."""
import argparse
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
VCF = os.path.join(ROOT, "variant_calling/svim_asm", "cohort_svim_asm_%s.vcf.gz" % SUF)
WU = os.path.join(ROOT, "analyses/cdr_qtl_full_panel")
OUT = os.path.join(WU, "output/resolution/loo")
HC = os.path.join(ROOT, "workflow/scripts/haplogroup_clustering.py")

R2T, MIN_CARRIERS, MIN_HG, MIN_TIPS, MIN_SAMPLES = 0.3, 5, 10, 5, 2

spec = importlib.util.spec_from_file_location("hc", HC)
hc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hc)

ap = argparse.ArgumentParser()
ap.add_argument("--chrom", required=True)
ap.add_argument("--limit", type=int, default=0)
a = ap.parse_args()
os.makedirs(OUT, exist_ok=True)
C = a.chrom


def branches_of(matrix, samples):
    res = hc.perform_clustering(matrix, samples, (2, 30), 5, "cosine_nj")
    bv = hc.extract_branch_vectors(res["nj_tree"], res["samples"], min_tips=MIN_TIPS)
    for b in bv:
        b["clade"] = frozenset(s for s, v in zip(res["samples"], b["vector"]) if v)
    return bv, res


def resolve(vec, branches):
    out, sv, it = [], vec.copy(), 0
    while True:
        ac = int(sv.sum())
        if ac < MIN_CARRIERS:
            return out
        r2 = np.array([hc._r2_binary(sv, b["vector"]) for b in branches])
        t = int(np.argmax(r2))
        if r2[t] < R2T:
            return out
        mono = sv & branches[t]["vector"]
        n = int(mono.sum())
        if (it == 0 and n < MIN_HG) or n < MIN_CARRIERS:
            return out
        out.append(dict(r2=float(r2[t]), n=n, tips=len(branches[t]["clade"]),
                        carriers=frozenset(np.where(mono > 0)[0].tolist())))
        sv = sv & (~branches[t]["vector"].astype(bool)).astype(np.int8)
        it += 1


matrix, variants, samples = hc.parse_vcf_region(VCF, C, MIN_SAMPLES, 1.0, None)
ids = [v["id"] for v in variants]
col = {v: j for j, v in enumerate(ids)}

hg = pd.read_csv("%s/%s/%s_haplogroups.tsv" % (DEC, C, C), sep="\t")
TAB = hg["sample"].astype(str).tolist()
assert set(TAB) <= set(samples), (
    "%s: %d haplotypes in the haplogroup table are absent from the VCF"
    % (C, len(set(TAB) - set(samples))))
at = {s: i for i, s in enumerate(samples)}
order = TAB + [s for s in samples if s not in set(TAB)]
matrix = matrix[[at[s] for s in order]]
samples = order

full, res = branches_of(matrix, samples)
X, HAPS = res["matrix"], res["samples"]
assert HAPS == TAB, (
    "%s: after permuting into the haplogroup table's order the clustered haplotypes still "
    "differ from it at %d positions"
    % (C, sum(1 for a, b in zip(HAPS, TAB) if a != b) + abs(len(HAPS) - len(TAB))))
print("%s: %d haplotypes x %d SVs, %d candidate branches"
      % (C, X.shape[0], X.shape[1], len(full)))

saved = TreeNode.read("%s/%s/%s_nj.nwk" % (DEC, C, C)).root_at_midpoint()
S = set(HAPS)
shipped = {}
for i, n in enumerate(saved.traverse(include_self=False)):
    if n.is_tip():
        continue
    u = n.subset() & S
    if MIN_TIPS <= len(u) <= len(HAPS) - MIN_TIPS:
        shipped["B%d" % i] = frozenset(u)
mine = {b["clade"]: b for b in full}
assert set(shipped.values()) == set(mine), (
    "%s: the rebuilt tree holds %d clades and the shipped newick %d, sharing %d; the rebuild "
    "is not the tree the components were defined on"
    % (C, len(mine), len(set(shipped.values())),
       len(set(shipped.values()) & set(mine))))

d = pd.read_csv("%s/%s/%s_decomposition.tsv" % (DEC, C, C), sep="\t")
d["carriers"] = d.carrier_indices.fillna("").apply(
    lambda s: frozenset(int(i) for i in s.split(",") if s))
mono = d[d.type == "monophyletic"]
worst = 0.0
for r in mono.itertuples():
    if r.sv_id not in col:
        continue
    sv = X[:, col[r.sv_id]].copy()
    for prev in mono[(mono.sv_id == r.sv_id) & (mono.iteration < r.iteration)].itertuples():
        sv = sv & (~mine[shipped[str(prev.branch_id)]]["vector"].astype(bool)).astype(np.int8)
    worst = max(worst,
                abs(hc._r2_binary(sv, mine[shipped[str(r.branch_id)]]["vector"]) - float(r.r2)))
assert worst < 1e-9, ("%s: the rebuilt tree reproduces recorded branch R2 only to %.3g" % (C, worst))
print("  precondition: %d clades and %d component R2 values reproduce exactly"
      % (len(mine), len(mono)))

recs, _, _, _ = hc.decompose_svs(X, full, variants, HAPS, r2_threshold=R2T,
                                 min_carriers=MIN_CARRIERS)
rebuilt = {}
for r in recs:
    if r["type"] == "monophyletic":
        rebuilt.setdefault(r["sv_id"], []).append(frozenset(r["carrier_indices"]))
nmono = sum(len(v) for v in rebuilt.values())
assert nmono == len(mono), ("%s: the rebuild yields %d monophyletic components against the "
                            "shipped %d" % (C, nmono, len(mono)))
same = sum(1 for r in mono.itertuples() if r.carriers in rebuilt.get(r.sv_id, []))
tied = 0
for r in mono.itertuples():
    if r.sv_id not in col:
        continue
    sv = X[:, col[r.sv_id]].copy()
    for prev in mono[(mono.sv_id == r.sv_id) & (mono.iteration < r.iteration)].itertuples():
        sv = sv & (~mine[shipped[str(prev.branch_id)]]["vector"].astype(bool)).astype(np.int8)
    vals = np.array([hc._r2_binary(sv, b["vector"]) for b in full])
    tied += int((vals >= vals.max() - 1e-12).sum() > 1)
print("  the same %d components with the same R2; %d of %d come back with identical carriers, "
      "and %d have more than one branch matching equally well"
      % (nmono, same, len(mono), tied))
TIE_RATE = tied / float(len(mono))

orig = {}
for r in mono.itertuples():
    orig.setdefault(r.sv_id, []).append(r)
multi = sorted(k for k, v in orig.items() if len(v) >= 2 and k in col)
single = sorted(k for k, v in orig.items() if len(v) == 1 and k in col)
if not multi:
    sys.exit("%s carries no multi-origin SV" % C)

size = {k: int(X[:, col[k]].sum()) for k in single}
ctrl, taken = [], set()
for k in multi:
    n = int(X[:, col[k]].sum())
    pick = min((c for c in single if c not in taken), key=lambda c: (abs(size[c] - n), c))
    taken.add(pick)
    ctrl.append(pick)

targets = [(k, "multi_origin") for k in multi] + [(k, "single_origin") for k in ctrl]
if a.limit:
    targets = targets[:a.limit]

rows = []
for n_done, (sv_id, arm) in enumerate(targets, 1):
    j = col[sv_id]
    vec = X[:, j].copy()
    loo, _ = branches_of(np.delete(X, j, axis=1), HAPS)
    found = resolve(vec, loo)
    was = sorted(orig[sv_id], key=lambda r: r.iteration)
    for i, r in enumerate(was):
        best = max(((len(r.carriers & f["carriers"]) / float(len(r.carriers | f["carriers"])), k)
                    for k, f in enumerate(found)), default=(0.0, -1))
        rows.append(dict(
            chrom=C, sv_id=sv_id, arm=arm, origin_rank=i,
            n_origins_before=len(was), n_origins_after=len(found),
            branch_before=str(r.branch_id), r2_before=float(r.r2),
            n_carriers_before=int(r.n_carriers),
            r2_after=(found[best[1]]["r2"] if best[1] >= 0 else np.nan),
            n_carriers_after=(found[best[1]]["n"] if best[1] >= 0 else 0),
            jaccard=best[0], recovered=bool(best[0] >= 0.9),
            tie_rate=TIE_RATE))
    if n_done % 10 == 0 or n_done == len(targets):
        print("  %3d/%d trees rebuilt" % (n_done, len(targets)))

out = pd.DataFrame(rows)
out.to_csv("%s/%s_loo.tsv" % (OUT, C), sep="\t", index=False)
print("\n--- carriers recovered on a tree that never saw the SV ---")
print(out.groupby("arm").jaccard.agg(["size", "mean", "median"]).to_string())
print(out.groupby("arm").recovered.mean().to_string())
print("\nwrote %s/%s_loo.tsv" % (OUT, C))
