#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""How often would a random carrier set of the same size be declared another origin?"""
import collections
import os
import sys

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

sys.path.insert(0, CENSV_ROOT + "/workflow/scripts")
from haplogroup_clustering import extract_branch_vectors  # noqa: E402

import skbio  # noqa: E402

ROOT = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
DEC = os.path.join(ROOT, "sv_extraction/decompose", SUF)
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
OUT = os.path.join(ROOT, "analyses/cdr_qtl_full_panel/output",
                   "origin_null_calibration.tsv")

R2_THRESHOLD = 0.3
MIN_CARRIERS = 5
NPERM = 20000
SEED = 20260911
TAB = "\t"


def r2_many(P, Ct, n):
    k = P.sum(axis=1).astype(np.float64)[:, None]
    m = Ct.sum(axis=0).astype(np.float64)[None, :]
    co = P.astype(np.float64) @ Ct.astype(np.float64)
    num = n * co - k * m
    den = np.sqrt(np.maximum((n * k - k * k) * (n * m - m * m), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(den > 0, num / den, 0.0)
    return r * r


def carriers_of(cell):
    s = str(cell)
    if s in ("nan", "", "None"):
        return []
    return [int(v) for v in s.split(",") if v.strip() != ""]


rows = []
orders_used = set()
rng = np.random.default_rng(SEED)
checked = 0

for chrom in CHROMS:
    nwk = os.path.join(DEC, chrom, "%s_nj.nwk" % chrom)
    log = os.path.join(DEC, chrom, "%s_decomposition.tsv" % chrom)
    if not (os.path.exists(nwk) and os.path.exists(log)):
        sys.exit("FAIL missing decomposition inputs for %s" % chrom)

    tree = skbio.TreeNode.read(nwk)
    tip_names = [str(t.name) for t in tree.tips()]
    n = len(tip_names)

    d = pd.read_csv(log, sep=TAB)
    per_sv = collections.defaultdict(list)
    for _, r in d.iterrows():
        per_sv[str(r.sv_id)].append((int(r.iteration), str(r.type), str(r.branch_id),
                                     carriers_of(r.carrier_indices),
                                     None if pd.isna(r.r2) else float(r.r2)))
    tasks = []
    for sv_id, rs in per_sv.items():
        rs.sort(key=lambda x: x[0])
        suffix = [set() for _ in rs]
        acc = set()
        for j in range(len(rs) - 1, -1, -1):
            acc = acc | set(rs[j][3])
            suffix[j] = set(acc)
        rank = 0
        for j, (it, typ, br, car, r2) in enumerate(rs):
            if typ != "monophyletic":
                continue
            tasks.append((br, suffix[j], r2, rank, sv_id, len(car)))
            rank += 1

    hgf = os.path.join(DEC, chrom, "%s_haplogroups.tsv" % chrom)
    cands = []
    if os.path.exists(hgf):
        h = pd.read_csv(hgf, sep=TAB)
        col = "sample" if "sample" in h.columns else h.columns[0]
        cands.append(("haplogroups.tsv order", [str(x) for x in h[col]]))
    cands += [("tree tip order", tip_names), ("sorted tip names", sorted(tip_names))]

    picked, tried = None, []
    for name, samples in cands:
        if sorted(samples) != sorted(tip_names):
            tried.append("%s (different tip set)" % name)
            continue
        bvs = extract_branch_vectors(tree, samples, min_tips=MIN_CARRIERS)
        Ct = np.column_stack([b["vector"] for b in bvs]).astype(np.int8)
        bid = {b["branch_id"]: j for j, b in enumerate(bvs)}
        worst, why = 0.0, None
        for br, res, want, _, _, _ in tasks:
            if br not in bid:
                why = "logged branch %s absent from the rebuilt set" % br
                break
            x = np.zeros((1, n), dtype=np.int8)
            x[0, sorted(res)] = 1
            worst = max(worst, abs(float(r2_many(x, Ct[:, [bid[br]]], n)[0, 0]) - want))
            if worst > 1e-8:
                why = "worst |dr2| %.3e" % worst
                break
        if why is None:
            picked = (name, samples, bvs, Ct, bid, worst)
            break
        tried.append("%s (%s)" % (name, why))

    if picked is None:
        sys.exit("FAIL %s: no candidate haplotype order reproduces the logged R2 values, so "
                 "the null would sit on a different instrument than the run. Tried: %s"
                 % (chrom, "; ".join(tried)))
    order_name, samples, bvs, Ct, bid, worst = picked
    checked += len(tasks)
    orders_used.add(order_name)
    clade_size = Ct.sum(axis=0)

    for br, res, want, rank, sv_id, n_comp in tasks:
        if rank == 0 or not res:
            continue
        k = len(res)
        P = np.zeros((NPERM, n), dtype=np.int8)
        for b in range(NPERM):
            P[b, rng.choice(n, size=k, replace=False)] = 1
        sc = r2_many(P, Ct, n)
        best_t = sc.argmax(axis=1)
        best = sc[np.arange(NPERM), best_t]
        inter = (P * Ct[:, best_t].T).sum(axis=1)
        declared = (best >= R2_THRESHOLD) & (inter >= MIN_CARRIERS)
        rows.append({
            "chrom": chrom, "sv_id": str(sv_id)[:60], "rank": rank, "branch_id": br,
            "n_residual": k, "n_component": n_comp, "clade_tips": int(clade_size[bid[br]]),
            "observed_r2": want, "n_branches": len(bvs), "n_tip": n,
            "null_median_max_r2": float(np.median(best)),
            "null_p95_max_r2": float(np.quantile(best, 0.95)),
            "p_random_declared": float(declared.mean()),
            "p_value": float((best >= want).mean()),
        })
    print("%s  order=%s  branches %d  tips %d  components %d (worst |dr2| %.1e)  "
          "secondary so far %d"
          % (chrom, order_name, len(bvs), n, len(tasks), worst, len(rows)), flush=True)

t = pd.DataFrame(rows)
t.to_csv(OUT, sep=TAB, index=False)
print()
print("instrument: %d logged monophyletic components recomputed, every one matching the log"
      % checked)
print("haplotype orders that reproduced the run: %s" % sorted(orders_used))
print("later origins tested: %d across %d chromosomes" % (len(t), t.chrom.nunique()))
print()
print("How often a RANDOM residual of the same size is DECLARED another origin")
print("(max R2 over the same branches >= 0.3 and the intersection holding >= 5 carriers):")
q = t["p_random_declared"]
print("  median %.3f   q1 %.3f   q3 %.3f   min %.3f   max %.3f"
      % (q.median(), q.quantile(.25), q.quantile(.75), q.min(), q.max()))
for thr in (0.05, 0.20, 0.50, 0.80):
    print("  later origins where a random residual is declared more than %3.0f%% of the time: "
          "%d of %d" % (100 * thr, int((q > thr).sum()), len(t)))
print()
print("p-value of the OBSERVED R2 against the same branch search:")
p = t["p_value"]
for thr in (0.05, 0.10, 0.50):
    print("  p > %.2f : %d of %d" % (thr, int((p > thr).sum()), len(t)))
print("  p <= 0.05 : %d of %d" % (int((p <= 0.05).sum()), len(t)))
print()
print("by the size of the component the origin produced:")
t["bucket"] = pd.cut(t.n_component, [0, 9, 19, 10 ** 9], labels=["5-9", "10-19", "20+"])
print(t.groupby("bucket", observed=True).agg(
    n=("p_value", "size"), median_obs_r2=("observed_r2", "median"),
    median_p=("p_value", "median"),
    median_random_declared=("p_random_declared", "median")).to_string())
print()
print("wrote %s" % OUT)
