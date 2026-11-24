#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""If the ribbon painted every origin, what quantity would it be showing?"""
import base64
import collections
import io
import json
import os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
DATA = os.path.join(P1, "cenSV_web/web/public/data")
DISC = os.path.join(P1, "analyses/cdr_qtl_full_panel/output")


def bits(b64, n):
    raw = base64.b64decode(b64)
    return {k for k in range(n) if raw[k >> 3] >> (7 - (k & 7)) & 1}


idx = json.load(io.open(os.path.join(DATA, "hap_index.json"), encoding="utf-8"))
cat = json.load(io.open(os.path.join(DATA, "catalog.json"), encoding="utf-8"))
ci = {k: i for i, k in enumerate(cat["columns"])}
G = lambda k, i: cat["data"][ci[k]][i]

by_coord = collections.defaultdict(dict)
meta = {}
for i in range(cat["rows"]):
    key = (G("chrom", i), G("pos", i), G("svtype", i), G("svlen", i))
    by_coord[key].setdefault(G("sv_category", i), []).append(i)
    meta[(G("chrom", i), G("sv_index", i))] = key

n_parent = n_no_parent = 0
gaps, exact = [], 0
for chrom in sorted(idx["multi"]):
    blk = json.load(io.open(os.path.join(DATA, "hap_%s.json" % chrom), encoding="utf-8"))
    nt = blk["tree"]["n_tip"]
    for sv, origins in blk["origins"].items():
        key = meta.get((chrom, int(sv)))
        if key is None:
            continue
        union = set()
        for o in origins:
            union |= bits(o[2], nt)
        af_union = len(union) / float(nt)
        par = by_coord[key].get("original_decomposed")
        if not par:
            n_no_parent += 1
            continue
        n_parent += 1
        af_par = G("svfreq", par[0])
        if af_par is None:
            continue
        gaps.append((abs(af_union - af_par), chrom, sv, af_union, af_par))
        if abs(af_union - af_par) < 1e-9:
            exact += 1

gaps.sort(reverse=True)
print("multi-origin rows with an original_decomposed parent in the catalog: %d" % n_parent)
print("rows with no parent row: %d" % n_no_parent)
print()
if gaps:
    print("union-of-origins frequency against the parent SV's frequency:")
    print("  identical: %d of %d" % (exact, len(gaps)))
    q = [g[0] for g in gaps]
    q.sort()
    print("  gap  median %.4f   q3 %.4f   worst %.4f" % (q[len(q) // 2], q[3 * len(q) // 4], q[-1]))
    print()
    print("  the five widest:")
    print("    %-7s %-6s %8s %8s %8s" % ("chrom", "sv", "union", "parent", "gap"))
    for g, c, s, a, b in gaps[:5]:
        print("    %-7s %-6s %8.4f %8.4f %8.4f" % (c, s, a, b, g))
