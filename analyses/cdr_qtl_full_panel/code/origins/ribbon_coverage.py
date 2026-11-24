#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Does the carrier ribbon mark the carriers of EVERY origin, or only of one?"""
import base64
import io
import json
import os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
DATA = os.path.join(P1, "cenSV_web/web/public/data")
PAIR = [0, 2, 1, 3]


def dosages(b64, n):
    raw = base64.b64decode(b64)
    return [PAIR[(raw[c >> 2] >> (6 - ((c & 3) << 1))) & 3] for c in range(n)]


def bits(b64, n):
    raw = base64.b64decode(b64)
    return {k for k in range(n) if raw[k >> 3] >> (7 - (k & 7)) & 1}


idx = json.load(io.open(os.path.join(DATA, "hap_index.json"), encoding="utf-8"))
si_of = {r["s"]: i for i, r in enumerate(idx["samples"])}
ns = len(idx["samples"])

cat = json.load(io.open(os.path.join(DATA, "catalog.json"), encoding="utf-8"))
ci = {k: i for i, k in enumerate(cat["columns"])}
G = lambda k, i: cat["data"][ci[k]][i]
row_of = {}
for i in range(cat["rows"]):
    row_of[(G("chrom", i), G("sv_index", i))] = (G("mapped_haplogroup", i), G("sv_id", i))

out, tot, marked, dead = [], 0, 0, 0
n_org = 0
for chrom in sorted(idx["multi"]):
    blk = json.load(io.open(os.path.join(DATA, "hap_%s.json" % chrom), encoding="utf-8"))
    nt = blk["tree"]["n_tip"]
    tips = blk["tree"]["tips"]
    packed_of = {m[0]: m[4] for m in blk["mono"]}
    for sv, origins in blk["origins"].items():
        pk = packed_of.get(int(sv))
        if pk is None:
            continue
        dos = dosages(pk, ns)
        hg, sv_id = row_of.get((chrom, int(sv)), ("?", "?"))
        per = []
        for o in origins:
            car = bits(o[2], nt)
            hit = 0
            for k in car:
                j = si_of.get(tips[k].rsplit("_hap", 1)[0])
                if j is not None and dos[j]:
                    hit += 1
            per.append((o[0], len(car), hit, o[3]))
            tot += len(car)
            marked += hit
            n_org += 1
            if hit == 0:
                dead += 1
        out.append((chrom, sv, sv_id, hg, sum(1 for d in dos if d), per))

print("%-6s %-6s %-16s %-7s %-6s  %s" %
      ("chrom", "sv", "sv_id", "row", "dos>0", "branch/carriers/marked in ribbon/tier"))
for chrom, sv, sv_id, hg, nd, per in out[:12]:
    print("%-6s %-6s %-16s %-7s %-6d  %s" % (chrom, sv, sv_id[:16], hg, nd,
          "   ".join("%s/%d/%d/%s" % (b, n, h, "imp" if t else "asm") for b, n, h, t in per)))

print()
print("carrier haplotypes across all drawn origins: %d" % tot)
print("of those, whose SAMPLE the ribbon marks as a carrier: %d (%.1f%%)"
      % (marked, 100.0 * marked / tot))
print("origins with NO carrier marked in the ribbon: %d of %d" % (dead, n_org))
print()
first = [(c, s, p[0]) for c, s, i, h, d, p in out]
later = [(c, s, q) for c, s, i, h, d, p in out for q in p[1:]]
for nm, g in [("first origins", first), ("later origins", later)]:
    t = sum(q[1] for _, _, q in g)
    m = sum(q[2] for _, _, q in g)
    z = sum(1 for _, _, q in g if q[2] == 0)
    print("%-15s n=%-4d carriers %-5d marked %-5d (%5.1f%%)  origins with none marked %d"
          % (nm, len(g), t, m, 100.0 * m / t if t else 0, z))
