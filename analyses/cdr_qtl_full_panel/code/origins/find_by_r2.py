#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Which row is the one on screen: two origins, 12 carriers each, r2 0.361 and 0.682?"""
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
row_of = {(G("chrom", i), G("sv_index", i)): (G("mapped_haplogroup", i), G("sv_id", i))
          for i in range(cat["rows"])}

WANT = {0.361, 0.682}
for chrom in sorted(idx["multi"]):
    blk = json.load(io.open(os.path.join(DATA, "hap_%s.json" % chrom), encoding="utf-8"))
    nt, tips = blk["tree"]["n_tip"], blk["tree"]["tips"]
    packed_of = {m[0]: m[4] for m in blk["mono"]}
    for sv, origins in blk["origins"].items():
        r2s = {round(o[1], 3) for o in origins if o[1] is not None}
        if not WANT <= r2s:
            continue
        pk = packed_of.get(int(sv))
        dos = dosages(pk, ns) if pk else None
        hg, sv_id = row_of.get((chrom, int(sv)), ("?", "?"))
        print("%s sv_index %s  %s   catalog row sits on %s" % (chrom, sv, sv_id, hg))
        for o in origins:
            car = bits(o[2], nt)
            hit = sum(1 for k in car
                      if si_of.get(tips[k].rsplit("_hap", 1)[0]) is not None
                      and dos[si_of[tips[k].rsplit("_hap", 1)[0]]])
            print("   branch %-6s %2d carriers  r2 %.3f  %-13s  marked in the ribbon: %d"
                  % (o[0], len(car), o[1], "imputed" if o[3] else "assembly only", hit))
