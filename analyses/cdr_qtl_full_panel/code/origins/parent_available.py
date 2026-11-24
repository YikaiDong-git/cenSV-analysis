#!/usr/bin/env python
# -*- coding: utf-8 -*-
import base64
import collections
import io
import json
import os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
DATA = os.path.join(P1, "cenSV_web/web/public/data")
DEMO = os.path.join(P1, "discovery/2026-08-30_censv_individual_views/output/demo",
                    "individual_demo.json")


def bits(b64, n):
    raw = base64.b64decode(b64)
    return {k for k in range(n) if raw[k >> 3] >> (7 - (k & 7)) & 1}


idx = json.load(io.open(os.path.join(DATA, "hap_index.json"), encoding="utf-8"))
cat = json.load(io.open(os.path.join(DATA, "catalog.json"), encoding="utf-8"))
ci = {k: i for i, k in enumerate(cat["columns"])}
G = lambda k, i: cat["data"][ci[k]][i]
coord = {(G("chrom", i), G("sv_index", i)):
         (G("pos", i), G("svtype", i), G("svlen", i)) for i in range(cat["rows"])}

demo = json.load(io.open(DEMO, encoding="utf-8"))

have = missing = 0
res_at = 0
gap_rows = []
for chrom in sorted(idx["multi"]):
    blk = demo["chroms"].get(chrom)
    hb = json.load(io.open(os.path.join(DATA, "hap_%s.json" % chrom), encoding="utf-8"))
    nt = hb["tree"]["n_tip"]
    sv = blk["sv"]
    at = collections.defaultdict(dict)
    for k in range(blk["n_sv"]):
        at[(sv["pos"][k], sv["t"][k], sv["len"][k])].setdefault(sv["cat"][k], []).append(k)
    for s, origins in hb["origins"].items():
        c = coord.get((chrom, int(s)))
        if c is None:
            continue
        grp = at.get(c, {})
        if "residual" in grp:
            res_at += 1
        if "original_decomposed" not in grp:
            missing += 1
            continue
        have += 1
        union = set()
        for o in origins:
            union |= bits(o[2], nt)
        pk = grp["original_decomposed"][0]
        gap_rows.append((sv["af"][pk], len(union) / float(nt), chrom, s))

print("multi-origin rows drawn: %d" % (have + missing))
print("  with the pre-decomposition parent in the same block: %d" % have)
print("  without one: %d" % missing)
print("  with a `residual` component at the same coordinate: %d" % res_at)
print()
if gap_rows:
    d = sorted(abs(a - b) for a, b, _, _ in gap_rows)
    print("parent allele frequency vs the union of drawn origins:")
    print("  union below parent in %d of %d rows"
          % (sum(1 for a, b, _, _ in gap_rows if b < a - 1e-12), len(gap_rows)))
    print("  gap  median %.4f   q3 %.4f   worst %.4f"
          % (d[len(d) // 2], d[3 * len(d) // 4], d[-1]))
    print()
    print("  So a ribbon painted from the union shows the carriers an origin EXPLAINS,")
    print("  and the parent row is the one that shows every carrier the SV has.")
