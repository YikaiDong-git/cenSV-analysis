#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""If the carrier ribbon were painted from the SV's own pre-decomposition row, what changes?"""
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

PAIR = [0, 2, 1, 3]


def bits(b64, n):
    raw = base64.b64decode(b64)
    return {k for k in range(n) if raw[k >> 3] >> (7 - (k & 7)) & 1}


def demo_row(b64, k, ncol):
    raw = base64.b64decode(b64)
    out = []
    for c in range(ncol):
        cell = k * ncol + c
        byte = raw[cell >> 2]
        hi = (byte >> (7 - ((cell & 3) << 1))) & 1
        lo = (byte >> (6 - ((cell & 3) << 1))) & 1
        out.append(hi + 2 * lo)
    return out


def ldjs_row(b64, n):
    raw = base64.b64decode(b64)
    return [PAIR[(raw[c >> 2] >> (6 - ((c & 3) << 1))) & 3] for c in range(n)]


idx = json.load(io.open(os.path.join(DATA, "hap_index.json"), encoding="utf-8"))
cat = json.load(io.open(os.path.join(DATA, "catalog.json"), encoding="utf-8"))
ci = {k: i for i, k in enumerate(cat["columns"])}
G = lambda k, i: cat["data"][ci[k]][i]
coord = {(G("chrom", i), G("sv_index", i)):
         (G("pos", i), G("svtype", i), G("svlen", i)) for i in range(cat["rows"])}

demo = json.load(io.open(DEMO, encoding="utf-8"))
snames = [r["s"] for r in demo["samples"]]
six = {s: i for i, s in enumerate(snames)}

n_origins = collections.Counter()
rows_parent = rows_noparent = 0
later_car = later_marked_parent = later_marked_now = 0
blank_parent = blank_now = 0
n_later = n_later_p = later_car_p = 0
unassigned = par_car = 0
contra = []
miss = []
af_gap = []

for chrom in sorted(idx["multi"]):
    blk = demo["chroms"][chrom]
    hb = json.load(io.open(os.path.join(DATA, "hap_%s.json" % chrom), encoding="utf-8"))
    tips = hb["tree"]["tips"]
    nt = hb["tree"]["n_tip"]
    mono = {m[0]: m for m in hb["mono"]}
    sv = blk["sv"]
    at = collections.defaultdict(dict)
    for k in range(blk["n_sv"]):
        at[(sv["pos"][k], sv["t"][k], sv["len"][k])].setdefault(sv["cat"][k], []).append(k)

    for s, origins in sorted(hb["origins"].items(), key=lambda kv: int(kv[0])):
        c = coord.get((chrom, int(s)))
        if c is None:
            continue
        n_origins[len(origins)] += 1
        grp = at.get(c, {})
        pk = grp["original_decomposed"][0] if "original_decomposed" in grp else None
        if pk is None:
            rows_noparent += 1
        else:
            rows_parent += 1
            par = demo_row(blk["gt"], pk, len(snames))
            called = [v for v in par if v != 3]
            af_gap.append((abs(sum(called) / (2.0 * len(called)) - sv["af"][pk]), chrom, s))
        now = ldjs_row(mono[int(s)][4], len(snames))

        for i, o in enumerate(origins):
            car = bits(o[2], nt)
            if i == 0:
                continue
            n_later += 1
            hit_p = hit_n = 0
            for k in car:
                j = six.get(tips[k].rsplit("_hap", 1)[0])
                if j is None:
                    continue
                later_car += 1
                if pk is not None:
                    later_car_p += 1
                if now[j] in (1, 2):
                    hit_n += 1
                if pk is not None:
                    if par[j] in (1, 2):
                        hit_p += 1
                    elif par[j] == 0:
                        contra.append((chrom, s, o[0], tips[k]))
                    else:
                        miss.append((chrom, s, o[0], tips[k]))
            later_marked_now += hit_n
            later_marked_parent += hit_p
            if not hit_n:
                blank_now += 1
            if pk is not None:
                if not hit_p:
                    blank_parent += 1
                n_later_p += 1

        if pk is not None:
            owned = set()
            for o in origins:
                for k in bits(o[2], nt):
                    j = six.get(tips[k].rsplit("_hap", 1)[0])
                    if j is not None:
                        owned.add(j)
            un = sum(1 for j in range(len(snames)) if par[j] in (1, 2) and j not in owned)
            unassigned += un
            par_car += sum(1 for j in range(len(snames)) if par[j] in (1, 2))

print("multi-origin rows drawn: %d  (parent present %d, absent %d)"
      % (rows_parent + rows_noparent, rows_parent, rows_noparent))
print("origins per row: %s"
      % ", ".join("%d origins x%d" % (k, v) for k, v in sorted(n_origins.items())))
print()
print("later-origin carrier haplotypes: %d over %d later origins" % (later_car, n_later))
print("  marked by the SHIPPED ribbon (the row's own component): %d (%.1f%%), %d origins blank"
      % (later_marked_now, 100.0 * later_marked_now / later_car, blank_now))
print()
print("on the %d rows that HAVE a parent: %d later origins, %d carrier haplotypes"
      % (rows_parent, n_later_p, later_car_p))
print("  marked by the PARENT row: %d (%.1f%%), %d origins blank"
      % (later_marked_parent, 100.0 * later_marked_parent / later_car_p, blank_parent))
print("  the remaining %d carrier haplotypes sit on the %d rows with no parent, where the"
      % (later_car - later_car_p, rows_noparent))
print("  ribbon keeps the shipped reading.")
print()
print("samples the parent row calls carriers: %d; of them %d (%.1f%%) belong to no drawn"
      % (par_car, unassigned, 100.0 * unassigned / par_car))
print("  origin - the decomposition's residual, which the origin column leaves uncoloured.")
print()
print("carrier haplotypes the parent row CONTRADICTS (its sample reads dosage 0): %d"
      % len(contra))
for r in contra[:6]:
    print("    %s sv_index %s origin %s tip %s" % r)
print("carrier haplotypes whose sample the parent row leaves uncalled: %d" % len(miss))
for r in miss[:6]:
    print("    %s sv_index %s origin %s tip %s" % r)
print()
if af_gap:
    af_gap.sort()
    print("parent dosage row vs the demo's own af for that row: worst %.6f (%s sv_index %s)"
          % (af_gap[-1][0], af_gap[-1][1], af_gap[-1][2]))
