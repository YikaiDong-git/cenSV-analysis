#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Two questions at once."""
import base64
import io
import json
import os
import subprocess
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
WEB = os.path.join(P1, "cenSV_web")
DATA = os.path.join(WEB, "web/public/data")
BEFORE = "da38889"

idx = json.load(io.open(os.path.join(DATA, "hap_index.json"), encoding="utf-8"))

print("=== 1. is the ribbon's carrier row the same as before the change? ===")
same = diff = 0
for chrom in sorted(idx["multi"]):
    now = json.load(io.open(os.path.join(DATA, "hap_%s.json" % chrom), encoding="utf-8"))
    raw = subprocess.run(["git", "-C", WEB, "show",
                          "%s:web/public/data/hap_%s.json" % (BEFORE, chrom)],
                         stdout=subprocess.PIPE, check=True).stdout
    old = json.loads(raw.decode("utf-8"))
    a = {m[0]: m for m in old["mono"]}
    b = {m[0]: m for m in now["mono"]}
    if set(a) != set(b):
        print("  %s: the mono rows themselves changed" % chrom)
        continue
    for k in a:
        if a[k] == b[k]:
            same += 1
        else:
            diff += 1
print("  mono rows identical before and after: %d; changed: %d" % (same, diff))
print("  (mono[4] is the packed row the ribbon paints, so a 0 here means the ribbon is")
print("   drawing exactly what it drew before this change)")
print()

print("=== 2. what orders origin 1 before origin 2? ===")
up = dn = tie = 0
r2_desc = r2_asc = 0
car_desc = car_asc = 0
examples = []


def bits(b64, n):
    raw = base64.b64decode(b64)
    return {k for k in range(n) if raw[k >> 3] >> (7 - (k & 7)) & 1}


def tips_under_y(tree, branch):
    node = tree["branch_node"][branch]
    kids = [[] for _ in range(tree["n_node"])]
    for i, p in enumerate(tree["parent"]):
        if p >= 0:
            kids[p].append(i)
    tip_of = {nd: k for k, nd in enumerate(tree["tip_node"])}
    out, stack = [], [node]
    while stack:
        i = stack.pop()
        if kids[i]:
            stack.extend(kids[i])
        elif i in tip_of:
            out.append(tip_of[i])
    return sum(out) / float(len(out)) if out else None


for chrom in sorted(idx["multi"]):
    blk = json.load(io.open(os.path.join(DATA, "hap_%s.json" % chrom), encoding="utf-8"))
    tree, nt = blk["tree"], blk["tree"]["n_tip"]
    for sv, origins in blk["origins"].items():
        if len(origins) < 2:
            continue
        o1, o2 = origins[0], origins[1]
        y1, y2 = tips_under_y(tree, o1[0]), tips_under_y(tree, o2[0])
        if y1 is None or y2 is None:
            continue
        if y1 < y2:
            up += 1
        elif y1 > y2:
            dn += 1
            if len(examples) < 4:
                examples.append("%s sv %s: origin 1 %s sits BELOW origin 2 %s"
                                % (chrom, sv, o1[0], o2[0]))
        else:
            tie += 1
        if o1[1] is not None and o2[1] is not None:
            if o1[1] > o2[1]:
                r2_desc += 1
            elif o1[1] < o2[1]:
                r2_asc += 1
        c1, c2 = len(bits(o1[2], nt)), len(bits(o2[2], nt))
        if c1 > c2:
            car_desc += 1
        elif c1 < c2:
            car_asc += 1

n = up + dn + tie
print("  by POSITION   : origin 1 higher in %d of %d, lower in %d" % (up, n, dn))
for e in examples:
    print("     %s" % e)
print("  by DECOMP r2  : origin 1 higher in %d, lower in %d" % (r2_desc, r2_asc))
print("  by CARRIERS   : origin 1 larger in %d, smaller in %d" % (car_desc, car_asc))
print()
print("  The payload keeps the resolution's own `rank`, which is the ITERATION that found the")
print("  component: iteration 0 scores the full carrier vector, iteration 1 scores what is")
print("  left after the first component's carriers are masked out.")
