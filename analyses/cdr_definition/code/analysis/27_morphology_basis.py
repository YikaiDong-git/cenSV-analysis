#!/usr/bin/env python3
from __future__ import print_function

import gzip
import os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

R = CENSV_ROOT + "/resources"
FAI = os.path.join(R, "chm13v2.0.fa.gz.fai")
ACTIVE = os.path.join(R, "chm13v2.0_censat_v2.1.active.bed")
CENSAT = os.path.join(R, "chm13v2.0_censat_v2.1.bed")
ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

sizes = {}
with open(FAI) as fh:
    for line in fh:
        f = line.split("\t")
        sizes[f[0]] = int(f[1])


def span(path):
    lo, hi = {}, {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path) as fh:
        for line in fh:
            if isinstance(line, bytes):
                line = line.decode()
            if not line.strip() or line.startswith(("track", "#")):
                continue
            f = line.split("\t")
            c, a, b = f[0], int(f[1]), int(f[2])
            lo[c] = min(lo.get(c, a), a)
            hi[c] = max(hi.get(c, b), b)
    return {c: (lo[c], hi[c]) for c in lo}


def klass(r):
    if r < 1.7:
        return "metacentric"
    if r < 3.0:
        return "submetacentric"
    if r < 7.0:
        return "subtelocentric"
    return "acrocentric"


act, cen = span(ACTIVE), span(CENSAT)

print("%-7s %10s %10s %10s %8s %7s %-16s | %-16s %s"
      % ("chrom", "length", "p arm", "q arm", "CI", "q/p", "class (active)",
         "class (censat)", "agree"))
print("-" * 118)
rows = []
for c in ORDER:
    if c not in act or c not in sizes:
        continue
    L = sizes[c]
    ma = (act[c][0] + act[c][1]) / 2.0
    mc = (cen[c][0] + cen[c][1]) / 2.0
    ra, rc = (L - ma) / ma, (L - mc) / mc
    ka, kc = klass(ra), klass(rc)
    ci = 100.0 * ma / L
    rows.append((c, ka, kc, ra, ci))
    print("%-7s %10d %10d %10d %7.1f%% %7.2f %-16s | %-16s %s"
          % (c, L, int(ma), int(L - ma), ci, ra, ka, kc,
             "" if ka == kc else "DIFFERS"))

print()
for k in ("metacentric", "submetacentric", "subtelocentric", "acrocentric"):
    m = [r[0] for r in rows if r[1] == k]
    if m:
        print("%-16s %2d: %s" % (k, len(m), ", ".join(m)))

dis = [r[0] for r in rows if r[1] != r[2]]
print("\nproxies disagree on: %s" % (", ".join(dis) if dis else "none"))

FEATURE_TABLE = (CENSV_ROOT + "/analyses/"
                 "imputability/code/01_build_feature_table.py")
if 'ACROCENTRIC = {"chr13", "chr14", "chr15", "chr21", "chr22"}' not in open(FEATURE_TABLE).read():
    raise SystemExit("FAIL the imputability feature table no longer defines the same "
                     "acrocentric set that Supplementary Table tab:cdr_side prints")

TRAD = {"chr13", "chr14", "chr15", "chr21", "chr22"}
got = set(r[0] for r in rows if r[1] == "acrocentric")
print("computed acrocentric set matches the conventional five: %s" % (got == TRAD))
if got != TRAD:
    print("  computed but not conventional: %s" % sorted(got - TRAD))
    print("  conventional but not computed: %s" % sorted(TRAD - got))
    print()
    print("This is the expected outcome and the reason tab:cdr_side prints the conventional")
    print("classification: CHM13v2.0 assembles the acrocentric short arms in full, so measured")
    print("arm ratios do not reproduce cut-offs calibrated on metaphase chromosomes.")
