#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Does the cv2>=0.1 candidate set actually contain any [0.1,0.3) predictor?"""
import csv
from collections import Counter
import os as _censv_os
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]

M = CENSV_SCRATCH + "/ssv_joint/manifest_joint_0p1_FULL.tsv"
S = CENSV_SCRATCH + "/ssv_joint/aggregate_0p1/survivors.tsv"

rows = list(csv.DictReader(open(M), delimiter="\t"))
band = Counter("[0.1,0.3)" if float(r["cv"]) < 0.3 else ">=0.3" for r in rows)
print("cv2>=0.1 run, candidate cells by predictor cv2 band: %s" % dict(band))
print("  total %d cells, %d distinct SVs" % (len(rows), len(set(r["sv_locus"] for r in rows))))
low = sorted(set(r["sv_locus"] for r in rows if float(r["cv"]) < 0.3))
print("  distinct SVs in [0.1,0.3): %d %s" % (len(low), low[:6]))
if rows:
    cvs = sorted(float(r["cv"]) for r in rows)
    print("  cv2 of candidates: min %.4f  1st pct %.4f  median %.4f"
          % (cvs[0], cvs[len(cvs) // 100], cvs[len(cvs) // 2]))

sc = Counter()
with open(S) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("survivor_gc") != "True":
            continue
        sc["[0.1,0.3)" if float(r["cv_r2"]) < 0.3 else ">=0.3"] += 1
print("scan-stage survivor_gc cells by band: %s" % dict(sc))
