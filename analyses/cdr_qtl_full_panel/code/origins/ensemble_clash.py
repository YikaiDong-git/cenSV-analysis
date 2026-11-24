#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Which chromosomes carry the sv_index clash that stops the catalog build?"""
import collections
import csv
import io
import os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
IMP = os.path.join(P1, "results/imputation/decompose/cosine_nj")
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

rows = collections.defaultdict(dict)
stamp = {}
for ens in ("best8", "best6"):
    for c in CHROMS:
        f = os.path.join(IMP, "%s_%s" % (ens, SUF), c, "%s_imputation_summary.tsv" % c)
        if not os.path.exists(f):
            continue
        stamp[(ens, c)] = os.path.getmtime(f)
        with io.open(f, encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                key = (c, r["sv_index"])
                rows[key][ens] = (r["pos"], r["svtype"], r["svlen"])

per = collections.Counter()
both = collections.Counter()
for (c, _i), d in rows.items():
    if len(d) < 2:
        continue
    both[c] += 1
    if len(set(d.values())) > 1:
        per[c] += 1

print("%-7s %10s %8s   %s" % ("chrom", "in both", "clashing", "best8 summary written"))
tot_b = tot_c = 0
for c in CHROMS:
    if not both[c]:
        continue
    tot_b += both[c]
    tot_c += per[c]
    when = stamp.get(("best8", c))
    import time
    w = time.strftime("%Y-%m-%d %H:%M", time.localtime(when)) if when else "-"
    mark = "  <-- clash" if per[c] else ""
    print("%-7s %10d %8d   %s%s" % (c, both[c], per[c], w, mark))
print("%-7s %10d %8d" % ("TOTAL", tot_b, tot_c))
