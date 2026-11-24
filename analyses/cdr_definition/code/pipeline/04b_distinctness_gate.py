#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""POSITIVE CONTROL: prove the two arms are actually different data."""
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition/output")
KEY = ["sample", "hap_local", "chm13_chrom"]

a = pd.read_csv(U + "/compat/compat_dip_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
b = pd.read_csv(U + "/compat/compat_gao_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
m = a.merge(b, on=KEY, suffixes=("_dip", "_gao"))
assert len(m) == len(a) == len(b), "arms cover different row sets"

print("rows compared: %d" % len(m))
print()
print("MUST DIFFER (merging changes these by construction):")
print("%-20s %12s %12s %14s %10s" % ("column", "dip median", "gao median", "rows differing", "verdict"))
print("-" * 74)
fails = []
for c, why in [("n_cdrs", "4 dips collapse to ~1 CDR"),
               ("median_cdr_width", "a merged span is far wider than a dip"),
               ("total_cdr_bp", "the span swallows the gaps between dips"),
               ("cdr_fraction", "hypomethylated fraction is computed on that span"),
               ("cdrs_per_mb", "density follows the count")]:
    nd = int((m[c + "_dip"] != m[c + "_gao"]).sum())
    ok = nd > 0.5 * len(m)
    print("%-20s %12.5g %12.5g %14d %10s"
          % (c, m[c + "_dip"].median(), m[c + "_gao"].median(), nd, "OK" if ok else "FAIL"))
    if not ok:
        fails.append((c, why, nd))

nd_pos = int((m.cdr_pos_offset_dip != m.cdr_pos_offset_gao).sum())
print()
print("cdr_pos_offset differs on %d of %d rows (%.1f%%) -- it SHOULD move a little, since a"
      % (nd_pos, len(m), 100.0 * nd_pos / len(m)))
print("merged span's midpoint is not the median of the dip midpoints inside it.")
if nd_pos == 0:
    fails.append(("cdr_pos_offset", "identical on every row -- the arms are the same data", 0))

import glob
import os
nd_calls = nc_calls = 0
for d in sorted(glob.glob(U + "/callsets/gao_cdrs/*")):
    tag = os.path.basename(d)
    src = (CENSV_ROOT + "/analyses/"
           "cdr_calls/output/calls/%s/cdrs.bed" % tag)
    if not os.path.exists(src):
        continue
    nd_calls += sum(1 for _ in open(src))
    nc_calls += sum(1 for _ in open(d + "/cdrs.bed"))
print()
print("call sets on disk: %d dips vs %d Gao CDRs (ratio %.2f)"
      % (nd_calls, nc_calls, nd_calls / float(nc_calls)))
assert nc_calls < 0.5 * nd_calls, \
    "the Gao call set is not materially smaller than the dip set (%d vs %d)" % (nc_calls, nd_calls)

if fails:
    for c, why, nd in fails:
        print("FAIL %s: %s (rows differing: %d)" % (c, why, nd))
    raise SystemExit("DISTINCTNESS GATE FAILED -- the two arms are not distinct data, so any "
                     "agreement between them is an artefact of the pipeline, not a result")
print()
print("[gate] DISTINCTNESS PASSED -- the arms are genuinely different call sets, so their")
print("       agreement on cdr_pos_offset is a finding rather than a copy.")
