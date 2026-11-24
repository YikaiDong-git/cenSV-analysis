#!/usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import importlib.util
import os
import sys
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
U = ROOT + "/analyses/cdr_definition"
MOD = U + "/code/pipeline/03_recompute_phenotypes.py"

spec = importlib.util.spec_from_file_location("rp", MOD)
rp = importlib.util.module_from_spec(spec)
sys.modules["rp"] = rp
spec.loader.exec_module(rp)

gmap = rp.build_globe_map()

owner, multi_contig = {}, set()
for key, by_contig in gmap.items():
    for contig in by_contig:
        owner.setdefault(contig, []).append(key)
        if len(by_contig) > 1:
            multi_contig.add(contig)

bins_by = {}
for tagdir in sorted(glob.glob(rp.CORR + "/*")):
    if os.path.isdir(tagdir):
        for ci, lst in rp.load_bed_by_contig(tagdir + "/binned_freq.bed", 4).items():
            bins_by.setdefault(ci, []).extend(lst)

blocks_by = {}
for key, by_contig in gmap.items():
    for contig, raw in by_contig.items():
        blocks_by[contig] = rp.merge(raw)

total = 0
cls = dict(positioned=0, outside_array=0, multi_contig=0, no_bins=0, contig_not_in_gmap=0)
chrY_pos = 0
for tagdir in sorted(glob.glob(rp.CORR + "/*")):
    if not os.path.isdir(tagdir):
        continue
    for contig, lst in rp.load_bed_by_contig(tagdir + "/cdrs.bed").items():
        for (cs, ce, _x) in lst:
            total += 1
            if contig not in owner:
                cls["contig_not_in_gmap"] += 1
                continue
            if contig in multi_contig:
                cls["multi_contig"] += 1
                continue
            if contig not in bins_by:
                cls["no_bins"] += 1
                continue
            inside = any(ce > a and cs < b for a, b in blocks_by[contig])
            if inside:
                cls["positioned"] += 1
                if owner[contig][0][2] == "chrY":
                    chrY_pos += 1
            else:
                cls["outside_array"] += 1

print("dip calls in cdrs.bed                          : %d" % total)
for k in ("positioned", "outside_array", "multi_contig", "no_bins", "contig_not_in_gmap"):
    print("  %-42s %d" % (k, cls[k]))
s = sum(cls.values())
print("  %-42s %d" % ("sum", s))
if s != total:
    raise SystemExit("FAIL the classes do not partition the call set: %d against %d, "
                     "shortfall %d" % (s, total, total - s))
print("CLOSES")
print("positioned calls on chrY                       : %d" % chrY_pos)
noY = cls["positioned"] - chrY_pos
print("positioned calls on chromosomes 1-22 and X     : %d" % noY)
if noY != 21515:
    raise SystemExit("FAIL the table caption states 21,515; this run gives %d" % noY)
