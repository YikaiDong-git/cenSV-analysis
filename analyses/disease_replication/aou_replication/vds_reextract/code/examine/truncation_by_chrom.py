#!/usr/bin/env python3
from collections import defaultdict
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RX = os.path.dirname(os.path.dirname(HERE))
JOSH = "/home/dataproc/workspace/working/jwein22/centromere_phewas/results/pericentromeric_snp_intervals/pericentromeric_intervals_hail.tsv"
LOCI = os.path.join(RX, "predictor_loci_all.tsv")

josh = {}
for l in open(JOSH):
    f = l.split()
    if not f or f[0] == "chrom" or len(f) < 3:
        continue
    josh[f[0]] = (int(f[1]), int(f[2]))

pos = defaultdict(list)
for l in open(LOCI):
    f = l.rstrip().split("\t")
    if f[0] == "chrom" or not f[0].startswith("chr"):
        continue
    pos[f[0]].append(int(f[1]))

print("%-6s %7s %11s %11s %6s %6s %6s  %-18s" %
      ("chrom", "n_pred", "josh_start", "josh_end", "in", "OUT", "%out", "outside_span_Mb"))
tot_out = 0
aff = []
for c in sorted(pos, key=lambda x: int(x[3:]) if x[3:].isdigit() else 99):
    lo, hi = josh.get(c, (None, None))
    ps = sorted(pos[c])
    if lo is None:
        print("%-6s %7d  (no Josh interval!)" % (c, len(ps)))
        continue
    inside = [p for p in ps if lo <= p <= hi]
    out = [p for p in ps if p < lo or p > hi]
    tot_out += len(out)
    span = ("%.1f-%.1f" % (min(out) / 1e6, max(out) / 1e6)) if out else "-"
    flag = "  <=AFFECTED" if out else ""
    if out:
        aff.append(c)
    print("%-6s %7d %11d %11d %6d %6d %5.1f%%  %-18s%s"
          % (c, len(ps), lo, hi, len(inside), len(out), 100 * len(out) / len(ps), span, flag))
print("\nTOTAL predictor loci outside Josh interval: %d / %d" % (tot_out, sum(len(v) for v in pos.values())))
print("AFFECTED chromosomes (%d): %s" % (len(aff), " ".join(aff)))
print("CLEAN chromosomes: %s" % " ".join(c for c in pos if c not in aff))
