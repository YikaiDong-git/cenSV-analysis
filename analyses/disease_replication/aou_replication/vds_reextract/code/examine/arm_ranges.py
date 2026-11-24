#!/usr/bin/env python3
from collections import defaultdict
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RX = os.path.dirname(os.path.dirname(HERE))
JOSH = "/home/dataproc/workspace/working/jwein22/centromere_phewas/results/pericentromeric_snp_intervals/pericentromeric_intervals_hail.tsv"
LOCI = os.path.join(RX, "predictor_loci_all.tsv")
PAD = 10000

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

rows = []
for c in ["chr1", "chr9", "chr16"]:
    lo, hi = josh[c]
    ps = sorted(pos[c])
    below = [p for p in ps if p < lo]
    above = [p for p in ps if p > hi]
    print("%s  Josh interval = %d - %d" % (c, lo, hi))
    for side, arr in [("PROXIMAL (below start)", below), ("DISTAL (above end)", above)]:
        if arr:
            s, e = min(arr) - PAD, max(arr) + PAD
            print("   %-22s %d loci  ->  %s:%d-%d  (%.3f-%.3fMb, width %.2fMb)"
                  % (side, len(arr), c, s, e, s / 1e6, e / 1e6, (e - s) / 1e6))
            rows.append((c, s, e, len(arr)))
    print()

print("=== intervals TSV rows (chrom  start_1based  end_1based_inclusive) ===")
for c, s, e, n in rows:
    print("%s\t%d\t%d" % (c, s, e))
print("\nTOTAL: %d loci, %d intervals, %.2f Mb" %
      (sum(n for *_, n in rows), len(rows), sum((e - s) for _, s, e, _ in rows) / 1e6))
