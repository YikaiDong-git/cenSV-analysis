#!/usr/bin/env python3
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RX = os.path.dirname(os.path.dirname(HERE))
JOSH = "/home/dataproc/workspace/working/jwein22/centromere_phewas/results/pericentromeric_snp_intervals/pericentromeric_intervals_hail.tsv"
BIMDIR = "/home/dataproc/workspace/working/jwein22/centromere_phewas/results/vds_pericentromeric_snps_min_ac1000"
LOCI = os.path.join(RX, "predictor_loci_all.tsv")

josh = {}
for l in open(JOSH):
    f = l.split()
    if not f or f[0] == "chrom" or len(f) < 3:
        continue
    josh[f[0]] = (int(f[1]), int(f[2]))

pred = defaultdict(list)
for l in open(LOCI):
    f = l.rstrip().split("\t")
    if f[0] == "chrom" or not f[0].startswith("chr"):
        continue
    pred[f[0]].append((int(f[1]), f[2], f[3]))

print("%-6s %8s %9s %11s %11s" % ("chrom", "in_ivl", "inBimPos", "DROP_pos", "DROP_allele"))
tot_in = tot_droppos = tot_dropall = 0
per_chrom_drop = {}
for c in sorted(pred, key=lambda x: int(x[3:]) if x[3:].isdigit() else 99):
    lo, hi = josh.get(c, (None, None))
    if lo is None:
        continue
    inside = [(p, r, a) for (p, r, a) in pred[c] if lo <= p <= hi]
    bim = "%s/%s/%s.pericentromeric_snps.min_ac1000.bim" % (BIMDIR, c, c)
    if not os.path.exists(bim):
        print("%-6s  (no bim: %s)" % (c, bim))
        continue
    bimpos = set()
    bimtrip = set()
    for l in open(bim):
        w = l.split()
        p = int(w[3])
        bimpos.add(p)
        bimtrip.add((p, w[5], w[4]))
        bimtrip.add((p, w[4], w[5]))
    droppos = [t for t in inside if t[0] not in bimpos]
    dropall = [t for t in inside if (t[0], t[1], t[2]) not in bimtrip]
    tot_in += len(inside)
    tot_droppos += len(droppos)
    tot_dropall += len(dropall)
    per_chrom_drop[c] = len(droppos)
    print("%-6s %8d %9d %11d %11d"
          % (c, len(inside), len(inside) - len(droppos), len(droppos), len(dropall)))
print("\nTOTAL inside-interval predictor loci: %d" % tot_in)
print("DROPPED (position absent) = min_ac + repeat: %d (%.1f%%)" % (tot_droppos, 100 * tot_droppos / tot_in))
print("DROPPED (allele-exact absent): %d (%.1f%%)" % (tot_dropall, 100 * tot_dropall / tot_in))
print("\nper-chrom drop:", dict(sorted(per_chrom_drop.items(), key=lambda kv: -kv[1])))
