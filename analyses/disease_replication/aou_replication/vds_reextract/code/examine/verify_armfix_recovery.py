#!/usr/bin/env python3
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RX = os.path.dirname(os.path.dirname(HERE))
LOCI = os.path.join(RX, "predictor_loci_all.tsv")
ROOT = "/home/dataproc/workspace/working/jwein22/centromere_phewas/results/vds_pericentromeric_snps_min_ac1000"
PRE = ROOT + "_pre_armfix_20260711"

pred = defaultdict(list)
for l in open(LOCI):
    f = l.rstrip().split("\t")
    if f[0] == "chrom" or not f[0].startswith("chr"):
        continue
    pred[f[0]].append((int(f[1]), f[2], f[3]))

JOSH_CLAIM = {"chr1": (1330, 1537), "chr9": (350, 415), "chr16": (1041, 1110)}


def bim_sets(path):
    pos, trip = set(), set()
    if not os.path.exists(path):
        return None, None
    for l in open(path):
        w = l.split()
        p = int(w[3])
        pos.add(p)
        trip.add((p, w[5], w[4]))
        trip.add((p, w[4], w[5]))
    return pos, trip


for c in ["chr1", "chr9", "chr16"]:
    newbim = "%s/%s/%s.pericentromeric_snps.min_ac1000.bim" % (ROOT, c, c)
    oldbim = "%s/%s/%s.pericentromeric_snps.min_ac1000.bim" % (PRE, c, c)
    npos, ntrip = bim_sets(newbim)
    opos, otrip = bim_sets(oldbim)
    ps = pred[c]
    tot = len(ps)
    rec_pos = sum(1 for (p, r, a) in ps if p in npos)
    rec_all = sum(1 for (p, r, a) in ps if (p, r, a) in ntrip)
    old_pos = sum(1 for (p, r, a) in ps if opos and p in opos) if opos else None
    jc, jt = JOSH_CLAIM[c]
    print("%s: our predictors=%d" % (c, tot))
    print("   recovered by POSITION in promoted bim: %d (%.2f%%)" % (rec_pos, 100 * rec_pos / tot))
    print("   recovered ALLELE-EXACT:                %d (%.2f%%)" % (rec_all, 100 * rec_all / tot))
    print("   Josh Section-3 claim:                  %d (%.2f%%)  -> match: %s"
          % (jc, 100 * jc / jt, "YES" if abs(rec_pos - jc) <= 3 or abs(rec_all - jc) <= 3 else "CHECK"))
    if old_pos is not None:
        print("   [before arm-fix, old bim]: %d (%.2f%%)  -> gain +%d loci"
              % (old_pos, 100 * old_pos / tot, rec_pos - old_pos))
    print()
