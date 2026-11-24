#!/usr/bin/env python3
from collections import defaultdict

JOSH = "/home/dataproc/workspace/working/jwein22/centromere_phewas/results/pericentromeric_snp_intervals/pericentromeric_intervals_hail.tsv"
WDIR = "/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/phewide_all_handoff"

josh = {}
for l in open(JOSH):
    f = l.split()
    if not f or f[0] == "chrom" or len(f) < 3:
        continue
    josh[f[0]] = (int(f[1]), int(f[2]))

HITS = {"chr1_sv650", "chr16_sv559", "chr16_sv37", "chr16_sv1029", "chr1_sv1029"}


def svkey(name):
    p = name.split("_")
    return p[0] + "_" + p[1]


for c in ["chr1", "chr9", "chr16"]:
    lo, hi = josh[c]
    tot = defaultdict(float)
    outw = defaultdict(float)
    totn = defaultdict(int)
    outn = defaultdict(int)
    for l in open("%s/sv_weights_long_%s.tsv" % (WDIR, c)):
        f = l.rstrip("\n").split("\t")
        if f[0] == "sv_locus" or len(f) < 7:
            continue
        sv = svkey(f[0])
        p = int(f[2])
        try:
            w = abs(float(f[6]))
        except ValueError:
            continue
        tot[sv] += w
        totn[sv] += 1
        if p < lo or p > hi:
            outw[sv] += w
            outn[sv] += 1
    fracs = {sv: (outw[sv] / tot[sv] if tot[sv] > 0 else 0) for sv in tot}
    g20 = [sv for sv in tot if fracs[sv] > 0.20]
    g50 = [sv for sv in tot if fracs[sv] > 0.50]
    anyout = [sv for sv in tot if outn[sv] > 0]
    print("==== %s  (Josh interval %d-%d) ====" % (c, lo, hi))
    print("  SVs on chrom: %d | with ANY outside locus: %d | >20%% weight lost: %d | >50%%: %d"
          % (len(tot), len(anyout), len(g20), len(g50)))
    for sv in sorted(tot):
        if sv in HITS:
            print("  HIT %-14s weight_out=%5.1f%%  (%d/%d loci outside)"
                  % (sv, 100 * fracs[sv], outn[sv], totn[sv]))
    top = sorted(tot, key=lambda s: -fracs[s])[:5]
    print("  top-truncated: " + ", ".join("%s %.0f%%" % (s, 100 * fracs[s]) for s in top))
    print()
