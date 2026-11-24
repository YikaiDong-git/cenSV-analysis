#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Mean 5mC along the active array, per chromosome -- the signal the CDRs are CALLED FROM."""
import glob
import importlib.util
import os
import sys
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
UNIT = ROOT + "/analyses/cdr_phenotypes"
MOD = UNIT + "/code/pipeline/compute_phenotypes.py"
OUT = UNIT + "/output/phenotypes"

spec = importlib.util.spec_from_file_location("rp2", MOD)
rp2 = importlib.util.module_from_spec(spec)
sys.modules["rp2"] = rp2
spec.loader.exec_module(rp2)

NB_REL, NB_ABS = 60, 60
ABS_MAX = 6.0
ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX"]


def _num(x):
    try:
        float(x)
        return True
    except Exception:
        return False


gmap = rp2.build_globe_map()
bins_by = {}
for tagdir in sorted(glob.glob(rp2.CORR + "/*")):
    if not os.path.isdir(tagdir):
        continue
    for ci, lst in rp2.load_bed_by_contig(tagdir + "/binned_freq.bed", 4).items():
        bins_by.setdefault(ci, []).extend(lst)
print("contigs with window data: %d" % len(bins_by))

acc = {c: dict(rel_s=np.zeros(NB_REL), rel_w=np.zeros(NB_REL), rel_n=np.zeros(NB_REL),
               abs_s=np.zeros(NB_ABS), abs_w=np.zeros(NB_ABS), abs_n=np.zeros(NB_ABS),
               nhap=0) for c in ORDER}
nwin = 0
for (sample, hap, chrom), by_contig in sorted(gmap.items()):
    if chrom not in acc or len(by_contig) > 1:
        continue
    contig, raw = next(iter(by_contig.items()))
    if contig not in bins_by:
        continue
    blocks = rp2.merge(raw)
    aa_s, aa_e = blocks[0][0], blocks[-1][1]
    S = sum(b - a for a, b in blocks)
    if S <= 0:
        continue
    acc[chrom]["nhap"] += 1
    for (ws, we, ex) in bins_by[contig]:
        if not ex or not _num(ex[0]):
            continue
        if not any(we > a and ws < b for a, b in blocks):
            continue
        m = (ws + we) / 2.0
        cov = float(ex[1]) if len(ex) > 1 and _num(ex[1]) else 1.0
        meth = float(ex[0])
        off = rp2.cum_covered(m, blocks)
        rel = off / S
        if not (0.0 <= rel <= 1.0):
            continue
        i = min(int(rel * NB_REL), NB_REL - 1)
        acc[chrom]["rel_s"][i] += meth
        acc[chrom]["rel_w"][i] += meth * cov
        acc[chrom]["rel_n"][i] += cov
        ab = off / 1e6
        if ab < ABS_MAX:
            j = min(int(ab / ABS_MAX * NB_ABS), NB_ABS - 1)
            acc[chrom]["abs_s"][j] += meth
            acc[chrom]["abs_w"][j] += meth * cov
            acc[chrom]["abs_n"][j] += cov
        nwin += 1
print("windows aggregated: %d" % nwin)

rows = []
for c in ORDER:
    a = acc[c]
    for i in range(NB_REL):
        if a["rel_n"][i] > 0:
            rows.append(dict(chrom=c, frame="rel", bin=i,
                             pos=(i + 0.5) / NB_REL,
                             meth_cw=a["rel_w"][i] / a["rel_n"][i],
                             cov=a["rel_n"][i], nhap=a["nhap"]))
    for j in range(NB_ABS):
        if a["abs_n"][j] > 0:
            rows.append(dict(chrom=c, frame="abs", bin=j,
                             pos=(j + 0.5) / NB_ABS * ABS_MAX,
                             meth_cw=a["abs_w"][j] / a["abs_n"][j],
                             cov=a["abs_n"][j], nhap=a["nhap"]))
P = pd.DataFrame(rows)
P.to_csv(OUT + "/methyl_profile.tsv", sep="\t", index=False, float_format="%.6g")
print("profile rows: %d" % len(P))

dp = pd.read_csv(OUT + "/per_dip_positions.tsv", sep="\t")
dp = dp[dp.chm13_chrom.isin(ORDER)]
dp["rel"] = (dp.pos_splice + 1.0) / 2.0
print("")
print("%-7s %14s %16s %10s" % ("chrom", "profile min at", "median dip at", "5mC min"))
agree = 0
tested = 0
for c in ORDER:
    pr = P[(P.chrom == c) & (P.frame == "rel")]
    d = dp[dp.chm13_chrom == c]
    if len(pr) < 20 or len(d) < 50:
        continue
    tested += 1
    lowpos = float(pr.loc[pr.meth_cw.idxmin(), "pos"])
    dipmed = float(d.rel.median())
    ok = abs(lowpos - dipmed) < 0.25
    agree += ok
    print("%-7s %14.3f %16.3f %10.1f  %s"
          % (c, lowpos, dipmed, pr.meth_cw.min(), "" if ok else "  <-- differ"))
print("")
print("profile minimum within 0.25 of the median called dip on %d of %d chromosomes"
      % (agree, tested))
assert agree >= tested * 0.5, (
    "the methylation profile minimum disagrees with the called dips on most chromosomes "
    "(%d of %d) -- the two should broadly coincide, so one of them is being computed wrong"
    % (agree, tested))
print("[check] the continuous profile and the discrete calls broadly agree")
print("[written] %s/methyl_profile.tsv" % OUT)
