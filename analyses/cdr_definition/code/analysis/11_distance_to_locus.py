#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The S15 claim -- \"the dip moves toward the locus\" -- recomputed under both CDR units."""
import numpy as np
import pandas as pd
from scipy.stats import kruskal
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
U = ROOT + "/analyses/cdr_definition"
WU = ROOT + "/analyses/cdr_qtl_full_panel"
LEAD = 44300221
SHIPPED_DIST = [0.571, 0.461, 0.297]
ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
TAB = "\t"

hc = pd.read_csv(U + "/output/phenotypes/per_hap_chrom_phenotypes_reproduced.tsv", sep=TAB)

calls, genomewide = {}, {}
for arm in ("dip", "gao"):
    d = pd.read_csv("%s/output/phenotypes/per_call_positions_%s.tsv" % (U, arm), sep=TAB)
    d = d[d.chm13_chrom.isin(ORDER)].copy()
    d["rel"] = (d.pos_splice + 1.0) / 2.0
    d["len_mb"] = d.active_len_sum / 1e6
    genomewide[arm] = d
    calls[arm] = d[d.chm13_chrom == "chr8"].copy()

ymax = float(np.percentile(genomewide["dip"].len_mb, 99.5))
c = hc[(hc.chm13_chrom == "chr8") & (hc.active_start <= LEAD) & (hc.active_end >= LEAD)].copy()
c["len_mb"] = c.active_len_sum / 1e6
c["rel"] = (LEAD - c.active_start) / (c.active_end - c.active_start)
c = c[c.len_mb <= ymax]
g = c.groupby(pd.qcut(c.len_mb, 5, duplicates="drop"), observed=True).agg(
    y=("len_mb", "median"), x=("rel", "median"))
svx = g.x.to_numpy()
track = float(np.median(svx))
print("locus track: %d haplotypes carry the locus inside the array; "
      "median relative position %.4f" % (len(c), track))
print("ymax (99.5th pct of dip-arm call array length, all chromosomes) = %.3f Mb" % ymax)

vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep=TAB)
loc = vt[(vt.chrom == "chr8") & (vt.pos == LEAD)].sv_locus.iloc[0]
G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
sm = open(WU + "/geno/full_samples.txt").read().split()
k = {l: j for j, l in enumerate(vt.sv_locus)}[loc]
dose = pd.Series(G[k], index=sm, name="dose")

for arm in calls:
    calls[arm] = calls[arm].join(dose, on="sample")
    calls[arm] = calls[arm][calls[arm].dose.notna()]

DOSE = sorted(calls["dip"].dose.unique())
print()
print("Median |relative call position - locus|, by lead-deletion dosage")
print("%-6s %8s %12s %8s %12s %10s" %
      ("dose", "n_dip", "dist_dip", "n_gao", "dist_gao", "shipped"))
print("-" * 62)
dd, dg = [], []
for i, dv in enumerate(DOSE):
    row = []
    for arm in ("dip", "gao"):
        s = calls[arm][(calls[arm].dose == dv) & (calls[arm].len_mb <= ymax)]
        row.append((len(s), float((s.rel - track).abs().median())))
    dd.append(row[0][1])
    dg.append(row[1][1])
    print("%-6.1f %8d %12.4f %8d %12.4f %10.3f"
          % (dv, row[0][0], row[0][1], row[1][0], row[1][1], SHIPPED_DIST[i]))

print()
worst = max(abs(round(x, 3) - s) for x, s in zip(dd, SHIPPED_DIST))
assert worst < 1e-9, ("the dip arm does not reproduce the published 0.571 / 0.461 / 0.297; "
                      "got %s" % [round(x, 3) for x in dd])
print("[gate] the dip arm reproduces the published distances %s exactly"
      % [round(x, 3) for x in dd])
print("       under the Gao CDR unit they are %s" % [round(x, 3) for x in dg])
print("       decrease from 0 to 2 copies: dip %.3f -> %.3f (%+.3f)   gao %.3f -> %.3f (%+.3f)"
      % (dd[0], dd[-1], dd[-1] - dd[0], dg[0], dg[-1], dg[-1] - dg[0]))
print("       monotone decreasing across dosage: dip %s   gao %s"
      % ("YES" if all(dd[i] > dd[i + 1] for i in range(len(dd) - 1)) else "NO",
         "YES" if all(dg[i] > dg[i + 1] for i in range(len(dg) - 1)) else "NO"))

print()
print("Kruskal-Wallis on distance-to-locus across dosage classes:")
for arm in ("dip", "gao"):
    s = calls[arm][calls[arm].len_mb <= ymax].copy()
    s["dist"] = (s.rel - track).abs()
    pc = kruskal(*[s[s.dose == d].dist for d in DOSE])
    h = s.groupby(["sample", "hap_local"]).agg(dist=("dist", "median"),
                                               dose=("dose", "first")).reset_index()
    ph = kruskal(*[h[h.dose == d].dist for d in DOSE])
    print("   %-4s per call      H = %8.2f  P = %.3g   (n = %d calls)"
          % (arm, pc[0], pc[1], len(s)))
    print("   %-4s per haplotype H = %8.2f  P = %.3g   (n = %d haplotypes)"
          % ("", ph[0], ph[1], len(h)))
print()
print("The per-call test counts a haplotype once per call it carries, so it is pseudoreplicated")
print("in the dip arm (median 4 calls per array) and nearly not at all in the Gao arm (median 1).")
print("The per-haplotype row is the one to quote.")
