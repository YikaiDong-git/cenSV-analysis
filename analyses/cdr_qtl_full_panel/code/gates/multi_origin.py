#!/usr/bin/env python
"""Which SVs turned out to have MORE THAN ONE origin on the centromere haplotype tree?"""
import os
import re
import sys
from collections import Counter

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200")
DEC = ROOT + "/sv_extraction/decompose/" + SUF
WU = ROOT + "/analyses/cdr_qtl_full_panel"
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

frames = []
for c in CHROMS:
    p = "%s/%s/%s_decomposition.tsv" % (DEC, c, c)
    if not os.path.exists(p):
        sys.exit("FAIL missing decomposition log for %s at %s -- the count would silently be "
                 "over 22 chromosomes while claiming 23" % (c, p))
    d = pd.read_csv(p, sep="\t")
    d["chrom"] = c
    frames.append(d)
dec = pd.concat(frames, ignore_index=True)
print("decomposition log: %d rows over %d chromosomes" % (len(dec), dec.chrom.nunique()))
print("row types: %s" % dict(Counter(dec.type)))
print("iterations present: %s" % sorted(dec.iteration.unique()))

mono = dec[dec.type == "monophyletic"]
per = mono.groupby(["chrom", "sv_id"]).size().rename("n_origins").reset_index()
print("\nSVs with at least one monophyletic component: %d" % len(per))
print("distribution of components per SV:")
for k, v in sorted(Counter(per.n_origins).items()):
    print("   %d origin%s : %6d SVs%s" % (k, "s" if k > 1 else " ", v,
                                          "   <-- multiple origins" if k > 1 else ""))

multi = per[per.n_origins >= 2].copy()
print("\n=== %d SVs have MULTIPLE origins (%.2f%% of the %d with any) ==="
      % (len(multi), 100.0 * len(multi) / len(per), len(per)))
print("\nby chromosome:")
bc = per.assign(multi=per.n_origins >= 2).groupby("chrom").agg(
    svs=("sv_id", "size"), multi=("multi", "sum"))
bc["pct"] = (100.0 * bc.multi / bc.svs).round(1)
bc = bc.reindex([c for c in CHROMS if c in bc.index])
print(bc.to_string())


def svtype(sv_id):
    m = re.match(r"svim_asm\.([A-Z_]+)\.", str(sv_id).split(";")[0])
    return m.group(1) if m else "?"


multi["svtype"] = multi.sv_id.map(svtype)
print("\nby SV type:")
allt = per.assign(svtype=per.sv_id.map(svtype)).groupby("svtype").size()
mt = multi.groupby("svtype").size()
for t in allt.sort_values(ascending=False).index:
    print("  %-12s %6d of %6d  (%.2f%%)"
          % (t, mt.get(t, 0), allt[t], 100.0 * mt.get(t, 0) / allt[t]))

key = set(zip(multi.chrom, multi.sv_id))
comp = mono[[(c, s) in key for c, s in zip(mono.chrom, mono.sv_id)]].copy()
comp = comp.sort_values(["chrom", "sv_id", "iteration"])
print("\ncomponent sizes across the %d multi-origin SVs:" % len(multi))
first = comp[comp.iteration == 0]
later = comp[comp.iteration > 0]
print("  first component  carriers median %.0f (IQR %.0f-%.0f), R2 median %.2f"
      % (first.n_carriers.median(), first.n_carriers.quantile(.25),
         first.n_carriers.quantile(.75), first.r2.median()))
print("  later components carriers median %.0f (IQR %.0f-%.0f), R2 median %.2f"
      % (later.n_carriers.median(), later.n_carriers.quantile(.25),
         later.n_carriers.quantile(.75), later.r2.median()))
print("  a later component is larger than the first in %d of %d SVs"
      % (sum(1 for _, g in comp.groupby(["chrom", "sv_id"])
             if g.n_carriers.iloc[1:].max() > g.n_carriers.iloc[0]), len(multi)))

bad = []
for (c, s), g in comp.groupby(["chrom", "sv_id"]):
    sets = [set(map(int, str(x).split(","))) for x in g.carrier_indices]
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            if sets[i] & sets[j]:
                bad.append((c, s, len(sets[i] & sets[j])))
print("\ncomponents of the same SV sharing a carrier: %d (must be 0; the mask is what "
      "separates them)" % len(bad))
assert not bad, bad[:5]

OUT = WU + "/output/multi_origin_svs.tsv"
comp.to_csv(OUT, sep="\t", index=False)
print("\nwrote %s" % OUT)

vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
vt["orig_id"] = vt.sv_id.astype(str).str.replace(r"_(B\d+|residual|original)$", "", regex=True)
sub_multi = vt[[(c, s) in key for c, s in zip(vt.chrom, vt.orig_id)]]
print("\n=== in the %d-locus QTL substrate ===" % vt.sv_locus.nunique())
print("loci derived from a multi-origin SV: %d, from %d distinct SVs on %d chromosomes"
      % (len(sub_multi), sub_multi.orig_id.nunique(), sub_multi.chrom.nunique()))
if len(sub_multi):
    print(sub_multi.groupby("chrom").agg(loci=("sv_locus", "size"),
                                         svs=("orig_id", "nunique")).to_string())
    print("\nlargest by carrier count:")
    top = (multi.merge(sub_multi.groupby(["chrom", "orig_id"]).svlen.first().reset_index(),
                       left_on=["chrom", "sv_id"], right_on=["chrom", "orig_id"])
           .merge(comp.groupby(["chrom", "sv_id"]).n_carriers.sum().rename("carriers")
                  .reset_index(), on=["chrom", "sv_id"]))
    for r in top.sort_values("carriers", ascending=False).head(15).itertuples():
        print("  %-6s %-8s %7d bp  %d origins, %3d carriers total"
              % (r.chrom, r.svtype, int(r.svlen), r.n_origins, r.carriers))
