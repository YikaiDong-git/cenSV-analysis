#!/usr/bin/env python3
"""Why GATE 3 returned almost nothing for the headline SVs: is it \"no channel\", or \"no data\"?"""
import glob
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
H = UNIT + "/output/headline_svs"

D = pd.concat([pd.read_csv(f, sep="\t",
                           usecols=["meth_tag", "cpg_key", "target_class", "hit_anchor",
                                    "n_valid", "lift_status", "plus_c_pos", "chrom"])
               for f in sorted(glob.glob(H + "/per_hap/*.tsv.gz"))], ignore_index=True)
T = pd.read_csv(H + "/null/headline_targets_lift.tsv", sep="\t")
D["sv_name"] = D.cpg_key.map(dict(zip(T.cpg_key, T.sv_name)))
D["measured"] = (D.hit_anchor == 1) & (D.n_valid > 0)
n_hap = D.meth_tag.nunique()
print("haplotypes " + str(n_hap) + "   CpGs " + str(D.cpg_key.nunique()))

print("\n=== lift status, by SV and arm ===")
print(pd.crosstab([D.sv_name, D.target_class], D.lift_status,
                  normalize="index").round(4).to_string())

per = (D.groupby(["sv_name", "target_class", "cpg_key"], as_index=False)
       .measured.mean().rename(columns={"measured": "frac_hap_measured"}))
print("\n=== per-CpG fraction of the 311 haplotypes with a methylation value ===")
print(per.groupby(["sv_name", "target_class"]).frac_hap_measured
      .describe()[["count", "mean", "25%", "50%", "75%", "max"]].round(4).to_string())

print("\n=== how many CpGs are measurable AT ALL? ===")
for (sv, arm), s in per.groupby(["sv_name", "target_class"]):
    n = len(s)
    print("  " + sv.ljust(13) + " " + arm.ljust(13)
          + "  never measured: " + str(int((s.frac_hap_measured == 0).sum())).rjust(5)
          + " / " + str(n).ljust(6)
          + "  >=10% of haps: " + str(int((s.frac_hap_measured >= 0.10).sum())).rjust(4)
          + "  >=50%: " + str(int((s.frac_hap_measured >= 0.50).sum())).rjust(4))

print("\n=== distance from the source SV, measured vs not ===")
SVPOS = {"chr12_DEL333": 37184003, "chr8_DEL113": 44300277}
p2 = per.merge(T[["cpg_key", "chrom", "plus_c_pos"]], on="cpg_key", how="left")
p2["dist_mb"] = [abs(int(r.plus_c_pos) - SVPOS[r.sv_name]) / 1e6 for r in p2.itertuples()]
p2["bin"] = pd.cut(p2.dist_mb, [-.01, 0.1, 1, 5, 25, 200],
                   labels=["<=100kb", "0.1-1Mb", "1-5Mb", "5-25Mb", ">25Mb"])
print(p2.groupby(["sv_name", "bin"], observed=True)
      .agg(n_cpg=("cpg_key", "size"),
           frac_ever_measured=("frac_hap_measured", lambda s: float((s > 0).mean())),
           median_frac=("frac_hap_measured", "median")).round(4).to_string())

print("\n=== VERDICT ===")
tot = per[per.target_class == "discovery"]
ever = int((tot.frac_hap_measured > 0).sum())
print("  discovery CpGs measurable in >=1 haplotype: " + str(ever) + " / " + str(len(tot))
      + "  (" + format(100.0 * ever / len(tot), ".1f") + "%)")
usable = int((tot.frac_hap_measured >= 0.10).sum())
print("  discovery CpGs measured in >=10% of haplotypes: " + str(usable)
      + " / " + str(len(tot)))
print("\n  If these are near zero, GATE 3 on these two SVs is UNTESTABLE, not clean.")
print("  The HPRC substrate is a pericentromeric window; these SVs' target CpGs are")
print("  mostly outside it. That is a statement about coverage, not about the channel.")
tmp = H + "/results/headline_testability.tsv.tmp"
p2.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
os.replace(tmp, H + "/results/headline_testability.tsv")
