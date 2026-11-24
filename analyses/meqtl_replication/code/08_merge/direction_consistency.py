#!/usr/bin/env python3
"""How strong is the direction evidence, once the honest unit is used?"""
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")

p = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
p["neg"] = (p.expected_rint_beta < 0).astype(int)
loc = p[p.source_target_distance_bp < 5e6].copy()

for label, dat in [("ALL pairs", p), ("within 5 Mb", loc)]:
    print("=" * 72)
    print(label)
    print("=" * 72)
    for t in ["DEL", "INS"]:
        s = dat[dat.source_type == t]
        if not len(s):
            continue
        per = s.groupby("sv_id").neg.agg(["size", "mean"]).rename(
            columns={"size": "n_pairs", "mean": "frac_neg"})
        n_sv = len(per)
        top = per.n_pairs.sort_values(ascending=False)
        share_top3 = float(top.head(3).sum()) / float(top.sum())
        maj = int((per.frac_neg >= 0.8).sum())
        min_ = int((per.frac_neg <= 0.2).sum())
        mixed = n_sv - maj - min_
        sv_neg = (per.frac_neg > 0.5).sum()
        bt = stats.binomtest(int(sv_neg), n_sv, 0.5) if n_sv else None
        print("\n  " + t + ":  " + str(len(s)) + " pairs from " + str(n_sv) + " SVs"
              + "   pair-level frac negative = " + format(float(s.neg.mean()), ".4f"))
        print("    top-3 SVs hold " + format(share_top3 * 100, ".1f") + "% of the pairs"
              + "   (max single SV: " + str(int(top.iloc[0])) + " pairs)")
        print("    SVs >=80% negative: " + str(maj) + " / " + str(n_sv)
              + "   <=20% negative: " + str(min_) + "   mixed: " + str(mixed))
        if bt is not None:
            print("    SV-level sign test (each SV = 1 unit): " + str(int(sv_neg)) + "/"
                  + str(n_sv) + " majority-negative, binomial p = " + format(bt.pvalue, ".3g"))
        print("    per-SV frac_neg: min " + format(float(per.frac_neg.min()), ".2f")
              + "  median " + format(float(per.frac_neg.median()), ".2f")
              + "  max " + format(float(per.frac_neg.max()), ".2f"))
    print()

print("=" * 72)
print("THE DISTAL 'FLIP': is it a flip, or is it attenuation into noise?")
print("=" * 72)
z = p.groupby("dist_class", observed=True).agg(
    n=("neg", "size"), frac_neg=("neg", "mean"),
    median_abs_beta=("expected_rint_beta", lambda s: float(s.abs().median()))).round(4)
print(z.to_string())
print("\n  A true sign FLIP means the distal fraction goes clearly BELOW 0.5 with")
print("  comparable effect magnitude. Attenuation into noise means it approaches 0.5")
print("  while |beta| shrinks. Which pattern is present decides whether 'direction")
print("  changes at the satellite boundary' is a finding or a restatement of decay.")
far = p[p.source_target_distance_bp > 5e6]
if len(far):
    k = int(far.neg.sum())
    n = int(len(far))
    bt = stats.binomtest(k, n, 0.5)
    print("\n  >5 Mb: " + str(k) + "/" + str(n) + " negative = " + format(k / n, ".4f")
          + "   binomial vs 0.5 p = " + format(bt.pvalue, ".3g"))
