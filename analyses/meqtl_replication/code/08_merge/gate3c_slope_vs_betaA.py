#!/usr/bin/env python3
"""GATE 3c: does cohort B's genotype-dependent MISSINGNESS predict cohort A's EFFECT SIZE?"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
OUT = UNIT + "/output/merge"

R = pd.read_csv(OUT + "/c4_gate3b_missingness_by_arm.tsv", sep="\t")
d = R[R.arm == "discovery"].copy()

pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
tgt = pd.read_csv(UNIT + "/output/targets/targets_v2_lift.tsv", sep="\t")
st = pd.read_csv(UNIT + "/output/targets/prereg_strata.tsv", sep="\t")
f2k = dict(zip(tgt.feature_id, tgt.cpg_key))
pairs["cpg_key"] = pairs.feature_id.map(f2k)
P = pairs[pairs.cpg_key.notna()].merge(st[["pair_id", "stratum", "priority84"]],
                                       on="pair_id", how="left")

m = d.merge(P[["sv_id", "cpg_key", "pair_id", "expected_rint_beta", "expected_rint_p",
               "source_target_distance_bp", "source_type", "dist_class", "cohortA_sign",
               "stratum", "priority84"]],
            on=["sv_id", "cpg_key"], how="inner")
print("discovery pairs with both a missingness slope and a cohort-A beta: " + str(len(m)))

m["absbetaA"] = m.expected_rint_beta.abs()
m["neglog10pA"] = -np.log10(m.expected_rint_p.clip(lower=1e-300))

print("\n=== does cohort B missingness slope track cohort A effect magnitude? ===")
for yname, y in [("|beta_A|", m.absbetaA), ("-log10 p_A", m.neglog10pA)]:
    rho, p = stats.spearmanr(m.slope, y)
    print("  Spearman(slope_B , " + yname.ljust(11) + ") = " + format(rho, "+.4f")
          + "   p = " + format(p, ".3g"))
    rho2, p2 = stats.spearmanr(m.slope.abs(), y)
    print("    using |slope_B|                    = " + format(rho2, "+.4f")
          + "   p = " + format(p2, ".3g"))

print("\n=== cohort A effect, split by whether cohort B shows the channel ===")
m["channel"] = np.where(m.slope < -0.05, "carriers_lose_measurement",
                        np.where(m.slope > 0.05, "carriers_gain", "flat"))
z = m.groupby("channel").agg(n=("slope", "size"),
                             median_absbetaA=("absbetaA", "median"),
                             median_neglog10pA=("neglog10pA", "median"),
                             frac_A_negative=("cohortA_sign", lambda s: float((s < 0).mean())),
                             median_slope=("slope", "median")).round(4)
print(z.to_string())

a = m[m.channel == "carriers_lose_measurement"].absbetaA
b = m[m.channel != "carriers_lose_measurement"].absbetaA
if len(a) > 5 and len(b) > 5:
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    print("\n  |beta_A| in channel-positive vs rest: median "
          + format(float(a.median()), ".4f") + " vs " + format(float(b.median()), ".4f")
          + "   Mann-Whitney p = " + format(u.pvalue, ".3g"))

print("\n=== the channel by SV type and distance ===")
print(m.groupby("source_type").slope.agg(["size", "median"]).round(4).to_string())
print(m.groupby("dist_class", observed=True).slope.agg(["size", "median"]).round(4).to_string())

print("\n=== and by pre-registered stratum ===")
print(m.groupby("stratum").agg(n=("slope", "size"), median_slope=("slope", "median"),
                               frac_channel=("channel", lambda s: float(
                                   (s == "carriers_lose_measurement").mean()))).round(4).to_string())
print("\n  The 84-priority stratum matters most here: it is the artefact-immune endpoint,")
print("  so if the channel is weaker there, the surviving claim is cleaner.")
print(m.groupby("priority84").agg(n=("slope", "size"), median_slope=("slope", "median"),
                                  frac_channel=("channel", lambda s: float(
                                      (s == "carriers_lose_measurement").mean()))).round(4).to_string())

tmp = OUT + "/c5_gate3c_slope_vs_betaA.tsv.tmp"
m.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
os.replace(tmp, OUT + "/c5_gate3c_slope_vs_betaA.tsv")
print("\nwrote c5_gate3c_slope_vs_betaA.tsv")
print("\nCAVEAT: correlational. A positive association is consistent with the collider")
print("reading but does not establish it; a null association weakens the collider reading")
print("without proving the effect is real.")
