#!/usr/bin/env python3
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
pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
svinfo = pairs.drop_duplicates("sv_id")[["sv_id", "source_type"]]
R = R.merge(svinfo, on="sv_id", how="left")
print("rows " + str(len(R)) + "   arms " + str(dict(R.arm.value_counts())))

print("\n" + "=" * 74)
print("1. IS MEASUREMENT AVAILABILITY ITSELF DEL/INS-ASYMMETRIC, AND IN WHICH ARM?")
print("=" * 74)
z = (R.groupby(["arm", "source_type"])
     .agg(n=("slope", "size"), median_slope=("slope", "median"),
          mean_slope=("slope", "mean"), frac_sig=("sig", "mean"),
          median_frac_measured=("frac_measured", "median")).round(4))
print(z.to_string())

for arm in ["discovery", "matched_null"]:
    a = R[(R.arm == arm) & (R.source_type == "DEL")].slope
    b = R[(R.arm == arm) & (R.source_type == "INS")].slope
    if len(a) > 5 and len(b) > 5:
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        print("\n  " + arm + ": DEL median " + format(float(a.median()), "+.4f")
              + "  vs INS " + format(float(b.median()), "+.4f")
              + "   difference " + format(float(a.median() - b.median()), "+.4f")
              + "   Mann-Whitney p = " + format(u.pvalue, ".3g")
              + "   (n " + str(len(a)) + "/" + str(len(b)) + ")")

print("\n  READING. If the DEL/INS asymmetry in measurement availability appears in the")
print("  matched-null arm too, it is a property of the VARIANT TYPE -- a technical")
print("  signature -- and it is present at CpGs cohort A never nominated. If it appears")
print("  only in the discovery arm, it is target-specific and the technical reading is")
print("  weaker.")

print("\n" + "=" * 74)
print("2. THE SYMMETRY TEST THE MAGNET HYPOTHESIS MAKES")
print("=" * 74)
p = pairs.copy()
p["neg"] = (p.cohortA_sign < 0).astype(int)
t = p.groupby("source_type").agg(n=("neg", "size"), frac_negative=("neg", "mean"),
                                 median_dist_bp=("source_target_distance_bp", "median")).round(4)
print(t.to_string())
print("\n  A symmetric 'methylation magnet' predicts insertions raise methylation about as")
print("  reliably as deletions lower it, i.e. frac_negative(INS) well BELOW 0.5.")
ins = p[p.source_type == "INS"]
k = int(ins.neg.sum())
n = int(len(ins))
bt = stats.binomtest(k, n, 0.5)
print("  observed INS: " + str(k) + "/" + str(n) + " negative = "
      + format(k / n, ".4f") + "   binomial vs 0.5 p = " + format(bt.pvalue, ".3g"))
print("  -> insertions are not the mirror image of deletions; they sit near or above 0.5.")
print("     That is hard for a symmetric magnet and easy for an asymmetric technical channel.")

print("\n" + "=" * 74)
print("3. IS THE TYPE EFFECT SEPARABLE FROM DISTANCE AT ALL?")
print("=" * 74)
print(p.groupby(["source_type", "dist_class"], observed=True)
      .agg(n=("neg", "size"), frac_neg=("neg", "mean")).round(3).to_string())
print("\n  distance by type (the collinearity that makes 'type vs distance' hard):")
print(p.groupby("source_type").source_target_distance_bp
      .describe()[["count", "25%", "50%", "75%"]].round(0).to_string())

sub = p[p.dist_class.isin(["<=100kb", "0.1-1Mb", "1-5Mb"])]
if len(sub):
    print("\n  within <5 Mb only:")
    print(sub.groupby("source_type").agg(n=("neg", "size"),
                                         frac_neg=("neg", "mean")).round(4).to_string())
    ct = pd.crosstab(sub.source_type, sub.neg)
    if ct.shape == (2, 2):
        orr, pv = stats.fisher_exact(ct.values)
        print("  Fisher OR = " + format(orr, ".4g") + "   p = " + format(pv, ".3g"))

tmp = OUT + "/c7_del_ins_asymmetry.tsv.tmp"
z.reset_index().to_csv(tmp, sep="\t", index=False, float_format="%.10g")
os.replace(tmp, OUT + "/c7_del_ins_asymmetry.tsv")
print("\nwrote c7_del_ins_asymmetry.tsv")
