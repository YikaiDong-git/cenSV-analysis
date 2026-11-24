#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
OUT = UNIT + "/output/merge"

m = pd.read_csv(OUT + "/c5_gate3c_slope_vs_betaA.tsv", sep="\t")
tgt = pd.read_csv(UNIT + "/output/targets/targets_v2_lift.tsv", sep="\t")
print("pairs: " + str(len(m)) + "   source SVs: " + str(m.sv_id.nunique()))

m["y_neg"] = (m.cohortA_sign < 0).astype(int)
m["log10d"] = np.log10(m.source_target_distance_bp.clip(lower=1))
m["is_del"] = (m.source_type == "DEL").astype(int)
m = m[m.log10d.notna() & m.slope.notna()].copy()
print("modelled: " + str(len(m)) + "   negative-sign rate " + str(round(float(m.y_neg.mean()), 4)))

print("\n=== marginal association (NOT evidence, shown for contrast) ===")
ct = pd.crosstab(m.slope < -0.05, m.y_neg)
print(ct.to_string())
if ct.shape == (2, 2):
    orr, p = stats.fisher_exact(ct.values)
    print("  Fisher OR = " + format(orr, ".4g") + "   p = " + format(p, ".3g"))

print("\n=== logistic: P(cohort-A effect negative), SE clustered on source SV ===")
base = "y_neg ~ is_del + log10d"
full = "y_neg ~ slope + is_del + log10d"
res = {}
for name, f in [("without slope_B", base), ("with slope_B", full)]:
    rob = smf.logit(f, data=m).fit(cov_type="cluster", cov_kwds={"groups": m.sv_id.values}, disp=0)
    res[name] = rob
    print("\n  " + name + "   (n=" + str(int(rob.nobs))
          + ", clusters=" + str(m.sv_id.nunique()) + ")")
    for term, coef, se, pv in zip(rob.model.exog_names, np.asarray(rob.params), np.asarray(rob.bse), np.asarray(rob.pvalues)):
        print("    " + term.ljust(12) + " beta " + format(coef, "+9.4f")
              + "   SE " + format(se, "7.4f") + "   p " + format(pv, ".3g"))

r = res["with slope_B"]
i = list(r.model.exog_names).index("slope")
print("\n  CONDITIONAL EFFECT OF slope_B: beta " + format(float(np.asarray(r.params)[i]), "+.4f")
      + "  p = " + format(float(np.asarray(r.pvalues)[i]), ".3g"))
print("  (negative beta on slope means: the MORE carriers lose measurement -- i.e. the more")
print("   negative slope_B -- the MORE likely cohort A reported a negative effect.)")

print("\n=== within the 84-priority stratum (the artefact-immune endpoint) ===")
p84 = m[m.priority84 == True]
print("  n = " + str(len(p84)) + "   SVs " + str(p84.sv_id.nunique())
      + "   negative-sign rate " + str(round(float(p84.y_neg.mean()), 4)))
if len(p84) > 12:
    print(pd.crosstab(p84.slope < -0.05, p84.y_neg).to_string())
    z = p84.groupby(p84.slope < -0.05).y_neg.agg(["size", "mean"]).round(4)
    print(z.to_string())

print("\n=== is the channel simply distance? ===")
rho, pv = stats.spearmanr(m.log10d, m.slope)
print("  Spearman(log10 distance, slope_B) = " + format(rho, "+.4f")
      + "   p = " + format(pv, ".3g"))
rho2, pv2 = stats.spearmanr(m.frac_measured, m.slope)
print("  Spearman(frac_measured, slope_B)  = " + format(rho2, "+.4f")
      + "   p = " + format(pv2, ".3g"))

tmp = OUT + "/c6_gate3d_conditional.tsv.tmp"
m.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
os.replace(tmp, OUT + "/c6_gate3d_conditional.tsv")
print("\nwrote c6_gate3d_conditional.tsv")
