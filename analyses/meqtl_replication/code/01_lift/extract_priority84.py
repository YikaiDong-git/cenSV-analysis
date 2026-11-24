#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
E = M + "/results/meqtl_mac8_elastic_tag_conditioning/summary"
UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
OUT = UNIT + "/output/targets"


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


cc = pd.read_csv(E + "/conditioning_consensus.tsv", sep="\t")
print("elastic-tag consensus: " + str(cc.shape))
cls_cols = [c for c in cc.columns if "class" in c.lower()]
print("classification columns: " + str(cls_cols))
for c in cls_cols:
    print("  " + c + ": " + str(dict(cc[c].value_counts(dropna=False))))

target_col = None
for c in cls_cols:
    if (cc[c].astype(str) == "primary_model_supported").sum() == 84:
        target_col = c
        break
assert target_col is not None, "no classification column yields 84 primary_model_supported"
print("\nusing column: " + target_col)

sup = cc[cc[target_col].astype(str) == "primary_model_supported"].copy()
print("supported pairs: " + str(len(sup))
      + "   by scope: " + str(dict(sup.scope.value_counts())) if "scope" in sup
      else str(len(sup)))

pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
th = pd.read_csv(UNIT + "/output/a1_lift/pilot_qc/target_hit_summary.tsv", sep="\t")
rec = set(th[th.n_hit > 0].target_id)

assert "pair_id" in sup.columns
p84 = set(sup.pair_id)
print("\n84-set INTERSECT the 465 pre-registered pair_ids: " + str(len(p84 & set(pairs.pair_id))))

pairs["priority84"] = pairs.pair_id.isin(p84)
pairs["recoverable_pilot8"] = pairs.target_id.isin(rec)

print("\n--- composition of the 84 vs the full 465 ---")
tab = []
for name, sub in [("all 465", pairs), ("priority 84", pairs[pairs.priority84])]:
    tab.append(dict(
        stratum=name, n=len(sub),
        n_cis=int((sub.scope == "cis").sum()) if "scope" in sub else -1,
        n_trans=int((sub.scope == "trans").sum()) if "scope" in sub else -1,
        frac_A_negative=round(float((sub.cohortA_sign < 0).mean()), 4),
        n_recoverable=int(sub.recoverable_pilot8.sum()),
        frac_recoverable=round(float(sub.recoverable_pilot8.mean()), 4),
        n_source_SV=int(sub.sv_id.nunique())))
print(pd.DataFrame(tab).to_string(index=False))

print("\n--- the 84 by distance class ---")
z = (pairs[pairs.priority84].groupby("dist_class", observed=True)
     .agg(n=("pair_id", "size"),
          frac_A_neg=("cohortA_sign", lambda x: float((x < 0).mean())),
          n_recoverable=("recoverable_pilot8", "sum")).reset_index())
print(z.to_string(index=False))

rec84 = pairs[pairs.priority84 & pairs.recoverable_pilot8]
rec465 = pairs[pairs.recoverable_pilot8]
print("\n=== the decision numbers ===")
print("  primary candidate  : 84-set AND recoverable = " + str(len(rec84))
      + " pairs from " + str(rec84.sv_id.nunique()) + " source SVs")
if len(rec84):
    print("    f_A on that set  = " + str(round(float((rec84.cohortA_sign < 0).mean()), 5)))
print("  secondary          : 465 AND recoverable  = " + str(len(rec465))
      + " pairs from " + str(rec465.sv_id.nunique()) + " source SVs")
print("    f_A on that set  = " + str(round(float((rec465.cohortA_sign < 0).mean()), 5)))

pairs["stratum"] = np.where(pairs.priority84 & pairs.recoverable_pilot8, "primary_84_recoverable",
                    np.where(pairs.recoverable_pilot8, "secondary_465_recoverable",
                             "not_testable_pilot8"))
print("\n  frozen strata: " + str(dict(pairs.stratum.value_counts())))
write_atomic(pairs[["pair_id", "sv_id", "target_id", "scope", "dist_class",
                    "cohortA_sign", "priority84", "recoverable_pilot8", "stratum"]],
             OUT + "/prereg_strata.tsv")
print("\nNOTE recoverable_pilot8 is PROVISIONAL (8 haplotypes). v2 recomputes it on 311;")
print("     the stratum LABELS are frozen, the recoverability filter is re-derived.")
