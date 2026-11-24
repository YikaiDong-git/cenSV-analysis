#!/usr/bin/env python3
"""GATE: is cohort A's finding an mQTL, or coordinate bookkeeping?"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
OUT = UNIT + "/output/arms/gate_inside_deletion"

FLANK = 500


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


os.makedirs(OUT, exist_ok=True)
pr = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
tg = pd.read_csv(UNIT + "/output/targets/targets_cpg.tsv", sep="\t")
pr = pr.merge(tg[["target_id", "chm13_pos"]].rename(
    columns={"chm13_pos": "chm13_pos_target"}), on="target_id", how="left")
assert pr.chm13_pos_target.notna().all(), "a pair lost its target coordinate"
print("pairs: " + str(len(pr)))

same_chrom = pr.source_chrom == pr.chm13_chrom_target
pr["sv_span"] = pr.source_end - pr.source_pos
pr["inside_deletion"] = (same_chrom & (pr.source_type == "DEL")
                         & (pr.chm13_pos_target >= pr.source_pos)
                         & (pr.chm13_pos_target < pr.source_end))
pr["inside_any_sv"] = (same_chrom
                       & (pr.chm13_pos_target >= pr.source_pos)
                       & (pr.chm13_pos_target < pr.source_end))
pr["within_flank"] = (same_chrom
                      & (pr.chm13_pos_target >= pr.source_pos - FLANK)
                      & (pr.chm13_pos_target < pr.source_end + FLANK))

print("\n=== GATE 1: do target CpGs sit inside the source SV? ===")
print("SV span (bp):")
print(pr.groupby("source_type").sv_span.describe()[["count", "min", "50%", "max"]].to_string())
print("\ntarget CpG inside the deleted interval : "
      + str(int(pr.inside_deletion.sum())) + " / " + str(len(pr)))
print("target CpG inside ANY source SV span   : "
      + str(int(pr.inside_any_sv.sum())) + " / " + str(len(pr)))
print("target CpG within +/-" + str(FLANK) + " bp of the SV  : "
      + str(int(pr.within_flank.sum())) + " / " + str(len(pr)))

print("\nby SV type and cohort-A sign:")
t = (pr.groupby(["source_type", "cohortA_sign"])
     .agg(n=("pair_id", "size"), inside_del=("inside_deletion", "sum"),
          within_flank=("within_flank", "sum"),
          median_abs_dist=("source_target_distance_bp", "median")).reset_index())
print(t.to_string(index=False))
write_atomic(t, OUT + "/g1_inside_by_type_sign.tsv")

print("\nDEL vs INS negative fraction (the signature the panel predicted):")
s = (pr.groupby("source_type")
     .agg(n=("pair_id", "size"),
          frac_negative=("cohortA_sign", lambda x: float((x < 0).mean())),
          median_dist=("source_target_distance_bp", "median")).reset_index())
print(s.to_string(index=False))
write_atomic(s, OUT + "/g1_type_sign_summary.tsv")

print("\nminimum |distance| from source SV, by type (bp):")
q = pr.groupby("source_type").source_target_distance_bp.quantile(
    [0, 0.01, 0.05, 0.25, 0.5]).unstack()
print(q.to_string())
write_atomic(q.reset_index(), OUT + "/g1_distance_quantiles.tsv")

th = pd.read_csv(UNIT + "/output/a1_lift/pilot_qc/target_hit_summary.tsv", sep="\t")
rec = set(th[th.n_hit > 0].target_id)
pr["recoverable"] = pr.target_id.isin(rec)

ct = pd.crosstab(pr.cohortA_sign, pr.recoverable)
print("\n=== GATE 2: does the pericentromeric window select on cohort-A sign? ===")
print(ct.to_string())
a = int(ct.loc[-1, True]); b = int(ct.loc[-1, False])
c = int(ct.loc[1, True]); d = int(ct.loc[1, False])
orr, p = stats.fisher_exact([[a, b], [c, d]])
print("\nA-negative retained: " + str(a) + "/" + str(a + b)
      + " = " + str(round(100.0 * a / (a + b), 1)) + "%")
print("A-positive retained: " + str(c) + "/" + str(c + d)
      + " = " + str(round(100.0 * c / (c + d), 1)) + "%")
print("differential       : " + str(round((a / (a + b)) / (c / (c + d)), 2)) + "x")
print("Fisher exact       : OR = " + str(round(orr, 2)) + ",  p = " + format(p, ".3g"))
print("\nf_A in the FULL set      : " + str(round(float((pr.cohortA_sign < 0).mean()), 5)))
print("f_A in the RECOVERABLE set: "
      + str(round(float((pr[pr.recoverable].cohortA_sign < 0).mean()), 5)))

print("\nwhy: cohort-A sign by distance class, and recovery by distance class")
z = (pr.groupby("dist_class", observed=True)
     .agg(n=("pair_id", "size"),
          frac_A_negative=("cohortA_sign", lambda x: float((x < 0).mean())),
          frac_recoverable=("recoverable", "mean")).reset_index())
print(z.to_string(index=False))
write_atomic(z, OUT + "/g2_sign_and_recovery_by_distance.tsv")

write_atomic(pd.DataFrame([dict(
    n_pairs=len(pr),
    inside_deletion=int(pr.inside_deletion.sum()),
    inside_any_sv=int(pr.inside_any_sv.sum()),
    within_flank=int(pr.within_flank.sum()),
    fA_full=float((pr.cohortA_sign < 0).mean()),
    fA_recoverable=float((pr[pr.recoverable].cohortA_sign < 0).mean()),
    Aneg_retained=a, Aneg_total=a + b, Apos_retained=c, Apos_total=c + d,
    fisher_or=float(orr), fisher_p=float(p),
)]), OUT + "/gate_summary.tsv")
