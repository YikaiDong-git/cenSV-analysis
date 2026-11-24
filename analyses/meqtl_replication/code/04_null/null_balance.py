#!/usr/bin/env python3
"""B3: does the matched null actually BALANCE, and what binds where it fails?"""
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
OUT = UNIT + "/output/null_design"

np_ = pd.read_csv(OUT + "/b2_matched_null/b2_null_pairs.tsv", sep="\t")
dg = pd.read_csv(OUT + "/b2_matched_null/b2_positivity_diagnostics.tsv", sep="\t")
st = pd.read_csv(UNIT + "/output/targets/prereg_strata.tsv", sep="\t")
pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")

dg = dg.merge(st[["pair_id", "stratum", "priority84", "recoverable_pilot8"]], on="pair_id",
              how="left").merge(pairs[["pair_id", "dist_class", "source_type"]],
                                on="pair_id", how="left")
OKCAL = ["tight", "medium", "loose", "nocensat"]
dg["ok"] = dg.caliper.isin(OKCAL)

print("=== caliper reached, by distance class ===")
print(pd.crosstab(dg.dist_class, dg.caliper).to_string())

print("\n=== failures, by reason and distance class ===")
print(pd.crosstab(dg[~dg.ok].dist_class, dg[~dg.ok].reason).to_string())

print("\n=== pool size at the loosest caliper, for pairs that FAILED there ===")
f = dg[dg.reason == "pool_below_K_at_loosest_caliper"]
if len(f):
    print("  n = " + str(len(f)) + "   pool: min " + str(int(f.pool.min()))
          + ", median " + str(int(f.pool.median())) + ", max " + str(int(f.pool.max())))
    print("  pairs with pool == 0 (hard zero, no candidate at all): "
          + str(int((f.pool == 0).sum())))
    print("  by stratum: " + str(dict(f.stratum.value_counts())))
    print("  by source type: " + str(dict(f.source_type.value_counts())))

print("\n=== COVARIATE BALANCE, discovery vs matched null ===")
u = np_.drop_duplicates("matched_to_pair_id")
rows = []
specs = [("log10_dist", lambda d: np.log10(np.maximum(d.disc_dist_to_sv, 1)),
          lambda d: np.log10(np.maximum(d.null_dist_to_sv, 1))),
         ("mean_meth", lambda d: d.disc_mean_meth, lambda d: d.null_mean_meth),
         ("n_valid", lambda d: d.disc_n_valid, lambda d: d.null_n_valid)]
for name, fa, fb in specs:
    a = fa(u).astype(float).values
    b = fb(np_).astype(float).values
    sd = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2.0)
    smd = (np.mean(a) - np.mean(b)) / sd if sd > 0 else np.nan
    rows.append(dict(covariate=name, disc_mean=round(float(np.mean(a)), 4),
                     null_mean=round(float(np.mean(b)), 4),
                     smd=round(float(smd), 4),
                     verdict=("balanced" if abs(smd) < 0.10 else
                              "marginal" if abs(smd) < 0.25 else "UNBALANCED")))
same_cls = float((np_.null_cls == np_.disc_cls).mean())
rows.append(dict(covariate="censat_class_match", disc_mean=1.0,
                 null_mean=round(same_cls, 4), smd=np.nan,
                 verdict=("balanced" if same_cls > 0.95 else "check")))
bal = pd.DataFrame(rows)
print(bal.to_string(index=False))

print("\n=== balance restricted to the stratum that NEEDS b_hat (secondary) ===")
sec = np_.merge(st[["pair_id", "stratum"]], left_on="matched_to_pair_id",
                right_on="pair_id", how="left")
sec = sec[sec.stratum == "secondary_465_recoverable"]
print("  null pairs in that stratum: " + str(len(sec))
      + "   discovery pairs: " + str(sec.matched_to_pair_id.nunique()))
if len(sec):
    us = sec.drop_duplicates("matched_to_pair_id")
    for name, fa, fb in specs:
        a = fa(us).astype(float).values
        b = fb(sec).astype(float).values
        sd = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2.0)
        smd = (np.mean(a) - np.mean(b)) / sd if sd > 0 else np.nan
        print("  " + name.ljust(12) + " disc " + format(float(np.mean(a)), "9.4f")
              + "  null " + format(float(np.mean(b)), "9.4f")
              + "  SMD " + format(float(smd), "7.4f")
              + ("  balanced" if abs(smd) < 0.10 else
                 "  marginal" if abs(smd) < 0.25 else "  UNBALANCED"))
    print("  censat class match: " + str(round(float((sec.null_cls == sec.disc_cls).mean()), 4)))

print("\n=== nulls per discovery pair ===")
c = np_.groupby("matched_to_pair_id").size()
print("  pairs matched: " + str(len(c)) + "   nulls each: min " + str(int(c.min()))
      + ", median " + str(int(c.median())) + ", max " + str(int(c.max())))
print("  distinct null CpGs: " + str(np_.null_feature_id.nunique())
      + "   reused across pairs: "
      + str(int((np_.null_feature_id.value_counts() > 1).sum())))

tmp = OUT + "/b3_balance/b3_balance.tsv.tmp"
bal.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
os.replace(tmp, OUT + "/b3_balance/b3_balance.tsv")
tmp2 = OUT + "/b3_balance/b3_positivity_by_stratum.tsv.tmp"
z = dg.groupby("stratum").ok.agg(["size", "sum", "mean"]).reset_index()
z.columns = ["stratum", "n_pairs", "n_with_pool", "frac_with_pool"]
z.to_csv(tmp2, sep="\t", index=False, float_format="%.10g")
os.replace(tmp2, OUT + "/b3_balance/b3_positivity_by_stratum.tsv")
print("\nwrote b3_balance.tsv, b3_positivity_by_stratum.tsv")
