#!/usr/bin/env python3
"""Dose-response test of the alpha-monomer \"methylation magnet\" hypothesis."""
import os
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
SUB = (CENSV_ROOT + "/analyses/"
       "cdr_qtl_full_panel")
OUT = UNIT + "/output/merge"

p = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
p["svlen"] = (p.source_end - p.source_pos).abs()
var = pd.read_csv(SUB + "/geno/full_variants.tsv", sep="\t", engine="python")
key = var.drop_duplicates(["chrom", "pos", "end", "svtype"])[
    ["chrom", "pos", "end", "svtype", "svlen"]].rename(columns={"svlen": "cat_svlen"})
p = p.merge(key, left_on=["source_chrom", "source_pos", "source_end", "source_type"],
            right_on=["chrom", "pos", "end", "svtype"], how="left")
p["length"] = p.cat_svlen.fillna(p.svlen)
p["absbeta"] = p.expected_rint_beta.abs()
p["neg"] = (p.expected_rint_beta < 0).astype(int)
p["log10d"] = np.log10(p.source_target_distance_bp.clip(lower=1))
print("pairs " + str(len(p)) + "   with a catalog length: " + str(int(p.cat_svlen.notna().sum())))

d = p[(p.source_type == "DEL") & p.length.notna()].copy()
i = p[(p.source_type == "INS") & p.length.notna()].copy()
print("\nSV length (bp), by type, over distinct SVs:")
print(p.drop_duplicates("sv_id").groupby("source_type").length
      .describe()[["count", "min", "25%", "50%", "75%", "max"]].round(1).to_string())

print("\n" + "=" * 70)
print("1. DOSE-RESPONSE: |beta_A| vs deletion length")
print("=" * 70)
for name, sub in [("all deletions", d), ("deletions <5 Mb", d[d.source_target_distance_bp < 5e6])]:
    if len(sub) < 20:
        continue
    rho, pv = stats.spearmanr(sub.length, sub.absbeta)
    print("  " + name.ljust(18) + " n=" + str(len(sub))
          + "  Spearman(length, |beta|) = " + format(rho, "+.4f") + "  p = " + format(pv, ".3g"))
    g = sub.groupby("sv_id").agg(length=("length", "first"),
                                 med_absbeta=("absbeta", "median")).reset_index()
    if len(g) > 5:
        rho2, pv2 = stats.spearmanr(g.length, g.med_absbeta)
        print("    per-SV (n=" + str(len(g)) + " SVs): " + format(rho2, "+.4f")
              + "  p = " + format(pv2, ".3g") + "   <- the honest unit")

print("\n" + "=" * 70)
print("2. DOES LENGTH PREDICT SIGN, BEYOND TYPE AND DISTANCE?")
print("=" * 70)
sub = p[p.length.notna() & p.log10d.notna()].copy()
sub["is_del"] = (sub.source_type == "DEL").astype(int)
sub["log10len"] = np.log10(sub.length.clip(lower=1))
try:
    import statsmodels.formula.api as smf
    m = smf.logit("neg ~ log10len + is_del + log10d", data=sub).fit(
        cov_type="cluster", cov_kwds={"groups": sub.sv_id.values}, disp=0)
    for t, c, s, pv in zip(m.model.exog_names, np.asarray(m.params),
                           np.asarray(m.bse), np.asarray(m.pvalues)):
        print("  " + t.ljust(11) + " beta " + format(c, "+9.4f")
              + "  SE " + format(s, "7.4f") + "  p " + format(pv, ".3g"))
    print("  (clusters = " + str(sub.sv_id.nunique()) + " source SVs)")
except Exception as e:
    print("  model failed: " + str(e))

print("\n" + "=" * 70)
print("3. MONOMER REGISTER: is length near a multiple of 171 bp?")
print("=" * 70)
u = p.drop_duplicates("sv_id")[["sv_id", "source_type", "length"]].dropna()
u["n_monomer"] = u.length / 171.0
u["dist_to_int"] = (u.n_monomer - u.n_monomer.round()).abs() * 171.0
u["register"] = u.dist_to_int <= 10
print(u.groupby("source_type").agg(n=("length", "size"),
                                   median_monomers=("n_monomer", "median"),
                                   frac_in_register=("register", "mean")).round(3).to_string())
print("\n  The MAC8 report tested this directly and found 79.7% of eligible SVs within 10 bp")
print("  of a 171-bp multiple vs 82.8% in type- and length-matched live-HOR controls")
print("  (difference -3.1 pp, empirical p = 0.590) -- i.e. NOT enriched. Register-preserving")
print("  lengths are simply common in this callset.")

print("\n" + "=" * 70)
print("4. THE SYMMETRY THE MAGNET REQUIRES, WITHIN A MATCHED DISTANCE STRATUM")
print("=" * 70)
for lo, hi, lab in [(0, 1e5, "<=100 kb"), (1e5, 1e6, "0.1-1 Mb"), (1e6, 5e6, "1-5 Mb")]:
    s = p[(p.source_target_distance_bp >= lo) & (p.source_target_distance_bp < hi)]
    if not len(s):
        continue
    row = s.groupby("source_type").neg.agg(["size", "mean"])
    line = "  " + lab.ljust(10)
    for t in ["DEL", "INS"]:
        if t in row.index:
            line += "  " + t + " " + str(int(row.loc[t, "size"])).rjust(4) + " pairs, " \
                    + format(float(row.loc[t, "mean"]), ".3f") + " negative"
    print(line)
print("\n  A symmetric magnet predicts INS well below 0.5 in every stratum. It is not:")
print("  insertions are also predominantly negative wherever they are close in.")

tmp = OUT + "/c8_magnet_dose_response.tsv.tmp"
u.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
os.replace(tmp, OUT + "/c8_magnet_dose_response.tsv")
print("\nwrote c8_magnet_dose_response.tsv")
