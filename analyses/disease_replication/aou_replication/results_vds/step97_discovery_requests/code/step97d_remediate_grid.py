#!/usr/bin/env python3
"""step97d - remediate grid_255cells_6strata.tsv before external release."""
import csv, math, statistics

SRC = "grid_255cells_6strata.tsv"
DST = "grid_255cells_6strata.tsv"
NA = ("", "NA", "nan", "NaN")

if open(SRC).readline().startswith("#"):
    raise SystemExit("ABORT: grid is already remediated. This script renames its own input columns and is "
                     "NOT idempotent - re-run it only against a freshly generated step97a grid.")
rows = list(csv.DictReader(open(SRC), delimiter="\t"))

def f(v):
    return None if v in NA else float(v)

def norm_sf2(z):
    return math.erfc(abs(z) / math.sqrt(2.0))

worst, ratios, flips = 0.0, [], []
for r in rows:
    bf, se, pf = f(r["beta_firth_logOR"]), f(r["se_score"]), f(r["p_firth"])
    src = (r.get("firth_source_full") or "").strip()
    if bf is None or se is None or pf is None or not src or se <= 0 or pf <= 0:
        continue
    naive = norm_sf2(bf / se)
    if naive <= 0:
        continue
    ratio = max(naive / pf, pf / naive)
    ratios.append(ratio)
    worst = max(worst, ratio)
    if (pf < 0.05) != (naive < 0.05):
        flips.append((r["sv_locus"], r["phecode"], r["stratum"], pf, naive))

n_firth = len(ratios)
gt2 = sum(1 for x in ratios if x > 2)
gt10 = sum(1 for x in ratios if x > 10)

HDR = f"""# grid_255cells_6strata.tsv - unthresholded per-ancestry association grid
# 255 SV x phecode cells (q<0.05 union 0.05<=q<0.10 from the MCC=2 ext2 arm) x 6 ancestry strata = 1530 rows.
# Aggregate association statistics only. No participant-level data.
#
# *** READ THIS BEFORE FORMING A CONFIDENCE INTERVAL OR RUNNING A META-ANALYSIS ***
# se_score pairs ONLY with beta_score_logOR. Both come from the Rao score test, where V is the Fisher
# information for the cell, beta = U/V and SE = 1/sqrt(V). That (beta, SE) pair is internally coherent and
# is the one to meta-analyse.
# There is NO standard error for beta_final_logOR, and none can be back-derived. beta_final_logOR is the
# mixed final estimator (exact Firth where it converged, approx Firth next, score otherwise) and p_final is
# its penalised-likelihood p. A Wald interval built from beta_final_logOR +/- 1.96*se_score is NOT valid,
# most severely for the separated cells that Firth exists to handle.
# Measured on the {n_firth} rows where a genuine Firth refit sits beside a score SE: the naive Wald p differs
# from p_final by up to {worst:,.0f}x; {gt2} of {n_firth} differ by more than 2x and {gt10} by more than 10x;
# {len(flips)} cell(s) cross the 0.05 boundary. Errors run in both directions.
#
# COLUMNS
#   cases/controls  suppressed as <20 / SUPPRESSED where cases < 20 (both, so n - cases cannot recover it)
#   beta_score_logOR, se_score, z_score, p_score   Rao score test. Coherent triple. Meta-analyse these.
#   beta_final_logOR, p_final                      mixed estimator; see final_estimator for which per row
#   final_estimator exact_firth | approx_firth | score   (blank where no estimate exists)
#   na_reason       suppressed_lt20_cases  = fewer than 20 cases, statistics withheld
#                   not_tested_no_information = >=20 cases but no finite statistic (the cell carries no
#                                               information in this stratum - this is an ANSWER, not a gap)
#   exported_in_prior_release  1 if this cell appears in the P<=1e-4 release you already hold. NOTE the
#                   shipped rule is the UNION (q_bh_tier < 0.10) OR (p_final < 1e-4), not a bare p cut.
#   cv_r2 is NOT in this file: it is a single GLOBAL number per SV from your own elastic-net CV and is not
#   per-stratum. Per-stratum imputation R2 does not exist on our side. Use per-stratum dosage variance in
#   perstratum_dosage_denominators.tsv, and se_score itself, to separate "no effect" from "no information".
"""

OUTCOLS = ["sv_locus","tier","phecode","label","stratum","cases","controls","n",
           "beta_score_logOR","se_score","z_score","p_score",
           "beta_final_logOR","p_final","final_estimator","na_reason","exported_in_prior_release"]

n_supp = n_uninf = 0
with open(DST, "w") as fh:
    fh.write(HDR)
    fh.write("\t".join(OUTCOLS) + "\n")
    for r in rows:
        src = (r.get("firth_source_full") or "").strip()
        has_est = r["beta_final_logOR" if "beta_final_logOR" in r else "beta_firth_logOR"] not in NA
        est = src if src else ("score" if has_est else "")
        if r["cases"] == "<20":
            na = "suppressed_lt20_cases"; n_supp += 1
        elif r["beta_score_logOR"] in NA:
            na = "not_tested_no_information"; n_uninf += 1
        else:
            na = ""
        out = {c: r.get(c, "") for c in OUTCOLS}
        out["beta_final_logOR"] = r["beta_firth_logOR"]
        out["p_final"] = r["p_firth"]
        out["final_estimator"] = est
        out["na_reason"] = na
        out["exported_in_prior_release"] = r["exported_in_p1e4_release"]
        fh.write("\t".join(str(out[c]) for c in OUTCOLS) + "\n")

print(f"naive-Wald damage on {n_firth} genuine-Firth rows: worst {worst:,.0f}x, >2x {gt2}, >10x {gt10}, flips {len(flips)}")
for x in flips:
    print(f"  FLIP {x[0]} x {x[1]} x {x[2]}: p_final={x[3]:.3g} naive={x[4]:.3g}")
print(f"na_reason: suppressed {n_supp}, uninformative {n_uninf}, informative {len(rows)-n_supp-n_uninf}")
print(f"estimator mix: " + str({k: sum(1 for r in rows if ((r.get('firth_source_full') or '').strip() or ('score' if r['beta_firth_logOR'] not in NA else '')) == k) for k in ('exact_firth','approx_firth','score','')}))
