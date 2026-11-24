#!/usr/bin/env python3
import collections
import csv
import os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

B = CENSV_ROOT + ""
SRC = (B + "/analyses/sv_catalog/output/ref_ukbb/admit/"
       "admit_decompose_cv01.tsv")
OUT = (B + "/analyses/disease_replication/phewide_all_handoff/"
       "excluded_predictors_2026-08-02.tsv")

rows = list(csv.DictReader(open(SRC), delimiter="\t"))
exc = [r for r in rows if r["status"] != "ADMIT"]


def tier(r):
    if r["status"] == "FAIL_chrX_no_LD_ref":
        return "A_chrX_reference_only"
    return ("B_ref_snps_present" if int(r["n_matched"]) >= 1
            else "C_zero_ref_match")


cols = ["tier", "sv_locus", "primary_chrom", "sv_index", "sv_category", "strategy",
        "arm", "cv_r2", "n_weight_snps", "n_matched", "ld_match_frac", "n_in_block",
        "n_svim_primitives", "exclusion_status"]
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(cols)
    for r in sorted(exc, key=lambda r: (tier(r), r["primary_chrom"],
                                        int(r["sv_index"]))):
        w.writerow([tier(r)] + [r[c] for c in cols[1:-2]]
                   + [len(r["sv_id"].split(";")), r["status"]])
print("wrote %s (%d rows, %.2f MB)"
      % (OUT, len(exc), os.path.getsize(OUT) / 1e6))


def q(v, p):
    s = sorted(v)
    return s[min(len(s) - 1, int(len(s) * p))]


print("\n%-24s %6s %8s %10s %10s %9s" % ("tier", "preds", "SVs", "cv_r2 med",
                                         "cv_r2>=0.3", "med snps"))
for t in sorted({tier(r) for r in exc}):
    g = [r for r in exc if tier(r) == t]
    cv = [float(r["cv_r2"]) for r in g]
    print("%-24s %6d %8d %10.3f %10d %9d"
          % (t, len(g), len({r["sv_id"] for r in g}), q(cv, 0.5),
             sum(1 for c in cv if c >= 0.3),
             q([int(r["n_weight_snps"]) for r in g], 0.5)))
cv = [float(r["cv_r2"]) for r in exc]
print("%-24s %6d %8d %10.3f %10d %9d"
      % ("TOTAL", len(exc), len({r["sv_id"] for r in exc}), q(cv, 0.5),
         sum(1 for c in cv if c >= 0.3),
         q([int(r["n_weight_snps"]) for r in exc], 0.5)))

print("\nper-chromosome, tiers A+B only (the scoreable-plausible set):")
ab = [r for r in exc if tier(r) != "C_zero_ref_match"]
c = collections.Counter(r["primary_chrom"] for r in ab)
print("  " + "  ".join("chr%s:%d" % (k, v) for k, v in
                       sorted(c.items(), key=lambda x: (len(x[0]), x[0]))))
print("  tier A+B total: %d predictors on %d distinct SVs"
      % (len(ab), len({r["sv_id"] for r in ab})))

strong = [r for r in ab if float(r["cv_r2"]) >= 0.3 and int(r["n_weight_snps"]) >= 5]
print("  of those, cv_r2>=0.3 AND >=5 weight SNPs: %d predictors on %d SVs"
      % (len(strong), len({r["sv_id"] for r in strong})))
