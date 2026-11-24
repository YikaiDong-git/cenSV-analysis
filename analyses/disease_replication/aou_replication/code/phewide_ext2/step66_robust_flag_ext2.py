#!/usr/bin/env python3
"""Step 66 (ext2 arm) — robust-vs-score-inflated class on the enriched ext2 disease hits (<- step30n)."""
import csv
import math
import os

RV = os.path.expanduser(
    "~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds"
)
IN = f"{RV}/step30_phewide_ext2/phewide_disease_hits_enriched.tsv"
OUT = f"{RV}/step30_phewide_ext2/phewide_disease_hits_enriched_flagged.tsv"

P_GW = 5e-8
P_WEAK = 1e-4
OR_HI, OR_LO = 10.0, 0.1
CVR2_FLOOR = 0.30


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def classify(orf, pf, conv, cvr2):
    sane_or = OR_LO < orf < OR_HI
    if (conv != 1) or (not sane_or) or (not (pf < P_WEAK)) or math.isnan(pf):
        return "score_inflated"
    if pf < P_GW:
        return "robust" if (cvr2 >= CVR2_FLOOR) else "robust_low_impute"
    return "suggestive"


def main():
    rows = list(csv.DictReader(open(IN), delimiter="\t"))
    fields = list(rows[0].keys())
    add = ["logOR_firth", "robust_class"]
    out_fields = fields + add
    counts = {}
    for r in rows:
        orf = fnum(r["OR_firth"]); pf = fnum(r["p_firth"])
        conv = int(fnum(r["firth_converged"])) if r["firth_converged"] not in ("", "NA") else 0
        cvr2 = fnum(r["cv_r2"])
        r["logOR_firth"] = f"{math.log(orf):.3f}" if orf > 0 else ""
        cls = classify(orf, pf, conv, cvr2)
        r["robust_class"] = cls
        counts[cls] = counts.get(cls, 0) + 1

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {OUT}\n")
    print(f"{len(rows)} score-BH hits -> Firth-calibrated classes:")
    for c in ("robust", "robust_low_impute", "suggestive", "score_inflated"):
        print(f"  {c:20s} {counts.get(c,0):3d}")
    robust_any = [r for r in rows if r["robust_class"].startswith("robust")]
    usv = sorted({r["sv_key"] for r in robust_any})
    uph = sorted({r["phecode"] for r in robust_any})
    print(f"\nrobust(+low_impute): {len(robust_any)} rows, {len(usv)} unique SVs, {len(uph)} unique phecodes")
    print(f"  SVs: {usv}")
    print(f"  trustworthy core (robust, cv_r2>={CVR2_FLOOR}): "
          f"{sorted({r['sv_key'] for r in rows if r['robust_class']=='robust'})}")
    print("\n-- robust / robust_low_impute rows --")
    print(f"  {'sv':<14}{'phecode':<12}{'cases':>7}{'OR':>7}{'p_firth':>10}{'cvr2':>7}  class")
    for r in sorted(robust_any, key=lambda x: fnum(x["p_firth"])):
        print(f"  {r['sv_key']:<14}{r['phecode']:<12}{r['cases']:>7}{fnum(r['OR_firth']):>7.2f}"
              f"{fnum(r['p_firth']):>10.1e}{fnum(r['cv_r2']):>7.3f}  {r['robust_class']}")


if __name__ == "__main__":
    main()
