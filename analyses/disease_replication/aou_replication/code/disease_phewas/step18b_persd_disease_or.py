#!/usr/bin/env python3
"""Step 18b: re-express the ENTIRE EUR disease PheWAS on a per-SD dosage axis (ANALYTIC)."""
import os, csv, math

HOME = os.path.expanduser("~")
RV = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds"
OUT = f"{RV}/step16_disease"
NUMPY = f"{OUT}/disease_phewas_vds_eur.tsv"
QC = [f"{RV}/step6_222/dosage_qc_222.tsv", f"{RV}/step2/dosage_qc.tsv"]

sd = {}
for qf in QC:
    for r in csv.DictReader(open(qf), delimiter="\t"):
        sd.setdefault(r["sv_locus"], float(r["dosage_sd"]))

rows = list(csv.DictReader(open(NUMPY), delimiter="\t"))
persd = f"{OUT}/disease_phewas_persd_eur.tsv"
n_missing_sd = 0
with open(persd, "w") as f:
    f.write("sv_locus\tdisease\tlabel\tcategory\tn\tn_cases\tOR_perSD\tp\tq_bh\tdosage_sd\n")
    for x in rows:
        b = x["beta_logOR"]
        if b in ("", "nan") or x["sv_locus"] not in sd:
            orsd = ""
            if x["sv_locus"] not in sd and b not in ("", "nan"):
                n_missing_sd += 1
        else:
            orsd = f"{math.exp(float(b) * sd[x['sv_locus']]):.4g}"
        f.write(f"{x['sv_locus']}\t{x['disease']}\t{x['label']}\t{x['category']}\t{x['n']}\t"
                f"{x['n_cases']}\t{orsd}\t{x['p']}\t{x['q_bh']}\t{sd.get(x['sv_locus'],''):}\n")
print(f"wrote {persd}  ({len(rows)} pairs; {n_missing_sd} lacked a QC dosage_sd)")

orsd_map = {}
for x in rows:
    b = x["beta_logOR"]
    if b not in ("", "nan") and x["sv_locus"] in sd:
        orsd_map[(x["sv_locus"], x["disease"])] = math.exp(float(b) * sd[x["sv_locus"]])
comb = f"{OUT}/disease_phewas_COMBINED_eur.tsv"
old = list(csv.DictReader(open(comb), delimiter="\t"))
with open(comb, "w") as f:
    f.write("sv_locus\tdisease\tcategory\tn\tOR_perSD\tp_plink2\tq_bh\tflag\n")
    for x in old:
        v = orsd_map.get((x["sv_locus"], x["disease"]))
        col = x.get("OR_perSD", x.get("OR_rescaled", ""))
        orstr = f"{v:.4g}" if v is not None else col
        f.write(f"{x['sv_locus']}\t{x['disease']}\t{x['category']}\t{x['n']}\t{orstr}\t"
                f"{x['p_plink2']}\t{x['q_bh']}\t{x['flag']}\n")
print(f"rewrote {comb} with OR_perSD (was OR_rescaled)")

sig = [x for x in old if x["q_bh"] not in ("", "nan") and float(x["q_bh"]) < 0.05]
for x in sorted(sig, key=lambda x: float(x["q_bh"])):
    v = orsd_map.get((x["sv_locus"], x["disease"]))
    print(f"  {x['sv_locus'].split('_')[0]}_{x['sv_locus'].split('_')[1]}  {x['disease']:<28} OR/SD={v:.3f}  p={x['p_plink2']}  q={x['q_bh']}  {x['flag']}")
