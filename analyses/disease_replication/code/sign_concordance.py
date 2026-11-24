#!/usr/bin/env python
import csv
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BASE = CENSV_ROOT + "/analyses/disease_replication"
DISC = f"{BASE}/aou_handoff/discovery_sign_table.tsv"
AOU  = f"{BASE}/aou_replication/results/step4/aou_results.tsv"
COV  = f"{BASE}/aou_replication/results/step1/step1_coverage.tsv"
OUT  = f"{BASE}/aou_replication/results/step4/sign_concordance.tsv"
NTEST = 28
BONF  = 0.05 / NTEST

def sgn(x):
    try: x = float(x)
    except: return "NA"
    return "+" if x > 0 else ("-" if x < 0 else "0")

disc = {r["sv_locus"]: r for r in csv.DictReader(open(DISC), delimiter="\t")}
cov  = {r["sv_locus"]: r for r in csv.DictReader(open(COV),  delimiter="\t")}
aou  = {}
for r in csv.DictReader(open(AOU), delimiter="\t"):
    aou.setdefault(r["sv_locus"], {})[r["ancestry_mode"]] = r

rows = []
for sv, modes in aou.items():
    d = disc.get(sv, {})
    dsign = d.get("discovery_sign", "NA")
    c = cov.get(sv, {})
    any_row = next(iter(modes.values()))
    rec = {}
    for mode in ("eur", "all"):
        m = modes.get(mode)
        if m:
            s = sgn(m["beta"])
            conc = "Y" if (s == dsign and dsign in "+-") else ("N" if dsign in "+-" else "NA")
            rec[mode] = (m["beta"], m["p"], s, conc)
        else:
            rec[mode] = ("", "", "", "")
    rows.append(dict(sv=sv, rep_trait=d.get("rep_trait", any_row["rep_trait"]),
                     trait=any_row["trait"], in_robust24=d.get("in_robust24", ""),
                     dsign=dsign, dbz=d.get("discovery_bz", ""),
                     cov=f"{c.get('n_recoverable','?')}/{c.get('n_total','?')}",
                     eur=rec["eur"], allm=rec["all"]))

cols = ["sv_locus","rep_trait","trait","in_robust24","cov_recoverable","discovery_sign","discovery_bz",
        "eur_beta","eur_p","eur_sign","eur_concordant","all_beta","all_p","all_sign","all_concordant"]
with open(OUT, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: (r["in_robust24"] != "True", r["sv"])):
        e, a = r["eur"], r["allm"]
        fh.write("\t".join(map(str, [r["sv"], r["rep_trait"], r["trait"], r["in_robust24"], r["cov"],
                 r["dsign"], r["dbz"], e[0], e[1], e[2], e[3], a[0], a[1], a[2], a[3]])) + "\n")

def tier(p):
    try: p = float(p)
    except: return "na"
    if p < BONF: return "bonf"
    if p < 0.05: return "nom"
    if p < 0.4:  return "sub"
    return "null"

def tag(r, key):
    return f"{'_'.join(r['sv'].split('_')[:2])}/{r['trait']}({r['dsign']}vs{r[key][2]},p={float(r[key][1]):.1e})"

print(f"Bonferroni threshold = 0.05/{NTEST} = {BONF:.2e}\n")
for mode, key in (("EUR", "eur"), ("ALL-ancestry", "allm")):
    tst = [r for r in rows if r[key][1] != ""]
    conc = sum(1 for r in tst if r[key][3] == "Y")
    print(f"===== {mode} ({len(tst)} testable) =====")
    print(f"  overall sign concordance vs discovery: {conc}/{len(tst)}")
    for t, lab in (("bonf", f"Bonferroni p<{BONF:.1e}"), ("nom", "nominal p<0.05"), ("sub", "suggestive 0.05<=p<0.4")):
        sub = [r for r in tst if tier(r[key][1]) == t]
        if sub:
            cc = sum(1 for r in sub if r[key][3] == "Y")
            print(f"  {lab}: {cc}/{len(sub)} concordant | " + "; ".join(tag(r, key) for r in sub))
    print()

watch = {"chr12_sv721_DEL_non_best6_alt_edge1000":"survivor:height",
         "chr11_sv263_DEL_non_best8_edge500":"survivor:BMI",
         "chr11_sv271_DEL_ori_best8_edge500":"survivor:albumin",
         "chr3_sv40_DEL_non_best6_alt_edge500":"POS_control:PROS1/VTE",
         "chr19_sv1949_INS_mon_best8_midpoint500":"NEG_control:SNP-equiv"}
print("===== survivors & controls (EUR) =====")
by = {r["sv"]: r for r in rows}
for sv, lab in watch.items():
    r = by.get(sv)
    if not r: print(f"  {lab:24s} {sv}: NOT TESTED in AoU"); continue
    e = r["eur"]
    print(f"  {lab:24s} disc={r['dsign']} AoU_eur_beta={e[0]} p={e[1]} sign={e[2]} concordant={e[3]}")
print(f"\nWROTE {OUT}")
print("DONE_CONCORDANCE")
