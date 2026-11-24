#!/usr/bin/env python3
"""A5: QTL and AoU hits restricted to the simulation's high-confidence stratum."""
import csv, sys, os, json
from collections import defaultdict
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)


def hc(svtype, svlen, freq):
    f = float(freq)
    return svtype in ("INS", "DEL") and abs(float(svlen)) <= 1000 and min(f, 1 - f) >= 0.05


g = list(csv.DictReader(open(P + "/results/analysis_qtl/gene_tables_cv01_maf005/genes_by_modality.tsv"), delimiter="\t"))
prim = [r for r in g if r["pick"] == "primary" and r["significant"] == "1"]
alt = {(r["ENSG"], r["modality"], r["scope"]): r for r in g if r["pick"] == "bestq"}
q = []
for r in prim:
    is_hc = hc(r["svtype"], r["svlen"], r["sv_carrier_freq"])
    a = alt.get((r["ENSG"], r["modality"], r["scope"]))
    a_hc = a is not None and a["variant_id"] != r["variant_id"] and hc(a["svtype"], a["svlen"], a["sv_carrier_freq"])
    same = a is not None and float(a["slope"]) * float(r["slope"]) > 0
    q.append(dict(symbol=r["symbol"], modality=r["modality"], scope=r["scope"], variant_id=r["variant_id"],
                  svtype=r["svtype"], svlen=r["svlen"], freq=r["sv_carrier_freq"], lead_hc=int(is_hc),
                  alt_variant=a["variant_id"] if a else "", alt_hc=int(a_hc), alt_same_sign=int(same) if a else "",
                  retained=int(is_hc or (a_hc and same))))
with open(OUT + "/qtl_hits_hc.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(q[0]), delimiter="\t"); w.writeheader(); w.writerows(q)
cis = [x for x in q if x["scope"] == "cis"]
qs = dict(n_hits=len(q), n_cis=len(cis),
          cis_lead_hc=sum(x["lead_hc"] for x in cis),
          cis_retained_via_alt=sum(x["retained"] and not x["lead_hc"] for x in cis),
          cis_retained=sum(x["retained"] for x in cis),
          cis_genes=len({x["symbol"] for x in cis}),
          cis_genes_retained=len({x["symbol"] for x in cis if x["retained"]}),
          lead_reason={"not_INS_DEL": sum(x["svtype"] not in ("INS", "DEL") for x in cis),
                       "gt_1kb": sum(x["svtype"] in ("INS", "DEL") and abs(float(x["svlen"])) > 1000 for x in cis),
                       "freq_lt_5pct": sum(min(float(x["freq"]), 1 - float(x["freq"])) < 0.05 for x in cis)},
          suv39h1=[x for x in q if x["symbol"] == "SUV39H1"])

cache = {}


def imp(chrom, strat, idx):
    k = (chrom, strat)
    if k not in cache:
        if chrom in ("chr15", "chr18"):
            f = ("%s/discovery/2026-09-11_regen_chr15_chr18/archive_pre_regen_2026-09-11/imputation/%s/%s/"
                 "%s_imputation_summary.tsv" % (P, strat, chrom, chrom))
        else:
            f = "%s/results/imputation/decompose/cosine_nj/%s_%s/%s/%s_imputation_summary.tsv" % (P, strat, SUF, chrom, chrom)
        cache[k] = {int(r["sv_index"]): r for r in csv.DictReader(open(f), delimiter="\t")}
    return cache[k][idx]


allrows = list(csv.DictReader(open(P + "/paper/figures/data/aou_ext2/disease_hits_final.tsv"), delimiter="\t"))
for r in allrows:
    t = r["sv_locus"].split("_")
    ir = imp(t[0], t[-1], int(t[1].replace("sv", "")))
    assert ir["svtype"] == t[2], r["sv_locus"]
    r["svtype"], r["svlen"], r["svfreq"] = ir["svtype"], ir["svlen"], ir["svfreq"]
    r["hc"] = hc(ir["svtype"], ir["svlen"], ir["svfreq"])
sig = [r for r in allrows if r["tier_label"] == "significant"]
assert len(sig) == 170
byph = defaultdict(list)
for r in allrows:
    byph[r["phecode"]].append(r)
ph = []
for code in sorted({r["phecode"] for r in sig}):
    rs = byph[code]
    s = [r for r in rs if r["tier_label"] == "significant"]
    top = min(s, key=lambda r: float(r["p_final"]))
    hcs = [r for r in rs if r["hc"]]
    besthc = min(hcs, key=lambda r: float(r["p_final"])) if hcs else None
    ph.append(dict(phecode=code, label=top["label"], n_sig=len(s), n_sig_hc=sum(r["hc"] for r in s),
                   top_locus=top["sv_locus"], top_svtype=top["svtype"], top_svlen=top["svlen"], top_freq=top["svfreq"],
                   top_hc=int(top["hc"]), top_beta=top["beta_final_logOR"], top_p=top["p_final"],
                   best_hc_locus=besthc["sv_locus"] if besthc else "", best_hc_tier=besthc["tier_label"] if besthc else "",
                   best_hc_p=besthc["p_final"] if besthc else "", best_hc_q=besthc["q_bh_tier"] if besthc else "",
                   best_hc_same_sign=int(float(besthc["beta_final_logOR"]) * float(top["beta_final_logOR"]) > 0) if besthc else ""))
with open(OUT + "/aou_phecodes_hc.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(ph[0]), delimiter="\t"); w.writeheader(); w.writerows(ph)
hc_sig = [r for r in sig if r["hc"]]
reason = dict(not_INS_DEL=sum(r["svtype"] not in ("INS", "DEL") for r in sig),
              gt_1kb=sum(r["svtype"] in ("INS", "DEL") and abs(float(r["svlen"])) > 1000 for r in sig),
              freq_lt_5pct=sum(min(float(r["svfreq"]), 1 - float(r["svfreq"])) < 0.05 for r in sig))
aous = dict(n_sig_assoc=len(sig), n_sig_assoc_hc=len(hc_sig), n_loci=len({r["sv_locus"] for r in sig}),
            n_loci_hc=len({r["sv_locus"] for r in hc_sig}), n_phecodes=len(ph),
            n_phecodes_with_sig_hc=sum(x["n_sig_hc"] > 0 for x in ph),
            n_phecodes_top_hc=sum(x["top_hc"] for x in ph),
            hc_same_sign_where_top_not_hc=sum(1 for x in ph if not x["top_hc"] and x["best_hc_same_sign"] == 1),
            n_top_not_hc=sum(1 for x in ph if not x["top_hc"]),
            reason_not_hc=reason,
            chr16_sig=sum(r["sv_locus"].startswith("chr16_") for r in sig),
            chr16_sig_hc=sum(r["sv_locus"].startswith("chr16_") and r["hc"] for r in sig))
cdr = dict(chr8_lead=dict(svtype="DEL", svlen=113, af=0.55, hc=hc("DEL", 113, 0.55)))
json.dump(dict(qtl=qs, aou=aous, cdr=cdr), open(OUT + "/summary.json", "w"), indent=1)
print(json.dumps(dict(qtl={k: v for k, v in qs.items() if k != "suv39h1"}, aou=aous), indent=1))
