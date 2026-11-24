#!/usr/bin/env python3
"""A3 step 0: resolve the headline SVs to (chrom, pos, svtype, length) from the paper's own tables."""
import csv, sys, os
from collections import defaultdict
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
summ_cache = {}


def imp_row(chrom, strategy, idx):
    key = (chrom, strategy)
    if key not in summ_cache:
        f = "%s/results/imputation/decompose/cosine_nj/%s_%s/%s/%s_imputation_summary.tsv" % (P, strategy, SUF, chrom, chrom)
        summ_cache[key] = {int(r["sv_index"]): r for r in csv.DictReader(open(f), delimiter="\t")}
    return summ_cache[key][idx]


def parse_locus(loc):
    t = loc.split("_")
    chrom, idx, typ, strategy = t[0], int(t[1].replace("sv", "")), t[2], t[-1]
    return chrom, idx, typ, strategy


targets = [dict(target="chr8_CDR_113DEL", chrom="chr8", pos=44300221, svtype="DEL", length=113,
                source="manuscript R8b lead", freq="")]
for r in csv.DictReader(open(P + "/results/analysis_qtl/gene_tables_cv01_maf005/genes_by_modality.tsv"), delimiter="\t"):
    if r["symbol"] == "SUV39H1" and r["pick"] == "primary":
        targets.append(dict(target="chrX_SUV39H1_%s" % r["variant_id"], chrom=r["sv_chrom"], pos=int(r["svpos"]),
                            svtype=r["svtype"], length=int(r["svlen"]), source="cis-eQTL %s %s" % (r["modality"], r["variant_id"]),
                            freq=r["sv_carrier_freq"]))
        break
else:
    raise SystemExit("SUV39H1 not found")

hits = [r for r in csv.DictReader(open(P + "/paper/figures/data/aou_ext2/disease_hits_final.tsv"), delimiter="\t")
        if r["tier_label"] == "significant"]
assert len(hits) == 170, len(hits)
best = {}
loci = defaultdict(list)
for r in hits:
    chrom, idx, typ, strat = parse_locus(r["sv_locus"])
    loci[r["sv_locus"]].append(r)
    if chrom not in best or float(r["p_final"]) < float(best[chrom]["p_final"]):
        best[chrom] = r
with open(OUT + "/aou_significant_loci.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sv_locus", "chrom", "pos", "svtype", "svlen", "svfreq", "cv_r2", "n_hits", "min_p", "top_label"])
    for loc, rs in sorted(loci.items()):
        chrom, idx, typ, strat = parse_locus(loc)
        ir = imp_row(chrom, strat, idx)
        assert ir["svtype"] == typ, (loc, ir["svtype"])
        top = min(rs, key=lambda x: float(x["p_final"]))
        w.writerow([loc, chrom, ir["pos"], ir["svtype"], ir["svlen"], ir["svfreq"], ir["cv_r2"], len(rs), top["p_final"], top["label"]])
for chrom, r in sorted(best.items(), key=lambda kv: float(kv[1]["p_final"]))[:4]:
    c, idx, typ, strat = parse_locus(r["sv_locus"])
    ir = imp_row(c, strat, idx)
    if ir["svtype"] not in ("INS", "DEL"):
        continue
    targets.append(dict(target="AoU_%s_sv%d_%s" % (c, idx, typ), chrom=c, pos=int(ir["pos"]), svtype=ir["svtype"],
                        length=abs(int(float(ir["svlen"]))), source="AoU top on %s: %s p=%s" % (c, r["label"], r["p_final"]),
                        freq=ir["svfreq"]))
with open(OUT + "/targets.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(targets[0]), delimiter="\t"); w.writeheader(); w.writerows(targets)
for t in targets:
    print(t)
