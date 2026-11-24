#!/usr/bin/env python3
"""A2 step 2: called net length versus assembled array-length change, per haplotype and array."""
import csv, glob, os, sys, json, importlib.util
import numpy as np
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

O = sys.argv[1]
P = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
spec = importlib.util.spec_from_file_location(
    "mal", CENSV_COLLAB_ROOT + "/exploration/"
           "simulate_svim_asm_accuracy/scripts/measure_array_lengths.py")
mal = importlib.util.module_from_spec(spec); spec.loader.exec_module(mal)
ARR = {**mal.main_arrays(), **mal.ARRAYS}

rows = []
for f in sorted(glob.glob(O + "/lengths/shard_*.tsv")):
    rows += list(csv.DictReader(open(f), delimiter="\t"))
assert len(rows) == 392 * 23, len(rows)

out = []
vcache = {}
for r in rows:
    hap, c = r["haplotype"], r["chrom"]
    a, b = ARR[c]
    assert int(r["array_bp"]) == b - a
    vp = "%s/variant_calling/svim_asm/%s_%s.vcf.gz" % (P, hap, SUF)
    if vp not in vcache:
        vcache = {vp: pysam.VariantFile(vp)}
    vf = vcache[vp]
    ins = dele = dup = 0
    n = {"INS": 0, "DEL": 0, "DUP:TANDEM": 0, "BND": 0, "INV": 0, "other": 0}
    big = 0
    if c in vf.header.contigs:
        for rec in vf.fetch(c, a, b):
            if not (a <= rec.pos - 1 < b):
                continue
            t = rec.info.get("SVTYPE")
            L = rec.info.get("SVLEN")
            L = abs(L[0] if isinstance(L, tuple) else L) if L is not None else 0
            if t == "INS":
                ins += L; n["INS"] += 1
            elif t == "DEL":
                dele += L; n["DEL"] += 1
            elif t == "DUP:TANDEM":
                dup += L; n["DUP:TANDEM"] += 1
            elif t in ("BND", "INV"):
                n[t] += 1
            else:
                n["other"] += 1
            if L >= 10000:
                big += 1
    called = ins + dup - dele
    d = dict(haplotype=hap, chrom=c, status=r["status"], array_bp=b - a,
             delta_bp=r["delta_bp"], left_gap=r["left_gap"], right_gap=r["right_gap"],
             called_ins_bp=ins, called_dup_bp=dup, called_del_bp=dele, called_net_bp=called,
             n_ins=n["INS"], n_del=n["DEL"], n_dup=n["DUP:TANDEM"], n_bnd=n["BND"], n_inv=n["INV"],
             n_ge10kb=big)
    if r["status"] == "spanned":
        d["residual_bp"] = float(r["delta_bp"]) - called
    out.append(d)

fields = []
for r in out:
    for k in r:
        if k not in fields:
            fields.append(k)
with open(O + "/netlength_per_haplotype.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fields, delimiter="\t", restval="")
    w.writeheader(); w.writerows(out)

summ = []
for c in sorted(set(r["chrom"] for r in out), key=lambda x: (len(x), x)):
    sp = [r for r in out if r["chrom"] == c and r["status"] == "spanned"]
    clean = [r for r in sp if float(r["left_gap"]) <= 10000 and float(r["right_gap"]) <= 10000]
    for lab, S in (("spanned", sp), ("spanned_gaps_le10kb", clean)):
        if len(S) < 3:
            summ.append(dict(chrom=c, set=lab, n=len(S))); continue
        x = np.array([float(r["delta_bp"]) for r in S]); y = np.array([r["called_net_bp"] for r in S], float)
        res = x - y
        slope = float(np.polyfit(x, y, 1)[0]) if np.std(x) > 0 else float("nan")
        rr = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else float("nan")
        summ.append(dict(chrom=c, set=lab, n=len(S), median_abs_delta_kb=round(float(np.median(np.abs(x))) / 1e3, 1),
                         median_abs_called_kb=round(float(np.median(np.abs(y))) / 1e3, 1),
                         pearson_r=round(rr, 3), slope_called_on_delta=round(slope, 3),
                         median_abs_residual_kb=round(float(np.median(np.abs(res))) / 1e3, 1),
                         frac_resid_le10kb=round(float(np.mean(np.abs(res) <= 10000)), 3),
                         frac_resid_le_10pct_of_delta=round(float(np.mean(np.abs(res) <= 0.1 * np.maximum(np.abs(x), 1e4))), 3),
                         frac_called_explains_half=round(float(np.mean(np.abs(res) <= 0.5 * np.abs(x))), 3)))
allsp = [r for r in out if r["status"] == "spanned"]
x = np.array([float(r["delta_bp"]) for r in allsp]); y = np.array([r["called_net_bp"] for r in allsp], float)
summ.append(dict(chrom="ALL", set="spanned", n=len(allsp), median_abs_delta_kb=round(float(np.median(np.abs(x))) / 1e3, 1),
                 median_abs_called_kb=round(float(np.median(np.abs(y))) / 1e3, 1),
                 pearson_r=round(float(np.corrcoef(x, y)[0, 1]), 3), slope_called_on_delta=round(float(np.polyfit(x, y, 1)[0]), 3),
                 median_abs_residual_kb=round(float(np.median(np.abs(x - y))) / 1e3, 1),
                 frac_resid_le10kb=round(float(np.mean(np.abs(x - y) <= 10000)), 3)))
keys = []
for s in summ:
    for k in s:
        if k not in keys:
            keys.append(k)
with open(O + "/netlength_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, keys, delimiter="\t", restval="")
    w.writeheader(); w.writerows(summ)
status = {}
for r in out:
    status.setdefault(r["chrom"], {}).setdefault(r["status"], 0)
    status[r["chrom"]][r["status"]] += 1
json.dump(status, open(O + "/span_status.json", "w"), indent=1)
print(summ[-1])
