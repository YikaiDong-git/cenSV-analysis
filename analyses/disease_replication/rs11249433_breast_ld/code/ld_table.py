#!/usr/bin/env python3
"""Pearson r between each SV score (plink2 .sscore SUM columns) and the two breast SNP dosages."""
import csv
import math
import sys

score_f, sscore_f, raw_f, cov_f, out_f = sys.argv[1:6]

with open(score_f) as f:
    svs = f.readline().rstrip("\n").split("\t")[2:]
    rs_w = {}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] == "1:121538815":
            rs_w = dict(zip(svs, map(float, p[2:])))

with open(sscore_f) as f:
    rows = list(csv.reader(f, delimiter="\t"))
hdr = rows[0]
sum_cols = [hdr.index(s + "_SUM") for s in svs]
iid_i = hdr.index("IID")
sv_dose = {r[iid_i]: [float(r[i]) for i in sum_cols] for r in rows[1:]}

with open(raw_f) as f:
    rr = list(csv.reader(f, delimiter="\t"))
rh = rr[0]
snp_cols = rh[6:]
snp_dose = {r[1]: [float(x) if x != "NA" else math.nan for x in r[6:]] for r in rr[1:]}

iids = [i for i in sv_dose if i in snp_dose]
assert len(iids) == len(sv_dose) == len(snp_dose), "sample mismatch"


def corr(x, y):
    pairs = [(a, b) for a, b in zip(x, y) if not (math.isnan(a) or math.isnan(b))]
    n = len(pairs)
    mx = sum(a for a, _ in pairs) / n
    my = sum(b for _, b in pairs) / n
    sxy = sum((a - mx) * (b - my) for a, b in pairs)
    sxx = sum((a - mx) ** 2 for a, _ in pairs)
    syy = sum((b - my) ** 2 for _, b in pairs)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else math.nan


cov = {}
with open(cov_f) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        cov[r["sv_locus"]] = r

with open(out_f, "w") as f:
    f.write("sv_locus\tweight_on_rs11249433\tfrac_abs_weight_in_eur\t"
            + "\t".join(f"r_{c}\tr2_{c}" for c in snp_cols) + "\n")
    for k, s in enumerate(svs):
        x = [sv_dose[i][k] for i in iids]
        out = []
        for j in range(len(snp_cols)):
            r = corr(x, [snp_dose[i][j] for i in iids])
            out += [f"{r:.3f}", f"{r * r:.3f}"]
        f.write(f"{s}\t{rs_w.get(s, 0.0):.4f}\t{cov[s]['frac_abs_weight_in_eur']}\t" + "\t".join(out) + "\n")
print(f"n={len(iids)} samples; SNP columns: {snp_cols}")
