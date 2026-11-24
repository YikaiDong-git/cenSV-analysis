#!/usr/bin/env python3
"""Wide plink2 --score file for the chr1 AoU breast-cancer SV hits (ext2 run), keyed chr:pos."""
import csv
import sys

hits_f, w_f, pvar_f, out_score, out_cov = sys.argv[1:6]

svs = []
with open(hits_f) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if "breast" in r["label"].lower() and r["chrom"] == "chr1" and r["sv_locus"] not in svs:
            svs.append(r["sv_locus"])
assert svs, "no chr1 breast hits found"

win = set()
with open(pvar_f) as f:
    for line in f:
        if not line.startswith("#"):
            win.add(line.split("\t")[2])

wt = {}
a1 = {}
tot = {s: [0, 0.0] for s in svs}
cov = {s: [0, 0.0] for s in svs}
with open(w_f) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        s = r["sv_locus"]
        if s not in tot:
            continue
        vid = r["chrom"].replace("chr", "") + ":" + r["pos"]
        b = float(r["beta_std"])
        tot[s][0] += 1
        tot[s][1] += abs(b)
        if vid not in win:
            continue
        assert a1.setdefault(vid, r["a1_effect"]) == r["a1_effect"], f"a1 conflict at {vid}"
        wt.setdefault(vid, {})[s] = b
        cov[s][0] += 1
        cov[s][1] += abs(b)

with open(out_score, "w") as f:
    f.write("ID\tA1\t" + "\t".join(svs) + "\n")
    for vid in sorted(wt, key=lambda v: int(v.split(":")[1])):
        f.write(f"{vid}\t{a1[vid]}\t" + "\t".join(str(wt[vid].get(s, 0.0)) for s in svs) + "\n")

with open(out_cov, "w") as f:
    f.write("sv_locus\tn_weight_snps\tn_in_eur_window\tfrac_abs_weight_in_eur\n")
    for s in svs:
        f.write(f"{s}\t{tot[s][0]}\t{cov[s][0]}\t{cov[s][1] / tot[s][1]:.3f}\n")
print(f"{len(svs)} SVs, {len(wt)} score variants")
