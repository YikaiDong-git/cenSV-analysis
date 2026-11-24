#!/usr/bin/env python3
"""A1 step 3: pool Truvari summary.json across the 69 bundles, per config, overall and per chromosome."""
import json, glob, os, sys, csv
from collections import defaultdict

O = sys.argv[1]
ctrl = {r["bundle"]: (int(r["calls_ge50"]), int(r["calls_ge50_le50k"]))
        for r in csv.DictReader(open(O + "/control_calls.tsv"), delimiter="\t")}
rows = []
for cfg in ("hg002", "seq0", "seq0_big"):
    per = defaultdict(lambda: defaultdict(int))
    nb = 0
    for sj in sorted(glob.glob("%s/bench/%s/*/summary.json" % (O, cfg))):
        bid = sj.split("/")[-2]
        s = json.load(open(sj))
        chrom = bid.split("_")[0]
        for k in ("TP-base", "TP-comp", "FP", "FN"):
            per[chrom][k] += int(s[k])
            per["ALL"][k] += int(s[k])
        cf = ctrl[bid][1] if cfg != "seq0_big" else ctrl[bid][0]
        per[chrom]["FP_ctrl"] += cf
        per["ALL"]["FP_ctrl"] += cf
        nb += 1
    assert nb == 69, (cfg, nb)
    for chrom, v in per.items():
        tpb, tpc, fp, fn, fpc = v["TP-base"], v["TP-comp"], v["FP"], v["FN"], v["FP_ctrl"]
        prec = tpc / (tpc + fp) if tpc + fp else float("nan")
        prec_c = tpc / (tpc + fp + fpc) if tpc + fp + fpc else float("nan")
        rec = tpb / (tpb + fn) if tpb + fn else float("nan")
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else float("nan")
        rows.append(dict(config=cfg, chrom=chrom, tp_base=tpb, tp_comp=tpc, fp=fp, fn=fn,
                         fp_controls=fpc, precision=round(prec, 4),
                         precision_with_controls=round(prec_c, 4), recall=round(rec, 4),
                         f1=round(f1, 4)))
with open(O + "/truvari_pooled.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(sorted(rows, key=lambda r: (r["config"], r["chrom"] != "ALL", r["chrom"])))
for r in rows:
    if r["chrom"] == "ALL":
        print(r)
