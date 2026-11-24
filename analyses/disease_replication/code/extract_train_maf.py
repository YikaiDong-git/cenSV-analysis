#!/usr/bin/env python
import numpy as np, csv, os
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P     = CENSV_WORKSPACE + ""
SSIMP = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
BUND  = f"{SSIMP}/output/bundles"
HAND  = f"{P}/Project1_Centromere/analyses/disease_replication/aou_handoff"
OUT   = f"{HAND}/train_maf.tsv"

META = {r["sv_locus"]: r for r in csv.DictReader(open(f"{HAND}/sv_metadata.tsv"), delimiter="\t")}
W = {}
for r in csv.DictReader(open(f"{HAND}/sv_weights_long.tsv"), delimiter="\t"):
    W.setdefault(r["sv_locus"], []).append(r)

cache = {}; rows = []; missing = 0
for sv, wl in W.items():
    region = META[sv]["region"]; chrom = region.split("_")[0].replace("chr", "")
    bp = f"{BUND}/chr{chrom}/{region}_bundle.npz"
    if bp not in cache:
        cache[bp] = np.load(bp, allow_pickle=True)
    d = cache[bp]
    pos = d["snp_pos"].astype(int)
    ref = np.array([str(x).upper() for x in d["snp_ref"]])
    alt = np.array([str(x).upper() for x in d["snp_alt"]])
    maf = d["snp_maf_panel"].astype(float)
    lut = {(int(pos[i]), ref[i], alt[i]): float(maf[i]) for i in range(len(pos))}
    for r in wl:
        m = lut.get((int(r["pos"]), r["ref"].upper(), r["alt"].upper()))
        if m is None:
            missing += 1
        rows.append([sv, r["chrom"], r["pos"], r["ref"], r["alt"], r["a1_effect"], r["beta_std"],
                     ("" if m is None else f"{m:.6f}")])

rows.sort(key=lambda r: (r[0], r[1], int(r[2])))
tmp = OUT + ".tmp"
with open(tmp, "w") as fh:
    fh.write("sv_locus\tchrom\tpos\tref\talt\ta1_effect\tbeta_std\ttrain_maf_panel\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
os.replace(tmp, OUT)
print(f"rows={len(rows)} SVs={len(W)} missing_maf={missing} -> {OUT}")
print("DONE_TRAINMAF")
