#!/usr/bin/env python3
import sys, glob, os
from collections import defaultdict, Counter

NEW   = sys.argv[1] if len(sys.argv) > 1 else "."
V1    = os.path.join(NEW, "manifests/full_tasks.tsv")
DIAGS = os.path.join(NEW, "output/full/*/*/*_diag.tsv")
OUT   = os.path.join(NEW, "manifests/full_tasks_v2.tsv")

N_TRAITS = 464
TARGET   = 10800.0
SAFETY   = 1.3
MIN_N    = 8
WALL     = "6:00:00"
MEM      = "6G"

def pct(vals, q):
    v = sorted(vals); n = len(v)
    if n == 1: return v[0]
    k = (n - 1) * q / 100.0; f = int(k); c = min(f + 1, n - 1)
    return v[f] + (v[c] - v[f]) * (k - f)

reg_bun, reg_zm = {}, {}
reg_total = defaultdict(int)
reg_allwin = set()
with open(V1) as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 8: continue
        mode, arm, cs, ce, region = p[3], p[4], int(p[5]), int(p[6]), p[7]
        reg_bun[region] = p[1]; reg_zm[region] = p[2]
        if mode == "hitbundle":
            reg_total[region] = max(reg_total[region], ce)
        elif mode == "allwindow":
            reg_allwin.add(region)

wall_by = defaultdict(list)
nfiles = 0
for f in glob.glob(DIAGS):
    nfiles += 1
    with open(f) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 17: continue
            if c[3] != "hitbundle": continue
            try: w = float(c[16])
            except ValueError: continue
            if w >= 0: wall_by[(c[0], c[2])].append(w)

allw = [w for v in wall_by.values() for w in v]
GLOBAL_FB = pct(allw, 95) if allw else 3.0

def perfit(region, arm):
    v = wall_by.get((region, arm))
    if v and len(v) >= MIN_N:
        return pct(v, 90)
    other = "matched" if arm == "impute" else "impute"
    vo = wall_by.get((region, other))
    if vo and len(vo) >= MIN_N:
        base = pct(vo, 90)
        return base * 2.0 if arm == "impute" else base * 0.5
    return GLOBAL_FB

hb_rows, aw_rows = [], []
for region in sorted(reg_total):
    total = reg_total[region]; bun = reg_bun[region]; zm = reg_zm[region]
    for arm in ("matched", "impute"):
        pf  = perfit(region, arm) * SAFETY
        blk = max(1, int(TARGET / (N_TRAITS * pf)))
        cs = 0
        while cs < total:
            ce = min(cs + blk, total)
            hb_rows.append(("hitbundle", bun, zm, arm, cs, ce, region, pf, blk))
            cs = ce
    if region in reg_allwin:
        for arm in ("matched", "impute"):
            aw_rows.append(("allwindow", bun, zm, arm, 0, total, region, 0.0, total))
rows = hb_rows + aw_rows
N_HB = len(hb_rows)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as out:
    out.write("task\tbundle\tzmat\tmode\tarm\tcstart\tcend\tregion\tmem\ttime\n")
    for i, (mode, bun, zm, arm, cs, ce, region, pf, blk) in enumerate(rows, 1):
        out.write(f"{i}\t{bun}\t{zm}\t{mode}\t{arm}\t{cs}\t{ce}\t{region}\t{MEM}\t{WALL}\n")

nh = sum(1 for r in rows if r[0] == "hitbundle")
na = sum(1 for r in rows if r[0] == "allwindow")
print(f"V2 TASKS: {len(rows)}  (hitbundle {nh}=lines 1..{nh}, allwindow {na}=lines {nh+1}..{len(rows)})   diag_files_calibrated={nfiles}")
print(f"ROLLING_HITBUNDLE_TASKS={nh}   ALLWINDOW_OFFSET={nh}   ALLWINDOW_N={na}")
print(f"GLOBAL_FALLBACK per-fit (p95, used where <{MIN_N} samples) = {GLOBAL_FB:.2f}s")
bad = 0
cov = defaultdict(list)
for (mode, bun, zm, arm, cs, ce, region, pf, blk) in rows:
    if mode == "hitbundle": cov[(region, arm)].append((cs, ce))
for (region, arm), segs in cov.items():
    segs.sort()
    exp = 0
    for cs, ce in segs:
        if cs != exp: bad += 1
        exp = ce
    if exp != reg_total[region]: bad += 1
print(f"COVERAGE_INVARIANT: {'OK (blocks tile [0,total) exactly)' if bad==0 else f'FAIL bad={bad}'}")
print("\nregion\tarm\tp90xSAFETY(s)\tblock_SV\tn_blocks\test_block_h@p90")
cnt = Counter()
meta = {}
for (mode, bun, zm, arm, cs, ce, region, pf, blk) in rows:
    if mode == "hitbundle":
        cnt[(region, arm)] += 1
        meta[(region, arm)] = (pf, blk)
for (region, arm) in sorted(cnt):
    pf, blk = meta[(region, arm)]
    est_h = blk * N_TRAITS * (pf / SAFETY) / 3600.0
    print(f"{region}\t{arm}\t{pf:.2f}\t{blk}\t{cnt[(region,arm)]}\t{est_h:.2f}")
