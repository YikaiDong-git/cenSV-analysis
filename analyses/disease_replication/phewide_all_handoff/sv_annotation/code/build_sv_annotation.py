#!/usr/bin/env python3
import sys, csv

SUB   = sys.argv[1]
MSIG  = sys.argv[2]
MEXT  = sys.argv[3]
OUT   = sys.argv[4]
OUT16 = sys.argv[5]

HIT16 = {
    "chr16_sv1029","chr16_sv559","chr16_sv37","chr16_sv1017",
    "chr11_sv200","chr11_sv400","chr11_sv148","chr11_sv389",
    "chr1_sv650","chr1_sv2028","chr10_sv697","chr10_sv698",
    "chr17_sv317","chr17_sv331",
}

def load_set(path):
    s = set()
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            s.add(row["sv_locus"])
    return s

sig = load_set(MSIG)
ext = load_set(MEXT)

def sv_key(sv_locus):
    p = sv_locus.split("_")
    return f"{p[0]}_{p[1]}"

OUT_COLS = ["sv_locus","sv_key","svtype","cv_r2","n_svim_primitives","strategy",
            "sv_category","arm","block_chrom","block_start","block_end","block_span_mb",
            "n_weight_snps","discovery_set"]

rows = []
with open(SUB) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if row["status"] != "ADMIT":
            continue
        sl = row["sv_locus"]
        svtype = sl.split("_")[2]
        nprim = len([t for t in row["sv_id"].split(";") if t.strip()])
        bs, be = int(row["primary_block_start"]), int(row["primary_block_end"])
        dset = "signal" if sl in sig else ("null_ext" if sl in ext else "admit_only")
        rows.append({
            "sv_locus": sl,
            "sv_key": sv_key(sl),
            "svtype": svtype,
            "cv_r2": row["cv_r2"],
            "n_svim_primitives": nprim,
            "strategy": row["strategy"],
            "sv_category": row["sv_category"],
            "arm": row["arm"],
            "block_chrom": "chr"+row["primary_chrom"] if not row["primary_chrom"].startswith("chr") else row["primary_chrom"],
            "block_start": bs,
            "block_end": be,
            "block_span_mb": round((be-bs)/1e6, 3),
            "n_weight_snps": row["n_weight_snps"],
            "discovery_set": dset,
        })

def write(path, rr):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        for x in rr:
            w.writerow(x)

write(OUT, rows)
sub16 = [x for x in rows if x["sv_key"] in HIT16]
write(OUT16, sub16)

nsig = sum(1 for x in rows if x["discovery_set"]=="signal")
next = sum(1 for x in rows if x["discovery_set"]=="null_ext")
nao  = sum(1 for x in rows if x["discovery_set"]=="admit_only")
print(f"total ADMIT rows written : {len(rows)}")
print(f"  signal / null_ext / admit_only : {nsig} / {next} / {nao}")
print(f"16-hit-SV subset rows    : {len(sub16)}  (keys matched: {sorted(set(x['sv_key'] for x in sub16))})")
print(f"missing hit keys         : {sorted(HIT16 - set(x['sv_key'] for x in sub16))}")
import statistics as st
cvs = [float(x['cv_r2']) for x in rows]
nps = [x['n_svim_primitives'] for x in rows]
from collections import Counter
print(f"cv_r2  min/med/max       : {min(cvs):.3f} / {st.median(cvs):.3f} / {max(cvs):.3f}")
print(f"n_svim_primitives min/med/max: {min(nps)} / {int(st.median(nps))} / {max(nps)}")
print(f"svtype counts            : {dict(Counter(x['svtype'] for x in rows))}")
