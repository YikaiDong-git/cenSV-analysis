#!/usr/bin/env python
import os, sys, gzip
from collections import Counter, defaultdict
import numpy as np, pandas as pd
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SUB = f"{P1}/subproject_methylation"
WU = f"{P1}/analyses/cdr_qtl_full_panel"
AUD = f"{P1}/analyses/cdr_calls"
OUT_DIR = f"{WU}/output/mqtl/per_hap_bam"
WINDOWS = [10_000, 50_000, 100_000]
MIN_MAPQ = 20
MIN_BLOCK = 50

TAG = sys.argv[1]
SAMPLE, HAP = TAG.rsplit("_", 1)
os.makedirs(OUT_DIR, exist_ok=True)

HMAP = pd.read_csv(f"{WU}/output/mqtl/bam_hap_map.tsv", sep="\t")
row = HMAP[HMAP.meth_tag == TAG]
if not len(row) or not isinstance(row.bam_tag.iloc[0], str) or not row.bam_tag.iloc[0]:
    print(f"[SKIP] {TAG}: haplotype identity unresolved (no published .fai) -- data gap, not error")
    sys.exit(0)
BAM_TAG = row.bam_tag.iloc[0]
FAI_LABEL = str(row.fai_label.iloc[0])
print(f"{TAG}: resolved by content -> BAM {BAM_TAG}, published assembly {SAMPLE}_{FAI_LABEL}"
      f"  (label {'AGREES' if BAM_TAG.rsplit('_', 1)[-1] == HAP else 'INVERTED'})")
BAM = f"{P1}/alignments/svim_asm/{BAM_TAG}_to_chm13.bam"
PAF = f"{P1}/alignments/dipcall_paf/{BAM_TAG}_to_chm13.bam"


def chain_path():
    for cand in (HAP, "mat" if HAP == "hap1" else "pat"):
        p = f"{SUB}/raw/chains_chm13/{SAMPLE}_{cand}_vs_CHM13.chain.gz"
        if os.path.exists(p):
            return p
    return None


qlen = {}
with gzip.open(PAF, "rt") as fh:
    for line in fh:
        f = line.split("\t", 3)
        if len(f) >= 2:
            qlen[f[0]] = int(f[1])
csize = {}
SRC = "fai"
fai = f"{SUB}/raw/asm_fai/{SAMPLE}_{FAI_LABEL}.fa.gz.fai"
for line in open(fai):
    f = line.split("\t")
    if len(f) >= 2:
        csize[f[0]] = int(f[1])
qc, cc = Counter(qlen.values()), Counter(csize.values())
qb = {v: k for k, v in qlen.items() if qc[v] == 1}
cb = {v: k for k, v in csize.items() if cc[v] == 1}
NAME = {qb[L]: cb[L] for L in set(qb) & set(cb)}
print(f"{TAG}: PAF contigs {len(qlen)}  assembly contigs {len(csize)} [{SRC}]  name map {len(NAME)}")

bed = f"{AUD}/output/calls/{TAG}/binned_freq.bed"
if not os.path.exists(bed):
    print(f"[SKIP] no binned_freq for {TAG}")
    sys.exit(0)
M = pd.read_csv(bed, sep="\t", header=None,
                names=["ctg_full", "start", "end", "frac", "covg"])
ARR = {c: (g["start"].values, g["end"].values, g["frac"].values, g["covg"].values)
       for c, g in (M.sort_values("start").groupby("ctg_full"))}
bed_ctgs = set(ARR)
mapped = {NAME[q] for q in NAME} & bed_ctgs
print(f"  methylation bins {len(M)} on {len(bed_ctgs)} contigs; name map covers "
      f"{len(mapped)}/{len(bed_ctgs)}")
if not mapped:
    print(f"[SKIP] {TAG}: no name map (no published .fai and no usable chain) -- data gap, not error")
    sys.exit(0)
assert len(mapped) >= 0.5 * len(bed_ctgs), (
    f"PARTIAL NAME MAP for {TAG}: only {len(mapped)}/{len(bed_ctgs)} methylation contigs resolved. "
    f"A partial join is worse than none -- it would place windows on the wrong contig.")

vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t", engine="python")
vt = vt[["sv_locus", "chrom", "pos", "end"]]
need = sorted(set(vt.chrom))
bam = pysam.AlignmentFile(BAM, "rb")
have = set(bam.references)

BLK = defaultdict(list)
CONSUME_REF = {0, 2, 3, 7, 8}
CONSUME_QRY = {0, 1, 4, 5, 7, 8}
for c in need:
    if c not in have:
        continue
    for a in bam.fetch(c):
        if a.is_unmapped or a.mapping_quality < MIN_MAPQ or a.cigartuples is None:
            continue
        ctg = NAME.get(a.query_name)
        if ctg is None or ctg not in ARR:
            continue
        r = a.reference_start
        q = 0
        for op, ln in a.cigartuples:
            if op in (0, 7, 8):
                if ln >= MIN_BLOCK:
                    if a.is_reverse:
                        qs = a.infer_read_length() - (q + ln)
                        BLK[c].append((r, r + ln, ctg, qs, qs + ln, True))
                    else:
                        BLK[c].append((r, r + ln, ctg, q, q + ln, False))
                r += ln; q += ln
            else:
                if op in CONSUME_REF:
                    r += ln
                if op in CONSUME_QRY:
                    q += ln
for c in BLK:
    BLK[c].sort(key=lambda x: x[0])
nrev = sum(1 for v in BLK.values() for b in v if b[5])
nblk_all = sum(len(v) for v in BLK.values())
bp_rev = sum(b[1] - b[0] for v in BLK.values() for b in v if b[5])
bp_all = sum(b[1] - b[0] for v in BLK.values() for b in v) or 1
print(f"  colinear blocks: {nblk_all:,} over {len(BLK)} CHM13 sequences; "
      f"REVERSE {nrev:,} blocks ({100*nrev/max(nblk_all,1):.1f}%), "
      f"{100*bp_rev/bp_all:.1f}% of aligned bp")


def window_stats(chrom, lo, hi):
    blocks = BLK.get(chrom)
    if not blocks:
        return np.nan, np.nan, 0, 0
    wsum = fsum = csum = 0.0
    nbin = nblk = 0
    for r0, r1, ctg, q0, q1, rev in blocks:
        if r0 >= hi:
            break
        if r1 <= lo:
            continue
        d0, d1 = max(r0, lo) - r0, min(r1, hi) - r0
        if rev:
            ln_ = r1 - r0
            s, e = q0 + (ln_ - d1), q0 + (ln_ - d0)
        else:
            s, e = q0 + d0, q0 + d1
        nblk += 1
        bs, be, bf, bc = ARR[ctg]
        i0, i1 = np.searchsorted(be, s, "right"), np.searchsorted(bs, e, "left")
        if i1 <= i0:
            continue
        ov = np.minimum(be[i0:i1], e) - np.maximum(bs[i0:i1], s)
        ok = ov > 0
        if not ok.any():
            continue
        w = ov[ok].astype(float)
        wsum += w.sum()
        fsum += float(np.sum(bf[i0:i1][ok] * w))
        csum += float(np.sum(bc[i0:i1][ok] * w))
        nbin += int(ok.sum())
    if wsum <= 0:
        return np.nan, np.nan, 0, nblk
    return fsum / wsum, csum / wsum, nbin, nblk


rows = []
for r in vt.itertuples():
    for W in WINDOWS:
        f, c, nbin, nblk = window_stats(r.chrom, max(0, int(r.pos) - W), int(r.end) + W)
        rows.append((r.sv_locus, r.chrom, W, f, c, nbin, nblk))
D = pd.DataFrame(rows, columns=["sv_locus", "chrom", "window_bp", "meth", "covg", "n_bins", "n_blocks"])
D.insert(0, "tag", TAG)
out = f"{OUT_DIR}/{TAG}.tsv.gz"
D.to_csv(out, sep="\t", index=False, float_format="%.6g", compression="gzip")
cover = D.groupby("window_bp")["meth"].apply(lambda s: 100 * s.notna().mean())
print(f"  wrote {out}  rows {len(D)}")
print("  % SVs with methylation, by window: " +
      "  ".join(f"{int(w)//1000}kb {v:.1f}%" for w, v in cover.items()))
