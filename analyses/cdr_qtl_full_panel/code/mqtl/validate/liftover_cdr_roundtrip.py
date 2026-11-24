#!/usr/bin/env python
"""Validation C -- CDR COORDINATE ROUND-TRIP, assembly -> CHM13, BAM route vs chain route."""
import os, sys, gzip
from collections import Counter, defaultdict
import numpy as np, pandas as pd
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SUB = f"{P1}/subproject_methylation"
WU = f"{P1}/analyses/cdr_qtl_full_panel"
OUT = f"{WU}/output/mqtl/validation/roundtrip" + ("" if os.environ.get("MIRROR", "1") == "1" else "_nomirror")
MIN_MAPQ, MIN_BLOCK = 20, 50
os.makedirs(OUT, exist_ok=True)

TAG = sys.argv[1]
MIRROR = os.environ.get("MIRROR", "1") == "1"
SAMPLE, HAP = TAG.rsplit("_", 1)

HMAP = pd.read_csv(f"{WU}/output/mqtl/bam_hap_map.tsv", sep="\t")
row = HMAP[HMAP.meth_tag == TAG]
if not len(row) or not isinstance(row.bam_tag.iloc[0], str) or not row.bam_tag.iloc[0]:
    print(f"[SKIP] {TAG}: haplotype identity unresolved")
    sys.exit(0)
BAM_TAG, FAI_LABEL = row.bam_tag.iloc[0], str(row.fai_label.iloc[0])
AGREES = bool(row.agrees_with_manifest.iloc[0])
BAM = f"{P1}/alignments/svim_asm/{BAM_TAG}_to_chm13.bam"
PAF = f"{P1}/alignments/dipcall_paf/{BAM_TAG}_to_chm13.bam"
print(f"{TAG} -> BAM {BAM_TAG}  (manifest label {'AGREES' if AGREES else 'INVERTED'})")

qlen = {}
with gzip.open(PAF, "rt") as fh:
    for line in fh:
        f = line.split("\t", 3)
        if len(f) >= 2:
            qlen[f[0]] = int(f[1])
csize = {}
for line in open(f"{SUB}/raw/asm_fai/{SAMPLE}_{FAI_LABEL}.fa.gz.fai"):
    f = line.split("\t")
    if len(f) >= 2:
        csize[f[0]] = int(f[1])
qc, cc = Counter(qlen.values()), Counter(csize.values())
qb = {v: k for k, v in qlen.items() if qc[v] == 1}
cb = {v: k for k, v in csize.items() if cc[v] == 1}
NAME = {qb[L]: cb[L] for L in set(qb) & set(cb)}
SHORT = {full.rsplit("#", 1)[-1]: full for full in csize}
print(f"  name map {len(NAME)}   assembly contigs {len(csize)}")

T = pd.read_csv(f"{SUB}/results/population_cdrs_chm13.tsv", sep="\t")
T = T[(T["sample"] == SAMPLE) & (T.hap_local == HAP)]
if not len(T):
    print(f"[SKIP] {TAG}: no CDR rows")
    sys.exit(0)
print(f"  CDR rows: {len(T)}  contigs {T.asm_contig.nunique()}")

bam = pysam.AlignmentFile(BAM, "rb")
BLK = defaultdict(list)
for c in [x for x in bam.references if x.startswith("chr")]:
    for a in bam.fetch(c):
        if a.is_unmapped or a.mapping_quality < MIN_MAPQ or a.cigartuples is None:
            continue
        ctg = NAME.get(a.query_name)
        if ctg is None:
            continue
        r, q = a.reference_start, 0
        for op, ln in a.cigartuples:
            if op in (0, 7, 8):
                if ln >= MIN_BLOCK:
                    if a.is_reverse:
                        qs = a.infer_read_length() - (q + ln)
                        BLK[ctg].append((qs, qs + ln, c, r, r + ln, True))
                    else:
                        BLK[ctg].append((q, q + ln, c, r, r + ln, False))
                r += ln; q += ln
            else:
                if op in {0, 2, 3, 7, 8}:
                    r += ln
                if op in {0, 1, 4, 5, 7, 8}:
                    q += ln
for k in BLK:
    BLK[k].sort(key=lambda x: x[0])
print(f"  blocks: {sum(len(v) for v in BLK.values()):,} over {len(BLK)} contigs")


def lift(ctg_full, qpos):
    for q0, q1, chrom, r0, r1, rev in BLK.get(ctg_full, ()):
        if q0 > qpos:
            break
        if qpos < q1:
            if rev and MIRROR:
                return chrom, r0 + (q1 - 1 - qpos)
            return chrom, r0 + (qpos - q0)
    return None, None


rows = []
for t in T.itertuples():
    full = SHORT.get(t.asm_contig)
    if full is None:
        rows.append(dict(asm_contig=t.asm_contig, reason="contig_not_in_fai")); continue
    c0, p0 = lift(full, int(t.asm_start))
    c1, p1 = lift(full, int(t.asm_end))
    if c0 is None or c1 is None:
        rows.append(dict(asm_contig=t.asm_contig, reason="no_block")); continue
    if c0 != c1:
        rows.append(dict(asm_contig=t.asm_contig, reason=f"split_{c0}_{c1}")); continue
    a, b = sorted((p0, p1))
    rows.append(dict(asm_contig=t.asm_contig, asm_start=t.asm_start, asm_end=t.asm_end,
                     chain_chrom=t.chm13_chrom, chain_start=t.chm13_start, chain_end=t.chm13_end,
                     bam_chrom=c0, bam_start=a, bam_end=b,
                     chrom_match=(c0 == t.chm13_chrom),
                     d_start=a - t.chm13_start, d_end=b - t.chm13_end,
                     chain_width=t.chm13_end - t.chm13_start, bam_width=b - a, reason="ok"))

D = pd.DataFrame(rows)
D.insert(0, "tag", TAG)
D.insert(1, "label_agrees", AGREES)
D.insert(2, "mirror", MIRROR)
D.to_csv(f"{OUT}/{TAG}.tsv", sep="\t", index=False)

ok = D[D.reason == "ok"]
print(f"\n  lifted {len(ok)}/{len(D)}  ({100*len(ok)/max(len(D),1):.1f}%)")
if len(D) > len(ok):
    print("  not lifted: " + ", ".join(f"{k} {v}" for k, v in
                                       Counter(D[D.reason != 'ok'].reason).most_common(5)))
if len(ok):
    print(f"  chromosome agreement: {100*ok.chrom_match.mean():.2f}%")
    m = ok[ok.chrom_match]
    if len(m):
        print(f"  d_start: median {m.d_start.median():,.0f}  "
              f"|d|<=1kb {100*(m.d_start.abs()<=1000).mean():.1f}%  "
              f"|d|<=10kb {100*(m.d_start.abs()<=10000).mean():.1f}%  "
              f"p95|d| {m.d_start.abs().quantile(0.95):,.0f}")
        print(f"  width ratio bam/chain: median {(m.bam_width/m.chain_width.replace(0,np.nan)).median():.4f}")
print("\ndone")
