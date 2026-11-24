#!/usr/bin/env python3
"""A3b: alignment-free local length test for whole-unit SVs."""
import sys, os, glob, csv, json
import numpy as np
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

tid, chrom, pos, svtype, L, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), sys.argv[6]
P = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
K, W, FL, R = 31, 10000, 5000, 200000
s0 = pos
e0 = pos + L if svtype == "DEL" else pos
fa = pysam.FastaFile(P + "/resources/chm13v2.0.fa.gz")
rlo, rhi = max(0, s0 - R), e0 + R
ref = fa.fetch(chrom, rlo, rhi).upper()


def unique_kmers(lo, hi):
    out = []
    for i in range(lo - rlo, hi - rlo - K):
        k = ref[i:i + K]
        if "N" in k:
            continue
        if ref.count(k) == 1:
            out.append((i, k))
    return out


left = unique_kmers(s0 - W - FL, s0 - W)
right = unique_kmers(e0 + W, e0 + W + FL)
assert left and right, "no unique anchors for %s (left %d, right %d)" % (tid, len(left), len(right))
D_ref = right[0][0] - left[-1][0]


def query_slice(rec, a, b):
    a, b = max(a, rec.reference_start), min(b, rec.reference_end)
    r, q = rec.reference_start, 0
    qa = qb = None
    for op, n in rec.cigartuples:
        if op == 4:
            q += n; continue
        if op == 5:
            continue
        if op in (0, 7, 8):
            if qa is None and r <= a < r + n:
                qa = q + (a - r)
            if qb is None and r < b <= r + n:
                qb = q + (b - r)
            r += n; q += n
        elif op == 1:
            q += n
        elif op in (2, 3):
            if qa is None and r <= a < r + n:
                qa = q
            if qb is None and r < b <= r + n:
                qb = q
            r += n
    if qa is None or qb is None:
        return None
    return rec.query_sequence[qa:qb].upper()


def svlen(rec):
    v = rec.info.get("SVLEN")
    return abs(v[0] if isinstance(v, tuple) else v) if v is not None else abs(len(rec.alts[0]) - len(rec.ref))


tol = max(10, int(0.05 * L))
haps = sorted(os.path.basename(b)[:-len("_to_chm13.bam")] for b in glob.glob(P + "/alignments/svim_asm/*_to_chm13.bam"))
rows = []
for h in haps:
    vf = pysam.VariantFile("%s/variant_calling/svim_asm/%s_%s.vcf.gz" % (P, h, SUF))
    car = 0
    if chrom in vf.header.contigs:
        for rec in vf.fetch(chrom, max(0, pos - 101), pos + 100):
            if rec.info.get("SVTYPE") == svtype and abs(svlen(rec) - L) <= tol and abs(rec.pos - pos) <= 100:
                car = 1
    vf.close()
    bam = pysam.AlignmentFile("%s/alignments/svim_asm/%s_to_chm13.bam" % (P, h))
    seq = None
    for rec in bam.fetch(chrom, s0, e0 + 1):
        if rec.is_unmapped or rec.is_secondary or rec.mapping_quality < 20 or rec.query_sequence is None:
            continue
        if rec.reference_start <= s0 - W - FL and rec.reference_end >= e0 + W + FL:
            seq = query_slice(rec, rlo, rhi)
            break
    bam.close()
    row = dict(haplotype=h, svim_carrier=car, covered=int(seq is not None))
    if seq is not None:
        lp = rp = None
        for i, k in reversed(left):
            if seq.count(k) == 1:
                lp = seq.index(k); break
        for i, k in right:
            if seq.count(k) == 1:
                rp = seq.index(k); break
        if lp is not None and rp is not None and rp > lp:
            row["delta"] = (rp - lp) - D_ref
    rows.append(row)

with open("%s/%s_unitcount.tsv" % (OUT, tid), "w", newline="") as fh:
    w = csv.DictWriter(fh, ["haplotype", "svim_carrier", "covered", "delta"], delimiter="\t", restval="")
    w.writeheader(); w.writerows(rows)
c = np.array([r["delta"] for r in rows if "delta" in r and r["svim_carrier"]], float)
n = np.array([r["delta"] for r in rows if "delta" in r and not r["svim_carrier"]], float)
exp = -L if svtype == "DEL" else L
summ = dict(target=tid, svtype=svtype, length=L, expected_shift=exp, n_left_anchors=len(left), n_right_anchors=len(right),
            anchor_span_ref=D_ref, carriers_measured=int(len(c)), noncarriers_measured=int(len(n)))
if len(c) and len(n):
    base = float(np.median(n))
    shift = float(np.median(c)) - base
    summ.update(median_shift=shift, shift_over_expected=round(shift / exp, 3),
                carriers_within_halfL_of_expected=round(float(np.mean(np.abs((c - base) - exp) <= L / 2)), 3),
                noncarriers_within_halfL_of_expected=round(float(np.mean(np.abs((n - base) - exp) <= L / 2)), 3),
                noncarrier_delta_iqr=[float(np.percentile(n, 25)), float(np.percentile(n, 75))])
json.dump(summ, open("%s/%s_unitcount_summary.json" % (OUT, tid), "w"), indent=1)
print(json.dumps(summ))
