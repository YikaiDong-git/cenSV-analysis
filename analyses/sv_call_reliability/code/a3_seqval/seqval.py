#!/usr/bin/env python3
"""A3: sequence validation of headline SVs in every HPRC haplotype."""
import sys, os, glob, csv, json
import edlib, pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

tid, chrom, pos, svtype, L, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), sys.argv[6]
P = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
REF = pysam.FastaFile(P + "/resources/chm13v2.0.fa.gz")
F, W, WIDE, POS_TOL, K = 60, 1000, 20000, 100, 31
os.makedirs(OUT, exist_ok=True)
s0 = pos
e0 = pos + L if svtype == "DEL" else pos


def query_slice(rec, r_lo, r_hi):
    if rec.reference_start > r_lo or rec.reference_end < r_hi:
        return None
    r, q = rec.reference_start, 0
    q_lo = q_hi = None
    for op, n in rec.cigartuples:
        if op in (4,):
            q += n; continue
        if op == 5:
            continue
        if op in (0, 7, 8):
            if q_lo is None and r <= r_lo < r + n:
                q_lo = q + (r_lo - r)
            if q_hi is None and r < r_hi <= r + n:
                q_hi = q + (r_hi - r)
            r += n; q += n
        elif op == 1:
            q += n
        elif op in (2, 3):
            if q_lo is None and r <= r_lo < r + n:
                q_lo = q
            if q_hi is None and r < r_hi <= r + n:
                q_hi = q
            r += n
    if q_lo is None or q_hi is None:
        return None
    return rec.query_sequence[q_lo:q_hi]


def svlen(rec):
    v = rec.info.get("SVLEN")
    if v is None:
        return abs(len(rec.alts[0]) - len(rec.ref)) if rec.alts else 0
    return abs(v[0] if isinstance(v, tuple) else v)


haps = sorted(os.path.basename(b)[:-len("_to_chm13.bam")]
              for b in glob.glob(P + "/alignments/svim_asm/*_to_chm13.bam"))
assert len(haps) == 392, len(haps)
tol = max(10, int(0.05 * L))
calls = {}
ins_seq = None
for h in haps:
    vf = pysam.VariantFile("%s/variant_calling/svim_asm/%s_%s.vcf.gz" % (P, h, SUF))
    hit = None
    if chrom in vf.header.contigs:
        for rec in vf.fetch(chrom, max(0, pos - POS_TOL - 1), pos + POS_TOL):
            if rec.info.get("SVTYPE") == svtype and abs(svlen(rec) - L) <= tol and abs(rec.pos - pos) <= POS_TOL:
                if hit is None or abs(rec.pos - pos) < abs(hit.pos - pos):
                    hit = rec
    calls[h] = hit
    if hit is not None and svtype == "INS" and ins_seq is None and len(hit.alts[0]) > len(hit.ref):
        ins_seq = hit.alts[0][len(hit.ref):].upper()
    vf.close()
ncar = sum(v is not None for v in calls.values())
assert ncar > 0, "no SVIM-asm carrier found for %s" % tid

refseq = lambda a, b: REF.fetch(chrom, a, b).upper()
if svtype == "DEL":
    alt_t = [refseq(s0 - F, s0) + refseq(e0, e0 + F)]
    ref_t = refseq(s0 - F, e0 + F)
    junctions = [alt_t[0]]
else:
    assert ins_seq is not None and len(ins_seq) >= F, "no literal insertion sequence in carriers"
    alt_t = [refseq(s0 - F, s0) + ins_seq[:F], ins_seq[-F:] + refseq(s0, s0 + F)]
    ref_t = refseq(s0 - F, s0 + F)
    junctions = alt_t
jk = set()
for j in junctions:
    for i in range(F - K + 1, F):
        jk.add(j[i:i + K])
rck = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
jk |= {rck(k) for k in jk}
chm_wide = refseq(s0 - WIDE, e0 + WIDE)
jk_spec = {k for k in jk if k not in chm_wide}

rows = []
for h in haps:
    bam = pysam.AlignmentFile("%s/alignments/svim_asm/%s_to_chm13.bam" % (P, h))
    win = wide = None
    for rec in bam.fetch(chrom, s0 - W, e0 + W):
        if rec.is_unmapped or rec.is_secondary or rec.mapping_quality < 20 or rec.query_sequence is None:
            continue
        x = query_slice(rec, s0 - W, e0 + W)
        if x is not None:
            win = x
            wide = query_slice(rec, s0 - WIDE, e0 + WIDE)
            break
    bam.close()
    row = dict(haplotype=h, svim_carrier=int(calls[h] is not None),
               call_pos=calls[h].pos if calls[h] is not None else "",
               call_len=svlen(calls[h]) if calls[h] is not None else "", callable=int(win is not None))
    if win is not None:
        win = win.upper()
        ed_alt = min(edlib.align(t, win, mode="HW", task="distance")["editDistance"] for t in alt_t)
        ed_ref = edlib.align(ref_t, win, mode="HW", task="distance")["editDistance"]
        la, lr = len(alt_t[0]), len(ref_t)
        if ed_alt < ed_ref and ed_alt <= 0.05 * la:
            sc = "ALT"
        elif ed_ref < ed_alt and ed_ref <= 0.05 * lr:
            sc = "REF"
        else:
            sc = "ambiguous"
        wseq = (wide or win).upper()
        nk = sum(wseq.count(k) for k in jk_spec)
        row.update(ed_alt=ed_alt, ed_ref=ed_ref, seq_call=sc, junction_kmer_hits=nk,
                   wide_window=int(wide is not None))
    rows.append(row)

keys = []
for r in rows:
    for k in r:
        if k not in keys:
            keys.append(k)
with open("%s/%s_haplotypes.tsv" % (OUT, tid), "w", newline="") as fh:
    w = csv.DictWriter(fh, keys, delimiter="\t", restval=""); w.writeheader(); w.writerows(rows)

cal = [r for r in rows if r["callable"]]
def frac(a, b):
    return round(a / b, 3) if b else None
car = [r for r in cal if r["svim_carrier"]]; non = [r for r in cal if not r["svim_carrier"]]
summ = dict(target=tid, chrom=chrom, pos=pos, svtype=svtype, length=L,
            n_haplotypes=len(rows), n_callable=len(cal), svim_carriers_total=ncar,
            carriers_callable=len(car), noncarriers_callable=len(non),
            carrier_seq_ALT=sum(r["seq_call"] == "ALT" for r in car),
            carrier_seq_REF=sum(r["seq_call"] == "REF" for r in car),
            carrier_seq_ambiguous=sum(r["seq_call"] == "ambiguous" for r in car),
            noncarrier_seq_ALT=sum(r["seq_call"] == "ALT" for r in non),
            noncarrier_seq_REF=sum(r["seq_call"] == "REF" for r in non),
            noncarrier_seq_ambiguous=sum(r["seq_call"] == "ambiguous" for r in non),
            n_specific_junction_kmers=len(jk_spec),
            carriers_with_junction_kmer=sum(r["junction_kmer_hits"] > 0 for r in car),
            noncarriers_with_junction_kmer=sum(r["junction_kmer_hits"] > 0 for r in non))
dec = [r for r in cal if r["seq_call"] != "ambiguous"]
agree = sum((r["seq_call"] == "ALT") == bool(r["svim_carrier"]) for r in dec)
summ.update(concordance_decided=frac(agree, len(dec)), n_decided=len(dec),
            sensitivity_seq=frac(summ["carrier_seq_ALT"], len(car)),
            specificity_seq=frac(summ["noncarrier_seq_REF"], len(non)),
            ambiguous_frac=frac(len(cal) - len(dec), len(cal)))
json.dump(summ, open("%s/%s_summary.json" % (OUT, tid), "w"), indent=1)
print(json.dumps(summ))
