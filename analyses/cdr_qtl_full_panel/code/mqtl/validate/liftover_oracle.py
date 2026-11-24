#!/usr/bin/env python
"""Validation B -- ORACLE TEST of the CIGAR walk against pysam's own reference implementation."""
import os, sys, gzip
from collections import Counter, defaultdict
import numpy as np, pandas as pd
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
WU = f"{P1}/analyses/cdr_qtl_full_panel"
OUT = f"{WU}/output/mqtl/validation"
MIN_MAPQ = 20
MIN_BLOCK = 50
os.makedirs(OUT, exist_ok=True)

TAG = sys.argv[1]
CHROM = sys.argv[2] if len(sys.argv) > 2 else "ALL"
PER_STRATUM = int(sys.argv[3]) if len(sys.argv) > 3 else 25
MAX_REFLEN = int(sys.argv[4]) if len(sys.argv) > 4 else 2_000_000

HMAP = pd.read_csv(f"{WU}/output/mqtl/bam_hap_map.tsv", sep="\t")
r_ = HMAP[HMAP.meth_tag == TAG]
BAM_TAG = r_.bam_tag.iloc[0]
BAM = f"{P1}/alignments/svim_asm/{BAM_TAG}_to_chm13.bam"
PAF = f"{P1}/alignments/dipcall_paf/{BAM_TAG}_to_chm13.bam"

qlen = {}
with gzip.open(PAF, "rt") as fh:
    for line in fh:
        f = line.split("\t", 3)
        if len(f) >= 2:
            qlen[f[0]] = int(f[1])

bam = pysam.AlignmentFile(BAM, "rb")
chroms = ([CHROM] if CHROM != "ALL"
          else [c for c in bam.references if c.startswith("chr") and c != "chrM"])

picked, seen = defaultdict(list), Counter()
for c in chroms:
    for a in bam.fetch(c):
        if a.is_unmapped or a.mapping_quality < MIN_MAPQ or a.cigartuples is None:
            continue
        if a.reference_length is None or not (1_000 <= a.reference_length <= MAX_REFLEN):
            continue
        key = ("rev" if a.is_reverse else "fwd", any(op == 5 for op, _ in a.cigartuples))
        seen[key] += 1
        if len(picked[key]) < PER_STRATUM:
            picked[key].append(a)
    if all(len(picked[k]) >= PER_STRATUM for k in
           [("fwd", True), ("fwd", False), ("rev", True), ("rev", False)]):
        break
print(f"{TAG} -> BAM {BAM_TAG}")
print("records available / sampled by stratum:")
for k in sorted(set(seen) | set(picked)):
    print(f"  strand={k[0]} hard_clip={k[1]}:  available {seen[k]:6d}  sampled {len(picked[k])}")

rows, badL = [], 0
for key, recs in picked.items():
    for a in recs:
        cig = a.cigartuples
        hard = any(op == 5 for op, _ in cig)
        L = a.infer_read_length()
        if a.query_name in qlen and qlen[a.query_name] != L:
            badL += 1
        lead_h = cig[0][1] if cig[0][0] == 5 else 0

        truth = {}
        for qp, rp in a.get_aligned_pairs(matches_only=True):
            if qp is None or rp is None:
                continue
            full = qp + lead_h
            truth[rp] = (L - 1 - full) if a.is_reverse else full

        blocks, r, q = [], a.reference_start, 0
        for op, ln in cig:
            if op in (0, 7, 8):
                if ln >= MIN_BLOCK:
                    if a.is_reverse:
                        qs = a.infer_read_length() - (q + ln)
                        blocks.append((r, r + ln, qs, qs + ln))
                    else:
                        blocks.append((r, r + ln, q, q + ln))
                r += ln; q += ln
            else:
                if op in {0, 2, 3, 7, 8}:
                    r += ln
                if op in {0, 1, 4, 5, 7, 8}:
                    q += ln

        for (r0, r1, q0, q1) in blocks:
            ln = r1 - r0
            if ln < 200:
                continue
            for frac in (0.0, 0.25, 0.5, 0.75, 0.99):
                rp = int(r0 + frac * ln)
                if rp not in truth:
                    continue
                walk = q0 + (rp - r0)
                mirror = q0 + ln - 1 - (rp - r0)
                fixed = mirror if a.is_reverse else walk
                rows.append(dict(strand="rev" if a.is_reverse else "fwd", hard_clip=hard,
                                 contig=a.query_name, ref=rp, block_len=ln, frac=frac,
                                 walk=walk, truth=truth[rp], delta=walk - truth[rp],
                                 delta_mirror=mirror - truth[rp],
                                 delta_fixed=fixed - truth[rp]))

D = pd.DataFrame(rows)
D.to_csv(f"{OUT}/oracle_{TAG}.tsv", sep="\t", index=False)
print(f"\nprobes {len(D)}   infer_read_length disagreeing with PAF: {badL}")
if len(D):
    print("\nAS SHIPPED (walk = q0 + d):")
    for keys, g in D.groupby(["strand", "hard_clip"]):
        print(f"  strand={keys[0]} hard_clip={keys[1]}:  n={len(g):6d}  "
              f"exact={100*(g.delta.abs() <= 1).mean():6.2f}%  "
              f"median|delta|={g.delta.abs().median():,.0f}  max|delta|={g.delta.abs().max():,.0f}")
    print("\nWITH THE MIRROR (walk = q0 + len - 1 - d), which is what the CHAIN route does:")
    for keys, g in D.groupby(["strand", "hard_clip"]):
        print(f"  strand={keys[0]} hard_clip={keys[1]}:  n={len(g):6d}  "
              f"exact={100*(g.delta_mirror.abs() <= 1).mean():6.2f}%  "
              f"median|delta|={g.delta_mirror.abs().median():,.0f}")
    print("\nSTRAND-CONDITIONAL COMPOSITE = what cis_mqtl_lift_bam.py ships after the fix.")
    print("This is the regression check: every stratum must be 100.00%.")
    worst = 0.0
    for keys, g in D.groupby(["strand", "hard_clip"]):
        ex = 100 * (g.delta_fixed.abs() <= 1).mean()
        worst = max(worst, float(g.delta_fixed.abs().max()))
        print(f"  strand={keys[0]} hard_clip={keys[1]}:  n={len(g):6d}  exact={ex:6.2f}%  "
              f"max|delta|={g.delta_fixed.abs().max():,.0f}")
    print(f"  OVERALL exact={100*(D.delta_fixed.abs() <= 1).mean():.4f}%   worst |delta| = {worst:,.0f} bp")
    print("\nby position within block (a DIRECTION error is exact only at the midpoint):")
    for (st, f_), g in D.groupby(["strand", "frac"]):
        print(f"  strand={st} frac={f_:.2f}:  exact={100*(g.delta.abs() <= 1).mean():6.2f}%  "
              f"median delta = {g.delta.median():+,.0f} bp  "
              f"median delta/block_len = {np.median(g.delta / g.block_len):+.3f}")
print("\ndone")
