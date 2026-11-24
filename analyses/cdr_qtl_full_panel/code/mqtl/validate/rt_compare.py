#!/usr/bin/env python
"""Measure the reverse-strand mirror defect on REAL intervals, mirror ON vs OFF, same CDRs."""
import glob
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")
OUT = f"{WU}/output/mqtl/validation"
rows = []
for lab, d in [("mirror_on_fixed", "roundtrip"), ("mirror_off_prefix", "roundtrip_nomirror")]:
    fs = sorted(glob.glob(f"{OUT}/{d}/*.tsv"))
    if not fs:
        print(f"{lab}: MISSING"); continue
    D = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
    ok = D[D.reason == "ok"]
    m = ok[ok.chrom_match]
    rows.append(dict(variant=lab, haps=D.tag.nunique(), rows=len(D),
                     pct_lifted=100 * len(ok) / len(D),
                     pct_chrom_match=100 * ok.chrom_match.mean(),
                     pct_exact=100 * (m.d_start == 0).mean(),
                     pct_within_1kb=100 * (m.d_start.abs() <= 1000).mean(),
                     pct_within_10kb=100 * (m.d_start.abs() <= 10000).mean(),
                     median_abs_d=m.d_start.abs().median(),
                     p99_abs_d=m.d_start.abs().quantile(0.99),
                     pct_width_exact=100 * (m.bam_width == m.chain_width).mean()))
S = pd.DataFrame(rows)
S.to_csv(f"{OUT}/roundtrip_mirror_comparison.tsv", sep="\t", index=False, float_format="%.6g")
print(S.to_string(index=False))

try:
    A = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{OUT}/roundtrip/*.tsv")])
    B = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{OUT}/roundtrip_nomirror/*.tsv")])
    k = ["tag", "asm_contig", "asm_start", "asm_end"]
    J = A[A.reason == "ok"][k + ["bam_start"]].merge(
        B[B.reason == "ok"][k + ["bam_start"]], on=k, suffixes=("_fix", "_bug"))
    moved = J.bam_start_fix != J.bam_start_bug
    print(f"\npaired CDRs {len(J):,}   moved by the fix: {int(moved.sum()):,} "
          f"({100*moved.mean():.2f}%)")
    if moved.any():
        dd = (J.bam_start_bug - J.bam_start_fix).abs()[moved]
        print(f"  displacement among those: median {dd.median():,.0f} bp   "
              f"p90 {dd.quantile(.9):,.0f}   max {dd.max():,.0f}")
except Exception as e:
    print(f"[paired comparison skipped: {e}]")
print("\ndone")
