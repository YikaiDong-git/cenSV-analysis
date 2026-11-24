#!/usr/bin/env python
"""Combine the genome-wide trans cells into one family-wise result."""
import glob, os
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")
GW = f"{WU}/output/trans/gw"
OUT = f"{WU}/output/trans"

summ = sorted(glob.glob(f"{GW}/*/*_summary.tsv"))
S = pd.concat([pd.read_csv(f, sep="\t") for f in summ], ignore_index=True)
tasks = [l.split() for l in open(f"{OUT}/gw_tasks.txt")]
ran = set(zip(S.phenotype, S.target_chrom))
missing = [t for t in tasks if (t[0], t[1]) not in ran]
print(f"cells: {len(S)} of {len(tasks)} tasks")
if missing:
    print(f"  NOT RUN ({len(missing)}) -- excluded from the family, not silently ignored:")
    for p, c in missing[:20]:
        print(f"    {p} {c}")

mx = sorted(glob.glob(f"{GW}/*/*_maxT.npy"))
vs = [np.load(f) for f in mx]
lens = {len(v) for v in vs}
assert len(lens) == 1, f"maxT vectors differ in length {lens} -- shuffles are NOT shared, cannot combine"
NPERM = lens.pop()
gmax = np.max(np.vstack(vs), axis=0)
print(f"global max-T null: {len(vs)} cells x {NPERM:,} shared permutations")
print(f"  null max|t|: median {np.median(gmax):.3f}  p95 {np.quantile(gmax,.95):.3f}  "
      f"max {gmax.max():.3f}")

M_TRANS = int(S.n_trans.sum())
print(f"\ntotal trans tests in the family: {M_TRANS:,}")
print(f"lambda across cells: median {S.lambda_trans.median():.3f}  "
      f"IQR {S.lambda_trans.quantile(.25):.3f}-{S.lambda_trans.quantile(.75):.3f}  "
      f"max {S.lambda_trans.max():.3f}  cells > 1.2: {int((S.lambda_trans > 1.2).sum())}")
print(f"candidates at p<=0.01: observed {int(S.n_cand_p01.sum()):,}  "
      f"expected {S.exp_cand_p01.sum():,.0f}  ratio {S.n_cand_p01.sum()/S.exp_cand_p01.sum():.3f}")
S.sort_values("lambda_trans", ascending=False).to_csv(f"{OUT}/gw_cells.tsv", sep="\t", index=False,
                                                     float_format="%.6g")

cand = sorted(glob.glob(f"{GW}/*/*_cand.tsv"))
C = pd.concat([pd.read_csv(f, sep="\t") for f in cand if os.path.getsize(f) > 60], ignore_index=True)
print(f"candidate rows pooled: {len(C):,}")
C["fwer"] = [(float((gmax >= a).sum()) + 1) / (NPERM + 1) for a in C.abs_t]
C = C.sort_values("p_param").reset_index(drop=True)
C["bh_q"] = np.minimum.accumulate(
    (C.p_param.values * M_TRANS / (np.arange(len(C)) + 1))[::-1])[::-1].clip(max=1.0)
C.to_csv(f"{OUT}/gw_results.tsv", sep="\t", index=False, float_format="%.6g")

sig = C[C.fwer < 0.05]
print(f"\nFWER < 0.05: {len(sig)} test(s)")
cols = ["sv_locus", "chrom", "pos", "svtype", "svlen", "af", "phenotype", "target_chrom",
        "n", "t", "p_param", "fwer", "bh_q"]
if len(sig):
    print(sig[cols].to_string(index=False))
else:
    print("  none. Strongest 12 by parametric p:")
    print(C.head(12)[cols].to_string(index=False))
print(f"\nBH q < 0.05: {int((C.bh_q < 0.05).sum())}  (secondary; discount if lambda >> 1)")

cis = sorted(glob.glob(f"{GW}/*/*_cis.tsv"))
if cis:
    CI = pd.concat([pd.read_csv(f, sep="\t") for f in cis if os.path.getsize(f) > 60],
                   ignore_index=True)
    top = CI.nlargest(8, "abs_t")
    print(f"\nCIS POSITIVE CONTROL (same code path, {len(CI):,} tests). Strongest:")
    print(top[["sv_locus", "chrom", "pos", "phenotype", "n", "t", "p_param"]].to_string(index=False))
    c8 = CI[(CI.chrom == "chr8") & (CI.pos.between(44_200_000, 44_500_000))]
    print(f"  chr8 44.2-44.5 Mb window: {len(c8)} tests, max |t| = "
          f"{c8.abs_t.max() if len(c8) else float('nan'):.3f}  "
          f"(the known leads must appear here or the scan is broken)")
print("\ndone")
