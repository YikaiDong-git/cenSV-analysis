#!/usr/bin/env python
"""Validation A -- CELL-BY-CELL CONCORDANCE of the two liftovers, on data already on disk."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
WU = f"{P1}/analyses/cdr_qtl_full_panel"
MAT = f"{WU}/output/mqtl/matrices"
OUT = f"{WU}/output/mqtl/validation"
BAMDIR = sys.argv[1] if len(sys.argv) > 1 else "bam"
TAG = sys.argv[2] if len(sys.argv) > 2 else ""
WINDOWS = [10, 50, 100]
os.makedirs(OUT, exist_ok=True)
print(f"chain vs {BAMDIR}")


def load(path):
    D = pd.read_csv(path, sep="\t")
    return D.set_index("sv_locus")


def melt(D, name):
    M = D.stack(dropna=True).rename(name).reset_index()
    M.columns = ["sv_locus", "sample", name]
    return M


rows, per_sample, per_chrom = [], [], []
for W in WINDOWS:
    C = load(f"{MAT}/chain/matrix_meth_{W}kb.tsv.gz")
    B = load(f"{MAT}/{BAMDIR}/matrix_meth_{W}kb_bam.tsv.gz")
    shared_sv = sorted(set(C.index) & set(B.index))
    shared_s = sorted(set(C.columns) & set(B.columns))
    J = melt(C.loc[shared_sv, shared_s], "chain").merge(
        melt(B.loc[shared_sv, shared_s], "bam"), on=["sv_locus", "sample"], how="inner")
    J["chrom"] = J.sv_locus.str.split("_").str[0]
    J["delta"] = J.bam - J.chain
    if not len(J):
        continue
    r = stats.pearsonr(J.chain, J.bam)
    rho = stats.spearmanr(J.chain, J.bam)
    rows.append(dict(window_kb=W, n_cells=len(J), n_sv=J.sv_locus.nunique(),
                     n_samples=J["sample"].nunique(), pearson_r=r.statistic,
                     spearman_rho=rho.statistic, median_delta=J.delta.median(),
                     mad_delta=float(np.median(np.abs(J.delta - J.delta.median()))),
                     frac_within_1pct=float((J.delta.abs() <= 1).mean()),
                     frac_within_5pct=float((J.delta.abs() <= 5).mean())))
    for s, g in J.groupby("sample"):
        if len(g) >= 20:
            per_sample.append(dict(window_kb=W, sample=s, n=len(g),
                                   r=stats.pearsonr(g.chain, g.bam).statistic,
                                   median_delta=g.delta.median()))
    for c, g in J.groupby("chrom"):
        if len(g) >= 20:
            per_chrom.append(dict(window_kb=W, chrom=c, n=len(g),
                                  r=stats.pearsonr(g.chain, g.bam).statistic,
                                  median_delta=g.delta.median()))
    J.to_csv(f"{OUT}/concordance_cells_{W}kb{TAG}.tsv", sep="\t", index=False, float_format="%.6g")

S = pd.DataFrame(rows)
PS = pd.DataFrame(per_sample)
PC = pd.DataFrame(per_chrom)
S.to_csv(f"{OUT}/concordance_summary{TAG}.tsv", sep="\t", index=False, float_format="%.6g")
PS.to_csv(f"{OUT}/concordance_by_sample{TAG}.tsv", sep="\t", index=False, float_format="%.6g")
PC.to_csv(f"{OUT}/concordance_by_chrom{TAG}.tsv", sep="\t", index=False, float_format="%.6g")

print(S.to_string(index=False))
print()
if len(PS):
    bad = PS[PS.r < 0.5].sort_values("r")
    print(f"per-sample r: median {PS.r.median():.4f}  min {PS.r.min():.4f}  "
          f"samples with r<0.5: {bad['sample'].nunique()}/{PS['sample'].nunique()}")
    if len(bad):
        print("  WORST (a broken name map for one haplotype would look exactly like this):")
        print(bad.head(12).to_string(index=False))
print()
if len(PC):
    print("per-chromosome:")
    print(PC.sort_values(["window_kb", "r"]).to_string(index=False))

try:
    tt = []
    for W in WINDOWS:
        a = pd.read_csv(f"{WU}/output/mqtl/scans/chain/scan_{W}kb.tsv", sep="\t")
        b = pd.read_csv(f"{WU}/output/mqtl/scans/{BAMDIR}/scan_{W}kb_bam.tsv", sep="\t")
        k = ["sv_locus", "window_bp"]
        m = a[k + ["t_meth", "n"]].merge(b[k + ["t_meth", "n"]], on=k, suffixes=("_chain", "_bam"))
        m = m[(m.n_chain >= 50) & (m.n_bam >= 50)]
        if len(m) >= 10:
            tt.append(dict(window_kb=W, n_sv=len(m),
                           r=stats.pearsonr(m.t_meth_chain, m.t_meth_bam).statistic,
                           sign_agree=float((np.sign(m.t_meth_chain) == np.sign(m.t_meth_bam)).mean())))
    if tt:
        T = pd.DataFrame(tt)
        T.to_csv(f"{OUT}/concordance_tstat{TAG}.tsv", sep="\t", index=False, float_format="%.6g")
        print("\nt-statistic concordance (SVs with n>=50 on BOTH routes):")
        print(T.to_string(index=False))
except FileNotFoundError as e:
    print(f"[t-stat concordance skipped: {e}]")
print("\ndone")
