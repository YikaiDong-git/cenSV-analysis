#!/usr/bin/env python
import os, glob, sys
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
WU = f"{P1}/analyses/cdr_qtl_full_panel"
ROUTE = sys.argv[1] if len(sys.argv) > 1 else "per_hap_bam"
SUF = "" if ROUTE == "per_hap" else "_bam"
RDIR = "bam" if SUF else "chain"

IN_DIR = f"{WU}/output/mqtl/{ROUTE}"
OUT = f"{WU}/output/mqtl"
WINDOWS = [10_000, 50_000, 100_000]

files = sorted(glob.glob(f"{IN_DIR}/*.tsv.gz"))
print(f"per-haplotype files: {len(files)}")
D = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
D["sample"] = D["tag"].str.rsplit("_", n=1).str[0]
print(f"rows {len(D):,}   haplotypes {D['tag'].nunique()}   individuals {D['sample'].nunique()}")

rows = []
for W in WINDOWS:
    d = D[D.window_bp == W]
    ind = d.groupby(["sv_locus", "chrom", "sample"], observed=True).agg(
        meth=("meth", "mean"), covg=("covg", "mean"), n_bins=("n_bins", "sum")).reset_index()
    ind = ind[ind.meth.notna()]
    Mth = ind.pivot_table(index="sv_locus", columns="sample", values="meth")
    Cov = ind.pivot_table(index="sv_locus", columns="sample", values="covg")
    Mth.to_csv(f"{OUT}/matrices/{RDIR}/matrix_meth_{W//1000}kb{SUF}.tsv.gz", sep="\t", float_format="%.6g",
               compression="gzip")
    Cov.to_csv(f"{OUT}/matrices/{RDIR}/matrix_cov_{W//1000}kb{SUF}.tsv.gz", sep="\t", float_format="%.6g",
               compression="gzip")
    n_ind = Mth.notna().sum(axis=1)
    for sv, n in n_ind.items():
        rows.append(dict(window_bp=W, sv_locus=sv, n_individuals=int(n)))
    tot_sv = D.sv_locus.nunique()
    print(f"\n=== window +/-{W//1000} kb ===")
    print(f"  SVs with >=1 individual : {len(n_ind):,} / {tot_sv:,} ({100*len(n_ind)/tot_sv:.1f}%)")
    for k in (30, 50, 100, 150):
        c = int((n_ind >= k).sum())
        print(f"  SVs with >={k:3d} individuals: {c:6,}  ({100*c/tot_sv:.1f}% of catalogue)")
    if len(n_ind):
        print(f"  individuals per SV: median {int(n_ind.median())}  "
              f"p25 {int(n_ind.quantile(.25))}  p75 {int(n_ind.quantile(.75))}  max {int(n_ind.max())}")

C = pd.DataFrame(rows)
C.to_csv(f"{OUT}/sv_coverage{SUF}.tsv", sep="\t", index=False)
print(f"\nwrote {OUT}/sv_coverage{SUF}.tsv and matrices/matrix_*{SUF}.tsv.gz  [route={ROUTE}]")

print("\n=== per-chromosome testable SVs at +/-100 kb (>=50 individuals) ===")
d100 = C[C.window_bp == 100_000].merge(
    D[["sv_locus", "chrom"]].drop_duplicates(), on="sv_locus", how="left")
g = d100.groupby("chrom").agg(sv_with_data=("sv_locus", "size"),
                              sv_ge50=("n_individuals", lambda s: int((s >= 50).sum())))
allsv = D[D.window_bp == 100_000].groupby("chrom").sv_locus.nunique().rename("sv_total")
print(g.join(allsv).fillna(0).astype(int).to_string())
