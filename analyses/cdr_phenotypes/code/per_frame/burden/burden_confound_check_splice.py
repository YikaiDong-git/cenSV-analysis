#!/usr/bin/env python
"""STEP 0 of the many-SV question: is \"SV burden\" substantially an assembly-quality readout?"""
import glob, os, re, sys
import numpy as np, pandas as pd
from collections import defaultdict
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
AUD = f"{ROOT}/analyses/cdr_calls"
CAN = CENSV_ROOT + "/analyses/cdr_phenotypes/output/compat/compat_root_splice"
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
OUT = f"{CENSV_ROOT}/analyses/cdr_phenotypes/output/per_frame/recompute/confound_splice"
os.makedirs(OUT, exist_ok=True)

G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
samples = np.array(open(f"{WU}/geno/full_samples.txt").read().split())
print(f"genotypes {G.shape} over {len(samples)} individuals, {len(vt)} loci")

rep = vt.reset_index().sort_values("index").drop_duplicates("geno_group", keep="first")
ri = rep["index"].values
Gu, vu = G[ri, :], vt.iloc[ri].reset_index(drop=True)
print(f"deduplicated to {len(vu)} genotype-groups for burden")

sign = np.where(vu.svtype.str.startswith("DEL"), -1.0, 1.0)
burden_n = np.nansum(Gu, axis=0)
burden_bp = np.nansum(Gu * (vu.svlen.values * sign)[:, None], axis=0)
burden_ch = np.nansum(Gu * np.abs(vu.svlen.values)[:, None], axis=0)
B = pd.DataFrame(dict(sample=samples, burden_n=burden_n, burden_bp=burden_bp,
                      burden_churn=burden_ch)).set_index("sample")

percc = {}
for c in sorted(set(vu.chrom)):
    m = (vu.chrom == c).values
    percc[c] = pd.DataFrame(dict(sample=samples, chrom=c,
                                 burden_n=np.nansum(Gu[m, :], axis=0),
                                 burden_churn=np.nansum(
                                     Gu[m, :] * np.abs(vu.svlen.values[m])[:, None], axis=0)))
PC = pd.concat(percc.values(), ignore_index=True)

coh = pd.read_csv(f"{AUD}/output/cohort_corrected_final.tsv", sep="\t")
print(f"cohort rows {len(coh)}  ({coh['sample'].nunique()} samples x hap)")

cov_rows = []
for f in sorted(glob.glob(f"{AUD}/output/calls/*/binned_freq.bed")):
    tag = os.path.basename(os.path.dirname(f))
    s, h = tag.rsplit("_", 1)
    try:
        a = pd.read_csv(f, sep="\t", header=None, usecols=[3, 4], names=["pct", "cov"])
    except Exception:
        continue
    cov_rows.append(dict(sample=s, hap=h, coverage=a["cov"].mean(), n_windows=len(a)))
CV = pd.DataFrame(cov_rows)
print(f"coverage parsed for {len(CV)} (sample,hap) pairs")

Q = coh.merge(CV, on=["sample", "hap"], how="outer")
Q["is_r1041"] = (Q.chemistry.astype(str).str.contains("1041")).astype(float)
qi = Q.groupby("sample").agg(n_contigs=("n_contigs", "sum"), coverage=("coverage", "mean"),
                             n_windows=("n_windows", "sum"), is_r1041=("is_r1041", "max"),
                             mean_methyl=("mean_methyl", "mean")).reset_index()

ph = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
ind_phen = ph.groupby("sample").agg(mean_5mC=("mean_5mC", "mean"),
                                    active_len=("active_len", "mean")).reset_index()

D = (B.reset_index().merge(qi, on="sample", how="inner")
     .merge(ind_phen, on="sample", how="inner")).dropna(
    subset=["burden_n", "n_contigs", "coverage", "is_r1041", "mean_5mC"])
print(f"\ncomplete individuals for the confound check: {len(D)}")


def r2(y, X, names):
    X = np.column_stack([np.ones(len(y))] + [np.asarray(x, float) for x in X])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    R2 = 1 - res.var() / y.var()
    n, k = len(y), X.shape[1]
    adj = 1 - (1 - R2) * (n - 1) / (n - k)
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X) * (res ** 2).sum() / (n - k)))
    t = b / se
    return R2, adj, dict(zip(names, zip(b[1:], t[1:])))


COV = [D.n_contigs.values, D.coverage.values, D.is_r1041.values]
CN = ["n_contigs", "coverage", "is_r1041"]
print("\n" + "=" * 88)
print("Q1  How much of SV BURDEN is explained by assembly quality?")
for lab, y in [("burden_n (alt-allele count)", D.burden_n.values),
               ("burden_churn (unsigned bp)", D.burden_churn.values),
               ("burden_bp (net bp)", D.burden_bp.values)]:
    R2, adj, terms = r2(y, COV, CN)
    print(f"  {lab:30s} R2={R2:.4f} (adj {adj:.4f})   "
          + "  ".join(f"{k} t={v[1]:+.2f}" for k, v in terms.items()))

print("\n" + "=" * 88)
print("Q2  How much of the PHENOTYPE is explained by assembly quality?")
for lab, y in [("mean_5mC", D.mean_5mC.values), ("active_len", D.active_len.values)]:
    R2, adj, terms = r2(y, COV, CN)
    print(f"  {lab:30s} R2={R2:.4f} (adj {adj:.4f})   "
          + "  ".join(f"{k} t={v[1]:+.2f}" for k, v in terms.items()))

print("\n" + "=" * 88)
print("Q3  Raw burden-phenotype association, before and after adjustment")
for bl, bv in [("burden_n", D.burden_n.values), ("burden_churn", D.burden_churn.values)]:
    for pl, pv in [("mean_5mC", D.mean_5mC.values), ("active_len", D.active_len.values)]:
        r_raw, p_raw = stats.pearsonr(bv, pv)
        R2f, _, tf = r2(pv, [bv] + COV, [bl] + CN)
        b_, t_ = tf[bl]
        df = len(pv) - 5
        p_adj = 2 * (1 - stats.t.cdf(abs(t_), df))
        print(f"  {bl:14s} -> {pl:12s} raw r={r_raw:+.3f} p={p_raw:.3g}   "
              f"adjusted t={t_:+.2f} p={p_adj:.3g}")

print("\n" + "=" * 88)
print("Q4  Local (per hap-chromosome) fragmentation vs local burden")
ROW_RE = re.compile(r"^(.+?)_(chr[0-9XY]+)_((.+?)#([12])#[^:]+):")
nint = defaultdict(int)
for bed in sorted(glob.glob(f"{AUD}/inputs/globe_aa/*_AS-HOR_stv_row.ort.bed")):
    for line in open(bed):
        p = line.split("\t")
        m = ROW_RE.match(p[0])
        if m:
            nint[(m.group(1), m.group(2))] += 1
NI = pd.DataFrame([dict(sample=k[0], chrom=k[1], n_array_intervals=v) for k, v in nint.items()])
L = PC.merge(NI, on=["sample", "chrom"], how="inner").merge(
    ph.groupby(["sample", "chm13_chrom"]).agg(
        mean_5mC=("mean_5mC", "mean"), active_len=("active_len", "mean")).reset_index()
      .rename(columns={"chm13_chrom": "chrom"}), on=["sample", "chrom"], how="inner").dropna()
print(f"  hap-chromosome-level rows: {len(L)} over {L.chrom.nunique()} chromosomes")
for a_, b_ in [("burden_n", "n_array_intervals"), ("burden_churn", "n_array_intervals"),
               ("mean_5mC", "n_array_intervals"), ("burden_n", "active_len"),
               ("burden_n", "mean_5mC")]:
    r_, p_ = stats.pearsonr(L[a_], L[b_])
    print(f"  r({a_:14s}, {b_:18s}) = {r_:+.3f}   p={p_:.2g}")
print("\n  within-chromosome (removes between-centromere scale differences):")
for a_, b_ in [("burden_n", "n_array_intervals"), ("burden_n", "mean_5mC")]:
    rs = [stats.pearsonr(g[a_], g[b_])[0] for _, g in L.groupby("chrom") if len(g) > 30]
    print(f"  r({a_}, {b_}): median {np.median(rs):+.3f}  range {min(rs):+.3f}..{max(rs):+.3f}")

D.to_csv(f"{OUT}/individual_burden_and_quality.tsv", sep="\t", index=False, float_format="%.6g")
L.to_csv(f"{OUT}/perchrom_burden_and_quality.tsv", sep="\t", index=False, float_format="%.6g")
print(f"\nwrote {OUT}/individual_burden_and_quality.tsv and perchrom_burden_and_quality.tsv")
