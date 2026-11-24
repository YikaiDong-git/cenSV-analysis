#!/usr/bin/env python
"""Stage 3: the centromeric cis-mQTL scan -- does an SV's dosage move 5mC in ITS OWN neighbourhood?"""
import os, sys, struct
import numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
WU = f"{P1}/analyses/cdr_qtl_full_panel"
AUD = f"{P1}/analyses/cdr_calls"
GRM = f"{P1}/subproject_methylation/results/greml/pericen_1m_all"
OUT = f"{WU}/output/mqtl"
MIN_IND = 50
NPERM = 10_000
RNG = np.random.default_rng(20260804)
W = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
SUF = sys.argv[2] if len(sys.argv) > 2 else "_bam"
RDIR = "bam" if SUF else "chain"




def rint(a):
    a = np.asarray(a, float)
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def read_grm(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]
                K[i, j] = v; K[j, i] = v
    return K, [l.split()[1] for l in open(f"{prefix}.grm.id")]


G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t", engine="python")
gsamp = np.array(open(f"{WU}/geno/full_samples.txt").read().split())
gidx = {s: i for i, s in enumerate(gsamp)}
locus_row = {l: i for i, l in enumerate(vt.sv_locus.values)}

Mth = pd.read_csv(f"{OUT}/matrices/{RDIR}/matrix_meth_{W//1000}kb{SUF}.tsv.gz", sep="\t", index_col=0)
Cov = pd.read_csv(f"{OUT}/matrices/{RDIR}/matrix_cov_{W//1000}kb{SUF}.tsv.gz", sep="\t", index_col=0)
print(f"window +/-{W//1000} kb: {Mth.shape[0]:,} SVs x {Mth.shape[1]} individuals")

coh = pd.read_csv(f"{AUD}/output/cohort_corrected_final.tsv", sep="\t")
coh["is_r1041"] = coh.chemistry.astype(str).str.contains("1041").astype(float)
QI = coh.groupby("sample").agg(n_contigs=("n_contigs", "sum"),
                               is_r1041=("is_r1041", "max")).reset_index()
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCS = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]],
                   columns=[f"gPC{i+1}" for i in range(10)], index=gids)
COV = QI.set_index("sample").join(PCS, how="inner")
print(f"covariate table: {COV.shape[0]} individuals x {COV.shape[1]} covariates")


def freedman_lane(y, g, Z, nperm=NPERM):
    Zi = np.linalg.pinv(Z)
    ry = y - Z @ (Zi @ y)
    rg = g - Z @ (Zi @ g)
    denom = float(rg @ rg)
    if denom <= 1e-12:
        return np.nan, np.nan
    n, k = Z.shape
    dof = n - k - 1

    def tstat(resid):
        b = float(rg @ resid) / denom
        e = resid - b * rg
        s2 = float(e @ e) / max(dof, 1)
        return b / np.sqrt(s2 / denom) if s2 > 0 else np.nan

    t_obs = tstat(ry)
    if not np.isfinite(t_obs):
        return np.nan, np.nan
    cnt = 0
    for _ in range(nperm):
        if abs(tstat(RNG.permutation(ry))) >= abs(t_obs):
            cnt += 1
    return t_obs, (cnt + 1) / (nperm + 1)


rows = []
n_tested = 0
for sv in Mth.index:
    if sv not in locus_row:
        continue
    m = Mth.loc[sv]
    ids = [s for s in m.index if pd.notna(m[s]) and s in gidx and s in COV.index]
    if len(ids) < MIN_IND:
        continue
    gvec = G[locus_row[sv], [gidx[s] for s in ids]]
    keep = np.isfinite(gvec)
    ids = [s for s, k in zip(ids, keep) if k]
    if len(ids) < MIN_IND:
        continue
    gvec = gvec[keep]
    if np.std(gvec) < 1e-9:
        continue
    cvals = Cov.loc[sv, ids].values.astype(float)
    Zc = COV.loc[ids].values.astype(float)
    Z = np.column_stack([np.ones(len(ids)), Zc, rint(cvals)])
    y = rint(m[ids].values.astype(float))
    t_m, p_m = freedman_lane(y, gvec, Z)
    Z0 = np.column_stack([np.ones(len(ids)), Zc])
    t_c, p_c = freedman_lane(rint(cvals), gvec, Z0)
    rows.append(dict(sv_locus=sv, window_bp=W, n=len(ids), t_meth=t_m, p_meth=p_m,
                     t_cov=t_c, p_cov=p_c))
    n_tested += 1

R = pd.DataFrame(rows)
if not len(R):
    sys.exit("no SV met the coverage threshold")
R = R.merge(vt[["sv_locus", "chrom", "pos", "svtype", "svlen", "af"]], on="sv_locus", how="left")
alpha = 0.05 / len(R)
R["pass_bonferroni"] = R.p_meth <= alpha
R["coverage_confounded"] = R.p_cov <= 0.05
R = R.sort_values("p_meth")
R.to_csv(f"{OUT}/scans/{RDIR}/scan_{W//1000}kb{SUF}.tsv", sep="\t", index=False, float_format="%.6g")

print(f"\ntested {n_tested:,} SVs at n>={MIN_IND}   alpha = 0.05/{len(R)} = {alpha:.3g}")
print(f"  p_meth <= alpha           : {int(R.pass_bonferroni.sum())}")
print(f"  of those, coverage-clean  : {int((R.pass_bonferroni & ~R.coverage_confounded).sum())}")
print(f"  SVs where dosage predicts COVERAGE (p<=0.05): {int(R.coverage_confounded.sum())} "
      f"({100*R.coverage_confounded.mean():.1f}%)  [the T2 artifact, measured here]")
print(f"  permutation floor = {1/(NPERM+1):.2e}")
print("\ntop 12 by p_meth:")
cols = ["sv_locus", "chrom", "pos", "svtype", "svlen", "af", "n", "t_meth", "p_meth", "t_cov",
        "p_cov", "coverage_confounded"]
print(R[cols].head(12).to_string(index=False))
print(f"\nwrote {OUT}/scan_{W//1000}kb.tsv")
