#!/usr/bin/env python
"""Stage 4: escalate the cis-mQTL permutation resolution where the 1e4 scan hit its floor."""
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
W = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
SUF = sys.argv[2] if len(sys.argv) > 2 else "_bam"
RDIR = "bam" if SUF else "chain"


NPERM = int(sys.argv[3]) if len(sys.argv) > 3 else 1_000_000
CHUNK = 50_000
P_ESCALATE = float(sys.argv[4]) if len(sys.argv) > 4 else 0.005
MIN_IND = 50
RNG = np.random.default_rng(20260804)


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


S = pd.read_csv(f"{OUT}/scans/{RDIR}/scan_{W//1000}kb{SUF}.tsv", sep="\t")
todo = S[S.p_meth <= P_ESCALATE].copy()
print(f"window +/-{W//1000} kb: {len(S):,} scanned, escalating {len(todo)} SVs at "
      f"p_meth <= {P_ESCALATE} to {NPERM:,} permutations (floor {1/(NPERM+1):.2e})")
if not len(todo):
    sys.exit(0)

G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t", engine="python")
gsamp = np.array(open(f"{WU}/geno/full_samples.txt").read().split())
gidx = {s: i for i, s in enumerate(gsamp)}
locus_row = {l: i for i, l in enumerate(vt.sv_locus.values)}

Mth = pd.read_csv(f"{OUT}/matrices/{RDIR}/matrix_meth_{W//1000}kb{SUF}.tsv.gz", sep="\t", index_col=0)
Cov = pd.read_csv(f"{OUT}/matrices/{RDIR}/matrix_cov_{W//1000}kb{SUF}.tsv.gz", sep="\t", index_col=0)
coh = pd.read_csv(f"{AUD}/output/cohort_corrected_final.tsv", sep="\t")
coh["is_r1041"] = coh.chemistry.astype(str).str.contains("1041").astype(float)
QI = coh.groupby("sample").agg(n_contigs=("n_contigs", "sum"),
                               is_r1041=("is_r1041", "max")).reset_index()
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCS = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]],
                   columns=[f"gPC{i+1}" for i in range(10)], index=gids)
COV = QI.set_index("sample").join(PCS, how="inner")


def fl_vec(y, g, Z, nperm=NPERM):
    Zi = np.linalg.pinv(Z)
    ry = y - Z @ (Zi @ y)
    rg = g - Z @ (Zi @ g)
    denom = float(rg @ rg)
    if denom <= 1e-12:
        return np.nan, np.nan
    n, k = Z.shape
    dof = max(n - k - 1, 1)

    def tvec(RY):
        b = RY @ rg / denom
        E = RY - np.outer(b, rg)
        s2 = np.einsum("ij,ij->i", E, E) / dof
        with np.errstate(divide="ignore", invalid="ignore"):
            return b / np.sqrt(s2 / denom)

    t_obs = float(tvec(ry[None, :])[0])
    if not np.isfinite(t_obs):
        return np.nan, np.nan
    cnt, done = 0, 0
    while done < nperm:
        b = min(CHUNK, nperm - done)
        idx = np.argsort(RNG.random((b, n)), axis=1)
        cnt += int((np.abs(tvec(ry[idx])) >= abs(t_obs)).sum())
        done += b
    return t_obs, (cnt + 1) / (nperm + 1)


rows = []
for r in todo.itertuples():
    sv = r.sv_locus
    m = Mth.loc[sv]
    ids = [s for s in m.index if pd.notna(m[s]) and s in gidx and s in COV.index]
    gvec = G[locus_row[sv], [gidx[s] for s in ids]]
    keep = np.isfinite(gvec)
    ids = [s for s, k in zip(ids, keep) if k]
    if len(ids) < MIN_IND:
        continue
    gvec = gvec[keep]
    cvals = Cov.loc[sv, ids].values.astype(float)
    Zc = COV.loc[ids].values.astype(float)
    Z = np.column_stack([np.ones(len(ids)), Zc, rint(cvals)])
    t_m, p_m = fl_vec(rint(m[ids].values.astype(float)), gvec, Z)
    t_c, p_c = fl_vec(rint(cvals), gvec, np.column_stack([np.ones(len(ids)), Zc]))
    rows.append(dict(sv_locus=sv, n=len(ids), t_meth_1e6=t_m, p_meth_1e6=p_m,
                     t_cov_1e6=t_c, p_cov_1e6=p_c))
    print(f"  {sv:44s} n={len(ids):3d}  t={t_m:+.2f}  p={p_m:.3e}  (was {r.p_meth:.1e})",
          flush=True)

E = pd.DataFrame(rows).merge(
    S[["sv_locus", "chrom", "pos", "svtype", "svlen", "af", "p_meth", "p_cov"]],
    on="sv_locus", how="left")
alpha = 0.05 / len(S)
E["pass_bonferroni"] = E.p_meth_1e6 <= alpha
E["coverage_confounded"] = E.p_cov_1e6 <= 0.05
E = E.sort_values("p_meth_1e6")
E.to_csv(f"{OUT}/scans/{RDIR}/escalated_{W//1000}kb{SUF}" + ("" if NPERM == 1_000_000 else f"_{NPERM//1_000_000}e6") + ".tsv", sep="\t", index=False, float_format="%.6g")
print(f"\nalpha = 0.05/{len(S)} = {alpha:.3g}   new floor = {1/(NPERM+1):.2e}")
print(f"  pass Bonferroni at 1e6 perms       : {int(E.pass_bonferroni.sum())}")
print(f"  of those, coverage-clean           : {int((E.pass_bonferroni & ~E.coverage_confounded).sum())}")
print(E[["sv_locus", "chrom", "pos", "svtype", "af", "n", "t_meth_1e6", "p_meth_1e6",
         "p_cov_1e6", "pass_bonferroni", "coverage_confounded"]].head(15).to_string(index=False))
print(f"\nwrote {OUT}/escalated_{W//1000}kb.tsv")
