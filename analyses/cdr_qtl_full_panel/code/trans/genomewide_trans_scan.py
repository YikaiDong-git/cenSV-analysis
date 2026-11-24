#!/usr/bin/env python
"""GENOME-WIDE per-SV TRANS scan -- one (phenotype x target chromosome) cell per task."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SUB = f"{P1}/subproject_methylation"
WU = f"{P1}/analyses/cdr_qtl_full_panel"
CANON = f"{P1}/analyses/cdr_qtl"

PHENO = sys.argv[1]
TARGET = sys.argv[2]
NPERM = int(sys.argv[3]) if len(sys.argv) > 3 else 100_000
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20260805
BATCH = 1000
OUT = f"{WU}/output/trans/gw/{PHENO}"
os.makedirs(OUT, exist_ok=True)


def rint(a):
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def read_grm(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    v = np.fromfile(f"{prefix}.grm.bin", dtype=np.float32)
    G = np.zeros((n, n)); k = 0
    for i in range(n):
        for j in range(i + 1):
            G[i, j] = G[j, i] = v[k]; k += 1
    return G, [l.strip().split()[1] for l in open(f"{prefix}.grm.id")]


PC = pd.read_csv(f"{CANON}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
PD_ = pd.read_csv(f"{CANON}/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t")
src = PD_ if PHENO == "cdr_pos_offset_deepest" else PC
if PHENO not in src.columns:
    sys.exit(f"[SKIP] {PHENO} not in its phenotype file")
BC = PC.groupby("sample")["basecaller"].first().to_dict()
Gm, gids = read_grm(f"{SUB}/results/greml/pericen_1m_all")
w, v = np.linalg.eigh(Gm)
PCMAP = {s: v[:, np.argsort(w)[::-1][:10]][i] for i, s in enumerate(gids)}

V = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t", engine="python")
Z = np.load(f"{WU}/geno/full_genotypes.npz")
G = Z[list(Z.keys())[0]]
FULL = open(f"{WU}/geno/full_samples.txt").read().split()
if G.shape[0] != len(FULL):
    G = G.T
G = G.astype(np.float64)
N_PANEL = len(FULL)
fpos = {s: i for i, s in enumerate(FULL)}

ind = src[src.chm13_chrom == TARGET].groupby("sample")[PHENO].mean().dropna()
samples = [s for s in ind.index if s in fpos and s in PCMAP]
if len(samples) < 40:
    sys.exit(f"[SKIP] {PHENO}/{TARGET}: only {len(samples)} individuals")
idx = np.array([fpos[s] for s in samples])
n = len(samples)

y = rint(ind.loc[samples].values)
bc = np.array([1.0 if str(BC.get(s, "")) == "R1041" else 0.0 for s in samples])
cols = [np.ones(n), np.array([PCMAP[s] for s in samples])]
if bc.std() > 1e-10:
    cols.append(bc[:, None])
X = np.column_stack(cols)
df = n - X.shape[1] - 1
Q, _ = np.linalg.qr(X)


def resid(M):
    return M - Q @ (Q.T @ M)


yr = resid(y)
yr = (yr - yr.mean()) / yr.std()

is_cis = (V.chrom == TARGET).values
Gsub = G[idx, :]
Gr = resid(Gsub)
sd = Gr.std(axis=0)
keep = sd > 1e-10
Gr = Gr[:, keep]
Gr = (Gr - Gr.mean(axis=0)) / Gr.std(axis=0)
Vk = V[keep].reset_index(drop=True)
cis_k = is_cis[keep]
trans_cols = np.where(~cis_k)[0]
print(f"{PHENO} / {TARGET}: n={n}  df={df}  SVs kept {int(keep.sum())}/{len(V)}  "
      f"trans {len(trans_cols)}  cis {int(cis_k.sum())}")


def tstat(r):
    return r * np.sqrt(df) / np.sqrt(np.maximum(1 - r ** 2, 1e-30))


r_obs = (Gr.T @ yr) / n
t_obs = tstat(r_obs)
p_par = 2 * stats.t.sf(np.abs(t_obs), df)

lam_all = float(np.median(t_obs[trans_cols] ** 2) / 0.4549)
print(f"  lambda (trans, from t) = {lam_all:.4f}")

Gt = np.asfortranarray(Gr[:, trans_cols].astype(np.float32))
yr32 = yr.astype(np.float32)
maxT = np.empty(NPERM, dtype=np.float32)
rng = np.random.default_rng(SEED)
done = 0
while done < NPERM:
    B = min(BATCH, NPERM - done)
    U = rng.random((B, N_PANEL))
    rank = np.argsort(np.argsort(U, axis=1), axis=1)
    perm = np.argsort(rank[:, idx], axis=1)
    Yp = yr32[perm]
    Yp -= Yp.mean(axis=1, keepdims=True)
    Yp /= Yp.std(axis=1, keepdims=True)
    R = (Yp @ Gt) / n
    np.clip(R, -0.999999, 0.999999, out=R)
    T = np.abs(R * np.sqrt(df) / np.sqrt(1 - R ** 2))
    maxT[done:done + B] = T.max(axis=1)
    done += B
    if done % 20000 == 0:
        print(f"  perms {done}/{NPERM}")
np.save(f"{OUT}/{TARGET}_maxT.npy", maxT)

def frame(cols_idx, tag):
    D = Vk.iloc[cols_idx][["sv_locus", "chrom", "pos", "svtype", "svlen", "af"]].copy()
    D["t"] = t_obs[cols_idx]
    D["abs_t"] = np.abs(t_obs[cols_idx])
    D["p_param"] = p_par[cols_idx]
    D["phenotype"] = PHENO
    D["target_chrom"] = TARGET
    D["n"] = n
    D["relation"] = tag
    return D.sort_values("p_param")


CAND_P = 0.01
tr = frame(trans_cols, "trans")
tr[tr.p_param <= CAND_P].to_csv(f"{OUT}/{TARGET}_cand.tsv", sep="\t", index=False,
                                float_format="%.6g")
frame(np.where(cis_k)[0], "cis").to_csv(f"{OUT}/{TARGET}_cis.tsv", sep="\t", index=False,
                                        float_format="%.6g")
q = np.quantile(np.abs(t_obs[trans_cols]), [.5, .9, .99, .999, 1.0])
pd.DataFrame([dict(phenotype=PHENO, target_chrom=TARGET, n=n, df=df,
                   n_trans=len(trans_cols), n_cis=int(cis_k.sum()),
                   lambda_trans=lam_all, nperm=NPERM,
                   n_cand_p01=int((tr.p_param <= CAND_P).sum()),
                   exp_cand_p01=CAND_P * len(trans_cols),
                   maxT_median=float(np.median(maxT)), maxT_p95=float(np.quantile(maxT, .95)),
                   obs_max_abs_t=float(np.abs(t_obs[trans_cols]).max()),
                   absT_q50=q[0], absT_q90=q[1], absT_q99=q[2], absT_q999=q[3])
              ]).to_csv(f"{OUT}/{TARGET}_summary.tsv", sep="\t", index=False, float_format="%.6g")
print(f"  candidates p<={CAND_P}: {int((tr.p_param <= CAND_P).sum())} "
      f"(expected {CAND_P*len(trans_cols):.0f})   observed max |t| = {np.abs(t_obs[trans_cols]).max():.3f}"
      f"   cell maxT p95 = {np.quantile(maxT, .95):.3f}")
print("done")
