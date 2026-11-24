#!/usr/bin/env python3
import json
import struct
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
CDR = ROOT / "analyses/cdr_qtl"
GRM_PREFIX = ROOT / "subproject_methylation/results/greml/pericen_1m_all"
OUT = CDR / "results/haplogroup_chr8"
CHROM, PHENO, SV_POS = "chr8", "cdr_pos_offset", 44300221
LD_EXCLUDE = 0.8
N_PC = 10


def read_gcta_grm_bin(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    with open(f"{prefix}.grm.bin", "rb") as f:
        vals = struct.unpack("f" * (n * (n + 1) // 2), f.read(4 * n * (n + 1) // 2))
    K = np.zeros((n, n)); k = 0
    for i in range(n):
        for j in range(i + 1):
            K[i, j] = K[j, i] = vals[k]; k += 1
    return K, [l.strip().split()[1] for l in open(f"{prefix}.grm.id")]


def fit(y, cols):
    design = np.column_stack([np.ones(len(y))] + [c if c.ndim == 2 else c[:, None] for c in cols])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    df = max(1, len(y) - np.linalg.matrix_rank(design))
    var = np.sum(resid ** 2) / df
    se = np.sqrt(np.maximum(0, np.diag(np.linalg.pinv(design.T @ design) * var)))
    return beta, se, df, resid


sv_dir = CDR / "canonical_geno/perchrom" / CHROM
G = np.load(sv_dir / f"{CHROM}_sv_genotypes.npz", allow_pickle=True)["genotypes"]
sv_samples = [s.strip() for s in open(sv_dir / f"{CHROM}_sv_samples.txt") if s.strip()]
sv_var = pd.read_csv(sv_dir / f"{CHROM}_sv_variants.tsv", sep="\t")
if G.shape[0] == len(sv_var):
    G = G.T

ph = pd.read_csv(CDR / "output/per_hap_chrom_phenotypes.tsv", sep="\t")
pcc = ph[ph["chm13_chrom"] == CHROM]
pheno_ind = pcc.groupby("sample")[PHENO].mean().reset_index()
bc_map = pcc.groupby("sample")["basecaller"].first().to_dict()
K, grm_samples = read_gcta_grm_bin(GRM_PREFIX)
ev, evec = np.linalg.eigh(K)
PCs = evec[:, np.argsort(ev)[::-1][:10]]
pcm = pd.DataFrame(PCs, columns=[f"PC{i+1}" for i in range(10)]); pcm["sample"] = grm_samples

keep = sorted(set(bc_map) & set(pheno_ind["sample"]) & set(grm_samples) & set(sv_samples))
pmap = pheno_ind.set_index("sample")[PHENO].to_dict()
pcmap = pcm.set_index("sample")[[f"PC{i+1}" for i in range(10)]].to_dict("index")
sidx = {s: i for i, s in enumerate(sv_samples)}
y = np.array([pmap[s] for s in keep], float)
GRMPC = np.array([[pcmap[s][f"PC{i+1}"] for i in range(10)] for s in keep], float)
bc = np.array([1.0 if bc_map.get(s) == "R1041" else 0.0 for s in keep], float)
Gk = G[[sidx[s] for s in keep], :].astype(float)
m = ~np.isnan(y)
y, GRMPC, bc, Gk = y[m], GRMPC[m], bc[m], Gk[m]
keep = [k for i, k in enumerate(keep) if m[i]]
n = len(y)
y_rint = stats.norm.ppf((stats.rankdata(y) - 0.5) / n)
cov = np.column_stack([GRMPC, bc[:, None]]) if bc.std() > 1e-10 else GRMPC
print(f"n = {n}; chr8 predictors = {Gk.shape[1]}")

lead_i = sv_var.index[(sv_var["chrom"] == CHROM) & (sv_var["pos"] == SV_POS) &
                      (sv_var["cv_r2"] > 0.9)].tolist()[0]
lead = Gk[:, lead_i]
b, se, df, _ = fit(y_rint, [cov, lead])
t_lead = b[-1] / se[-1]
print(f"lead {sv_var.iloc[lead_i]['sv_locus']}: t = {t_lead:+.3f} "
      f"P = {2*stats.t.sf(abs(t_lead), df):.3e}   [gate: expect t = -6.17]")
assert abs(t_lead + 6.17) < 0.05, "lead does not reproduce; design mis-assembled"

sd = Gk.std(axis=0)
poly = sd > 0
r_lead = np.zeros(Gk.shape[1])
lc = lead - lead.mean()
for j in np.where(poly)[0]:
    r_lead[j] = abs(np.corrcoef(lc, Gk[:, j])[0, 1])
block = (r_lead >= LD_EXCLUDE) | ~poly
resid_cols = np.where(~block)[0]
print(f"lead LD block (|r| >= {LD_EXCLUDE}) or monomorphic: {block.sum()} predictors excluded")
print(f"residual chr8 centromeric SV panel: {len(resid_cols)} predictors")

R = Gk[:, resid_cols]
Rz = (R - R.mean(axis=0)) / R.std(axis=0)
U, S, Vt = np.linalg.svd(Rz - Rz.mean(axis=0), full_matrices=False)
varexp = (S ** 2) / np.sum(S ** 2)
svpc = U[:, :N_PC] * S[:N_PC]
print(f"centromeric SV PCs: var explained {np.round(varexp[:N_PC]*100,1).tolist()}")

print("\nQ1  dip position ~ residual chr8 centromeric SV PCs")
rows = []
for k in range(N_PC):
    b, se, df, _ = fit(y_rint, [cov, svpc[:, k]])
    t = b[-1] / se[-1]; p = 2 * stats.t.sf(abs(t), df)
    rows.append({"component": f"svPC{k+1}", "var_explained": float(varexp[k]),
                 "beta": float(b[-1]), "t": float(t), "p": float(p)})
    print(f"  svPC{k+1:<2d} varexp={varexp[k]*100:5.1f}%  t={t:+6.2f}  P={p:.3g}")
pmin = min(r["p"] for r in rows)
print(f"  Bonferroni over {N_PC} PCs: alpha = {0.05/N_PC:.4f}; min P = {pmin:.3g} "
      f"-> {'SIGNIFICANT' if pmin < 0.05/N_PC else 'not significant'}")

link = linkage(pdist(Rz, metric="correlation"), method="average")
clust_rows = []
for k in (2, 3, 4, 5):
    lab = fcluster(link, k, criterion="maxclust")
    D = np.column_stack([(lab == c).astype(float) for c in range(2, k + 1)])
    b0, se0, df0, r0 = fit(y_rint, [cov])
    b1, se1, df1, r1 = fit(y_rint, [cov, D])
    F = ((np.sum(r0**2) - np.sum(r1**2)) / D.shape[1]) / (np.sum(r1**2) / df1)
    p = stats.f.sf(F, D.shape[1], df1)
    sizes = np.bincount(lab)[1:].tolist()
    clust_rows.append({"k": k, "sizes": sizes, "F": float(F), "p": float(p)})
    print(f"  k={k} sizes={sizes}  F={F:.2f}  P={p:.3g}")

print("\nQ2  lead SV conditioned on the residual centromeric structure")
q2 = {}
for label, extra in (("svPC1-5", svpc[:, :5]), ("svPC1-10", svpc[:, :N_PC])):
    b, se, df, _ = fit(y_rint, [cov, extra, lead])
    t = b[-1] / se[-1]; p = 2 * stats.t.sf(abs(t), df)
    q2[label] = {"beta": float(b[-1]), "t": float(t), "p": float(p),
                 "retained_frac": float(abs(b[-1]) / abs(fit(y_rint, [cov, lead])[0][-1]))}
    print(f"  + {label:9s}  lead t={t:+6.2f} P={p:.3g}  "
          f"({q2[label]['retained_frac']*100:.0f}% of marginal beta)")

OUT.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(OUT / "svpc_associations.tsv", sep="\t", index=False)
pd.DataFrame(clust_rows).to_csv(OUT / "cluster_associations.tsv", sep="\t", index=False)
summary = {"n": n, "n_chr8_predictors": int(Gk.shape[1]),
           "n_excluded_lead_block": int(block.sum()),
           "n_residual_panel": int(len(resid_cols)),
           "ld_exclude_threshold": LD_EXCLUDE,
           "lead_marginal_t": float(t_lead),
           "svpc_min_p": float(pmin), "svpc_bonferroni_alpha": 0.05 / N_PC,
           "lead_conditioned": q2, "clusters": clust_rows}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print("\nwrote", OUT / "summary.json")
