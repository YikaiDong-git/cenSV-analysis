#!/usr/bin/env python3
"""SNP-conditioning arm for the chr8 CDR dip-position lead."""
import json
import struct
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
CDR = ROOT / "analyses/cdr_qtl"
GRM_PREFIX = ROOT / "subproject_methylation/results/greml/pericen_1m_all"
OUT = CDR / "results/snp_conditioning_chr8"

CHROM = "chr8"
PHENO = "cdr_pos_offset"
SV_POS = 44300221
PUBLISHED_T = -6.17

MAX_SUSIE_SNPS = 300
MAX_CONDITION_SNPS = 20
LD_THRESHOLD = 0.8
RIDGE_MAX_SNPS = 300
RIDGE_FOLDS = 5
SEED = 104729


def residualize(matrix, covariates):
    design = np.column_stack([np.ones(len(covariates)), covariates])
    q, _ = np.linalg.qr(design, mode="reduced")
    return matrix - q @ (q.T @ matrix)


def fit_sv(y, sv, extras=None):
    columns = [np.ones(len(y))]
    if extras is not None and np.size(extras):
        columns.append(extras if extras.ndim == 2 else extras[:, None])
    columns.append(sv[:, None])
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ beta
    df = max(1, len(y) - np.linalg.matrix_rank(design))
    variance = np.sum(residual ** 2) / df
    se = np.sqrt(max(0.0, np.linalg.pinv(design.T @ design)[-1, -1] * variance))
    t = beta[-1] / se if se > 0 else np.nan
    return {"beta": float(beta[-1]), "se": float(se), "t": float(t),
            "pvalue": float(2 * stats.t.sf(abs(t), df)) if np.isfinite(t) else np.nan,
            "df": int(df)}


def safe_abs_correlation(x, y):
    xc = x - x.mean(axis=0)
    yc = y - y.mean()
    num = np.abs(xc.T @ yc)
    den = np.sqrt(np.sum(xc ** 2, axis=0) * np.sum(yc ** 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, num / den, 0.0)


def select_snps(y, sv, genotypes, maximum):
    y_score = np.abs(genotypes.T @ y) / np.sqrt(np.sum(genotypes ** 2, axis=0) * np.sum(y ** 2))
    sv_score = np.abs(genotypes.T @ sv) / np.sqrt(np.sum(genotypes ** 2, axis=0) * np.sum(sv ** 2))
    half = max(1, maximum // 2)
    selected = np.unique(np.r_[np.argsort(y_score)[-half:], np.argsort(sv_score)[-half:]])
    return selected[np.argsort(np.maximum(y_score[selected], sv_score[selected]))[-maximum:]]


def prune(genotypes, order, maximum, threshold=LD_THRESHOLD):
    selected = []
    sd = genotypes.std(axis=0)
    sd[sd == 0] = 1.0
    standardized = (genotypes - genotypes.mean(axis=0)) / sd
    for index in order[::-1]:
        if not selected or np.all(np.abs(
                standardized[:, selected].T @ standardized[:, index] / len(genotypes)) < threshold):
            selected.append(int(index))
        if len(selected) == maximum:
            break
    return selected


def ridge_residual_test(y, sv, genotypes, covariates, maximum, folds):
    if maximum <= 0 or genotypes.shape[1] == 0:
        return {k: np.nan for k in ("ridge_oof_r2", "ridge_y_oof_r2", "ridge_resid_beta",
                                    "ridge_resid_se", "ridge_resid_p")}
    sv_pred = np.zeros(len(sv)); y_pred = np.zeros(len(y))
    sv_obs = np.zeros(len(sv)); y_obs = np.zeros(len(y))
    splitter = KFold(n_splits=min(folds, len(sv)), shuffle=True, random_state=SEED)
    alphas = np.logspace(-3, 4, 16)
    for train, test in splitter.split(genotypes):
        train_design = np.column_stack([np.ones(len(train)), covariates[train]])
        test_design = np.column_stack([np.ones(len(test)), covariates[test]])
        gcoef, *_ = np.linalg.lstsq(train_design, genotypes[train], rcond=None)
        ocoef, *_ = np.linalg.lstsq(train_design, np.column_stack([sv[train], y[train]]), rcond=None)
        x_tr = genotypes[train] - train_design @ gcoef
        x_te = genotypes[test] - test_design @ gcoef
        o_tr = np.column_stack([sv[train], y[train]]) - train_design @ ocoef
        o_te = np.column_stack([sv[test], y[test]]) - test_design @ ocoef
        scores = np.maximum(safe_abs_correlation(x_tr, o_tr[:, 0]),
                            safe_abs_correlation(x_tr, o_tr[:, 1]))
        sel = np.argsort(scores)[-min(maximum, genotypes.shape[1]):]
        scaler = StandardScaler()
        xtr = scaler.fit_transform(x_tr[:, sel])
        xte = scaler.transform(x_te[:, sel])
        sv_pred[test] = RidgeCV(alphas=alphas).fit(xtr, o_tr[:, 0]).predict(xte)
        y_pred[test] = RidgeCV(alphas=alphas).fit(xtr, o_tr[:, 1]).predict(xte)
        sv_obs[test] = o_te[:, 0]; y_obs[test] = o_te[:, 1]
    sv_r = sv_obs - sv_pred; y_r = y_obs - y_pred
    den = np.sum(sv_r ** 2)
    beta = float(np.sum(sv_r * y_r) / den) if den > 0 else np.nan
    score = sv_r * (y_r - beta * sv_r)
    moment = np.mean(sv_r ** 2)
    se = float(np.sqrt(np.mean(score ** 2) / len(sv)) / moment) if moment > 0 else np.nan
    z = beta / se if np.isfinite(se) and se > 0 else np.nan
    svd = np.sum((sv_obs - sv_obs.mean()) ** 2); yd = np.sum((y_obs - y_obs.mean()) ** 2)
    return {"ridge_oof_r2": float(1 - np.sum(sv_r ** 2) / svd) if svd > 0 else np.nan,
            "ridge_y_oof_r2": float(1 - np.sum(y_r ** 2) / yd) if yd > 0 else np.nan,
            "ridge_resid_beta": beta, "ridge_resid_se": se,
            "ridge_resid_p": float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan}


def read_gcta_grm_bin(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    with open(f"{prefix}.grm.bin", "rb") as f:
        vals = struct.unpack("f" * (n * (n + 1) // 2), f.read(4 * n * (n + 1) // 2))
    K = np.zeros((n, n))
    k = 0
    for i in range(n):
        for j in range(i + 1):
            K[i, j] = K[j, i] = vals[k]; k += 1
    ids = [line.strip().split() for line in open(f"{prefix}.grm.id")]
    return K, [i[1] for i in ids]


print("=" * 74)
print("BUILDING THE PUBLISHED DESIGN")
sv_dir = CDR / "canonical_geno/perchrom" / CHROM
sv_geno = np.load(sv_dir / f"{CHROM}_sv_genotypes.npz", allow_pickle=True)["genotypes"]
sv_samples = [s.strip() for s in open(sv_dir / f"{CHROM}_sv_samples.txt") if s.strip()]
sv_var = pd.read_csv(sv_dir / f"{CHROM}_sv_variants.tsv", sep="\t")
if sv_geno.shape[0] == len(sv_var):
    sv_geno = sv_geno.T
print(f"  SV matrix {sv_geno.shape}  samples={len(sv_samples)}  variants={len(sv_var)}")

ph = pd.read_csv(CDR / "output/per_hap_chrom_phenotypes.tsv", sep="\t")
pc_chrom = ph[ph["chm13_chrom"] == CHROM]
pheno_ind = pc_chrom.groupby("sample")[PHENO].mean().reset_index()
bc_map = pc_chrom.groupby("sample")["basecaller"].first().to_dict()

K, grm_samples = read_gcta_grm_bin(GRM_PREFIX)
eigvals, eigvecs = np.linalg.eigh(K)
PCs = eigvecs[:, np.argsort(eigvals)[::-1][:10]]
pc_df = pd.DataFrame(PCs, columns=[f"PC{i+1}" for i in range(10)]); pc_df["sample"] = grm_samples

keep = sorted(set(bc_map) & set(pheno_ind["sample"]) & set(grm_samples) & set(sv_samples))
pheno_map = pheno_ind.set_index("sample")[PHENO].to_dict()
pc_map = pc_df.set_index("sample")[[f"PC{i+1}" for i in range(10)]].to_dict("index")
sv_idx = {s: i for i, s in enumerate(sv_samples)}

y = np.array([pheno_map[s] for s in keep], float)
PC_mat = np.array([[pc_map[s][f"PC{i+1}"] for i in range(10)] for s in keep], float)
bc = np.array([1.0 if bc_map.get(s) == "R1041" else 0.0 for s in keep], float)
sv_all = sv_geno[[sv_idx[s] for s in keep], :].astype(float)

m = ~np.isnan(y)
y, PC_mat, bc, sv_all = y[m], PC_mat[m], bc[m], sv_all[m]
keep = [k for i, k in enumerate(keep) if m[i]]
n = len(y)
print(f"  n after NaN drop = {n}")

rint = lambda a: stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))
y_rint = rint(y)
cov = np.column_stack([PC_mat, bc[:, None]]) if bc.std() > 1e-10 else PC_mat
X_pc = np.column_stack([np.ones(n), cov])
print(f"  covariates (incl. intercept) = {X_pc.shape[1]}")

print("=" * 74)
print(f"GATE -- reproduce the published marginal fit (t = {PUBLISHED_T})")
cands = sv_var.index[(sv_var["chrom"] == CHROM) & (sv_var["pos"] == SV_POS)].tolist()
print(f"  predictor rows at {CHROM}:{SV_POS:,} -> {len(cands)}")
y_resid = y_rint - X_pc @ np.linalg.lstsq(X_pc, y_rint, rcond=None)[0]
lead = None
for i in cands:
    g = sv_all[:, i]
    gr = g - X_pc @ np.linalg.lstsq(X_pc, g, rcond=None)[0]
    r = float(np.corrcoef(y_resid, gr)[0, 1])
    df = n - X_pc.shape[1] - 1
    t = r * np.sqrt(df / (1 - r ** 2))
    row = sv_var.iloc[i]
    flag = "  <== LEAD" if abs(t - PUBLISHED_T) < 0.05 else ""
    print(f"    {row['sv_locus']:44s} cv_r2={row['cv_r2']:.4f} t={t:+.3f}{flag}")
    if abs(t - PUBLISHED_T) < 0.05:
        lead = (i, row, t)
if lead is None:
    sys.exit(f"[FATAL] no predictor at {CHROM}:{SV_POS} reproduces t = {PUBLISHED_T}. "
             "The design is mis-assembled; refusing to run conditioning on it.")
lead_i, lead_row, lead_t = lead
print(f"  GATE PASSED: {lead_row['sv_locus']} reproduces t = {lead_t:+.4f}")

sv = sv_all[:, lead_i]
marginal = fit_sv(y_rint, sv, cov)
print(f"  marginal (covariate-adjusted): beta={marginal['beta']:+.4f} "
      f"p={marginal['pvalue']:.3e} df={marginal['df']}")

print("=" * 74)
print("CONDITIONING")
OUT.mkdir(parents=True, exist_ok=True)
panels = sorted(p.name for p in (ROOT / "snp_extraction").iterdir()
                if (p / f"{CHROM}_snp_genotypes.npz").exists())
rows = []
for panel in panels:
    g = np.load(ROOT / "snp_extraction" / panel / f"{CHROM}_snp_genotypes.npz",
                allow_pickle=True)["genotypes"]
    snp_samples = sv_samples if g.shape[1] == len(sv_samples) else None
    if snp_samples is None:
        print(f"  {panel}: sample axis {g.shape} unrecognised -- SKIPPED"); continue
    G = g.T[[sv_idx[s] for s in keep], :].astype(float)
    G = G[:, np.nanstd(G, axis=0) > 0]
    if np.isnan(G).any():
        G = np.where(np.isnan(G), np.nanmean(G, axis=0), G)
    n_snps = G.shape[1]
    if n_snps == 0:
        print(f"  {panel}: no polymorphic SNPs -- SKIPPED"); continue

    ridge = ridge_residual_test(y_rint, sv, G, cov, RIDGE_MAX_SNPS, RIDGE_FOLDS)
    Gr = residualize(G, cov)
    yv = residualize(y_rint[:, None], cov).ravel()
    svv = residualize(sv[:, None], cov).ravel()
    sel = select_snps(yv, svv, Gr, MAX_SUSIE_SNPS)
    Gs = Gr[:, sel]
    yscore = np.abs(Gs.T @ yv) / np.sqrt(np.sum(Gs ** 2, axis=0) * np.sum(yv ** 2))
    cond_idx = prune(Gs, np.argsort(yscore), MAX_CONDITION_SNPS)
    post = fit_sv(yv, svv, Gs[:, cond_idx])
    pre = fit_sv(yv, svv)
    att = abs(post["beta"] / pre["beta"]) if pre["beta"] else np.nan

    pd.DataFrame(np.column_stack([yv, svv, Gs]),
                 columns=["phenotype", "SV"] + [f"snp{j}" for j in range(Gs.shape[1])]
                 ).to_csv(OUT / f"susie_input_{panel}.tsv", sep="\t", index=False)

    rows.append({"panel": panel, "n": n, "n_snps_panel": n_snps,
                 "n_susie_snps": len(sel), "n_condition_snps": len(cond_idx),
                 "beta_pre": pre["beta"], "p_pre": pre["pvalue"],
                 "beta_post": post["beta"], "p_post": post["pvalue"],
                 "attenuation_ratio": att, **ridge})
    print(f"  {panel:15s} snps={n_snps:5d} cond={len(cond_idx):2d}  "
          f"p_pre={pre['pvalue']:.2e} -> p_post={post['pvalue']:.3g} "
          f"(retained {att*100:.1f}%)  ridge_p={ridge['ridge_resid_p']:.3g} "
          f"beta_ridge={ridge['ridge_resid_beta']:+.4f}")

res = pd.DataFrame(rows)
res.to_csv(OUT / "conditioning_by_panel.tsv", sep="\t", index=False)
meta = {"chrom": CHROM, "phenotype": PHENO, "sv_locus": str(lead_row["sv_locus"]),
        "sv_pos": int(SV_POS), "svlen": int(lead_row["svlen"]),
        "cv_r2": float(lead_row["cv_r2"]), "af": float(lead_row["af"]), "n": n,
        "published_t": PUBLISHED_T, "reproduced_t": float(lead_t),
        "marginal_beta": marginal["beta"], "marginal_p": marginal["pvalue"],
        "n_panels": len(rows)}
(OUT / "summary.json").write_text(json.dumps(meta, indent=2) + "\n")
print("=" * 74)
print(json.dumps(meta, indent=2))
print(f"\nwrote {OUT}/conditioning_by_panel.tsv  (+ per-panel SuSiE input matrices)")
