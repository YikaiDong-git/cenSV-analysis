#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Per-allele effect size of the chr8 113-bp deletion on CDR position, in kb."""
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
U = ROOT / "analyses/cdr_definition"
WU = ROOT / "analyses/cdr_qtl_full_panel"
CANON = ROOT / "analyses/cdr_qtl"
sys.path.insert(0, str(CANON / "code/scan"))

CHROM = "chr8"
PHENO = "cdr_pos_offset"
LEAD_POS = 44300221
PUBLISHED_T = -6.17302
GRM = str(ROOT / "subproject_methylation/results/greml/pericen_1m_all")

from b6_cdr_pos_qtl_perchrom import read_gcta_grm_bin  # noqa: E402

SV_DIR = WU / "geno/perchrom" / CHROM
_npz = np.load(SV_DIR / ("%s_sv_genotypes.npz" % CHROM))
_key = "arr_0" if "arr_0" in _npz else list(_npz.keys())[0]
sv_dosages = _npz[_key]
sv_samples = open(SV_DIR / ("%s_sv_samples.txt" % CHROM)).read().split()
sv_variants = pd.read_csv(SV_DIR / ("%s_sv_variants.tsv" % CHROM), sep="\t")
if sv_dosages.shape[0] != len(sv_samples):
    sv_dosages = sv_dosages.T

PF = U / "output/compat/compat_dip_splice/per_hap_chrom_phenotypes.tsv"
pc = pd.read_csv(PF, sep="\t")
pc_chrom = pc[pc["chm13_chrom"] == CHROM]
pheno_ind = pc_chrom.groupby("sample")[PHENO].mean().reset_index()
bc_map = pc_chrom.groupby("sample")["basecaller"].first().to_dict()

K, grm_samples = read_gcta_grm_bin(GRM)
eigvals, eigvecs = np.linalg.eigh(K)
PCs = eigvecs[:, np.argsort(eigvals)[::-1][:10]]
pc_df = pd.DataFrame(PCs, columns=["PC%d" % (i + 1) for i in range(10)])
pc_df["sample"] = grm_samples

keep = sorted(set(bc_map) & set(pheno_ind["sample"]) & set(grm_samples) & set(sv_samples))
pheno_map = pheno_ind.set_index("sample")[PHENO].to_dict()
pc_map = pc_df.set_index("sample")[["PC%d" % (i + 1) for i in range(10)]].to_dict("index")
sv_idx = {s: i for i, s in enumerate(sv_samples)}

y = np.array([pheno_map[s] for s in keep], float)
PC_mat = np.array([[pc_map[s]["PC%d" % (i + 1)] for i in range(10)] for s in keep], float)
bc_vec = np.array([1.0 if bc_map.get(s) == "R1041" else 0.0 for s in keep], float)
sv_sub = sv_dosages[[sv_idx[s] for s in keep], :].astype(float)

m = ~np.isnan(y)
y, PC_mat, bc_vec, sv_sub = y[m], PC_mat[m], bc_vec[m], sv_sub[m]
keep = [k for i, k in enumerate(keep) if m[i]]
n = len(y)

lead_j = int(np.where(sv_variants["pos"].values == LEAD_POS)[0][0])
g = sv_sub[:, lead_j]

X = (np.column_stack([np.ones(n), PC_mat, bc_vec])
     if len(set(bc_vec)) > 1 else np.column_stack([np.ones(n), PC_mat]))
df_resid = n - X.shape[1] - 1


def resid(v):
    return v - X @ np.linalg.lstsq(X, v, rcond=None)[0]


def t_of(yv, gv):
    yr, gr = resid(yv), resid(gv)
    yz = (yr - yr.mean()) / yr.std()
    gz = (gr - gr.mean()) / gr.std()
    r = float(np.dot(yz, gz) / n)
    return r * np.sqrt(df_resid / max(1.0 - r * r, 1e-300)), r, yr, gr


def rint(a):
    from scipy.stats import norm
    return norm.ppf((pd.Series(a).rank(method="average").values - 0.5) / len(a))


t_pub, r_pub, _, _ = t_of(rint(y), g)
print("GATE  published t = %.5f   refit t = %.5f   delta = %.3e"
      % (PUBLISHED_T, t_pub, abs(t_pub - PUBLISHED_T)))
if abs(t_pub - PUBLISHED_T) > 1e-4:
    sys.exit("[FATAL] refit does not reproduce the published model; raw-scale numbers withheld")
print("GATE PASSED  n=%d  df_resid=%d  covariates=%d\n" % (n, df_resid, X.shape[1]))

t_raw, r_raw, y_res, g_res = t_of(y, g)
beta_splice = float(np.dot(g_res, y_res) / np.dot(g_res, g_res))
se_splice = beta_splice / t_raw
print("RAW SPLICE SCALE (no RINT)")
print("  beta = %+.5f SPLICE units per allele   t = %.4f   r = %+.4f" % (beta_splice, t_raw, r_raw))
print("  95%% CI %+.5f .. %+.5f" % (beta_splice - 1.96 * se_splice, beta_splice + 1.96 * se_splice))

_LEN = "active_len" if "active_len" in pc_chrom.columns else "active_len_sum"
if _LEN not in pc_chrom.columns:
    sys.exit("[FATAL] no active-array length column; kb not derivable")
al = pc_chrom.groupby("sample")[_LEN].mean()
print("  length column used: %s" % _LEN)
L = np.array([al.get(s, np.nan) for s in keep], float)
ok = ~np.isnan(L)
half = L[ok] / 2.0
kb = beta_splice * half / 1000.0
print("\nKB CONVERSION  (1 SPLICE unit = half the active array)")
print("  chr8 active array  n=%d  min %.0f kb  median %.0f kb  max %.0f kb"
      % (ok.sum(), L[ok].min() / 1e3, np.median(L[ok]) / 1e3, L[ok].max() / 1e3))
print("  per-allele shift   median %+.1f kb   IQR %+.1f .. %+.1f   range %+.1f .. %+.1f"
      % (np.median(kb), np.percentile(kb, 25), np.percentile(kb, 75), kb.min(), kb.max()))
print("  mean per-allele shift %+.1f kb" % kb.mean())

print("\nMODEL-FREE CHECK  mean raw pos_splice by rounded genotype")
for gv in sorted(set(np.round(g))):
    s = np.round(g) == gv
    if s.sum():
        print("  dosage %d  n=%3d  mean pos_splice %+.4f  mean kb-from-centre %+.1f"
              % (gv, s.sum(), y[s].mean(), y[s].mean() * np.median(half) / 1000.0))
