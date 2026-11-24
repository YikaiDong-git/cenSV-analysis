#!/usr/bin/env python3
"""Head-to-head: the chr8 deletion versus its single best tag SNP, on dip position."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import struct
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
CDR = ROOT / "analyses/cdr_qtl"
GRM_PREFIX = ROOT / "subproject_methylation/results/greml/pericen_1m_all"
OUT = CDR / "results/snp_conditioning_chr8"
CHROM, PHENO, SV_POS = "chr8", "cdr_pos_offset", 44300221
PANEL = "midpoint1000"


def read_gcta_grm_bin(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    with open(f"{prefix}.grm.bin", "rb") as f:
        vals = struct.unpack("f" * (n * (n + 1) // 2), f.read(4 * n * (n + 1) // 2))
    K = np.zeros((n, n)); k = 0
    for i in range(n):
        for j in range(i + 1):
            K[i, j] = K[j, i] = vals[k]; k += 1
    return K, [l.strip().split()[1] for l in open(f"{prefix}.grm.id")]


def fit(y, X):
    design = np.column_stack([np.ones(len(y))] + [x if x.ndim == 2 else x[:, None] for x in X])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    df = max(1, len(y) - np.linalg.matrix_rank(design))
    var = np.sum(resid ** 2) / df
    cov = np.linalg.pinv(design.T @ design) * var
    return beta, np.sqrt(np.maximum(0, np.diag(cov))), df


sv_dir = CDR / "canonical_geno/perchrom" / CHROM
sv_geno = np.load(sv_dir / f"{CHROM}_sv_genotypes.npz", allow_pickle=True)["genotypes"]
sv_samples = [s.strip() for s in open(sv_dir / f"{CHROM}_sv_samples.txt") if s.strip()]
sv_var = pd.read_csv(sv_dir / f"{CHROM}_sv_variants.tsv", sep="\t")
if sv_geno.shape[0] == len(sv_var):
    sv_geno = sv_geno.T
print("SV genotype value set:", sorted(set(np.unique(sv_geno).tolist()))[:10],
      "dtype", sv_geno.dtype, "-> assembly genotypes, not imputed dosages")

ph = pd.read_csv(CDR / "output/per_hap_chrom_phenotypes.tsv", sep="\t")
pcc = ph[ph["chm13_chrom"] == CHROM]
pheno_ind = pcc.groupby("sample")[PHENO].mean().reset_index()
bc_map = pcc.groupby("sample")["basecaller"].first().to_dict()
K, grm_samples = read_gcta_grm_bin(GRM_PREFIX)
ev, evec = np.linalg.eigh(K)
PCs = evec[:, np.argsort(ev)[::-1][:10]]
pc_df = pd.DataFrame(PCs, columns=[f"PC{i+1}" for i in range(10)]); pc_df["sample"] = grm_samples

keep = sorted(set(bc_map) & set(pheno_ind["sample"]) & set(grm_samples) & set(sv_samples))
pmap = pheno_ind.set_index("sample")[PHENO].to_dict()
pcmap = pc_df.set_index("sample")[[f"PC{i+1}" for i in range(10)]].to_dict("index")
sidx = {s: i for i, s in enumerate(sv_samples)}
y = np.array([pmap[s] for s in keep], float)
PC_mat = np.array([[pcmap[s][f"PC{i+1}"] for i in range(10)] for s in keep], float)
bc = np.array([1.0 if bc_map.get(s) == "R1041" else 0.0 for s in keep], float)
svrow = sv_var.index[(sv_var["chrom"] == CHROM) & (sv_var["pos"] == SV_POS) &
                     (sv_var["cv_r2"] > 0.9)].tolist()[0]
sv = sv_geno[[sidx[s] for s in keep], svrow].astype(float)
m = ~np.isnan(y)
y, PC_mat, bc, sv = y[m], PC_mat[m], bc[m], sv[m]
keep = [k for i, k in enumerate(keep) if m[i]]
n = len(y)
y_rint = stats.norm.ppf((stats.rankdata(y) - 0.5) / n)
cov = np.column_stack([PC_mat, bc[:, None]]) if bc.std() > 1e-10 else PC_mat
print(f"n = {n}; SV genotype counts {dict(zip(*np.unique(sv, return_counts=True)))}")

G = np.load(ROOT / "snp_extraction" / PANEL / f"{CHROM}_snp_genotypes.npz",
            allow_pickle=True)["genotypes"].T[[sidx[s] for s in keep], :].astype(float)
var_tsv = ROOT / "snp_extraction" / PANEL / f"{CHROM}_snp_variants.tsv"
snpmeta = pd.read_csv(var_tsv, sep="\t") if var_tsv.exists() else None
keepcol = np.nanstd(G, axis=0) > 0
Gk = G[:, keepcol]
orig = np.where(keepcol)[0]

svc = sv - sv.mean(); r_raw = np.array([
    abs(np.corrcoef(svc, Gk[:, j])[0, 1]) for j in range(Gk.shape[1])])
best = int(np.nanargmax(r_raw))
print(f"\nbest tag SNP: panel column {orig[best]}  RAW |r| with the deletion = {r_raw[best]:.4f}")
print(f"  SNPs at raw |r| >= 0.9: {int(np.sum(r_raw >= 0.9))} / {Gk.shape[1]}")
print(f"  SNPs at raw |r| >= 0.8: {int(np.sum(r_raw >= 0.8))} / {Gk.shape[1]}")
if snpmeta is not None and len(snpmeta) > orig[best]:
    print("  metadata:", snpmeta.iloc[orig[best]].to_dict())

tag = Gk[:, best]
b_sv, se_sv, df1 = fit(y_rint, [cov, sv])
b_tag, se_tag, df2 = fit(y_rint, [cov, tag])
b_both, se_both, df3 = fit(y_rint, [cov, sv, tag])
t = lambda b, s: b / s
p = lambda b, s, d: 2 * stats.t.sf(abs(b / s), d)

print("\nHEAD-TO-HEAD (RINT dip position ~ covariates + ...)")
print(f"  SV alone       beta={b_sv[-1]:+.4f} se={se_sv[-1]:.4f} "
      f"t={t(b_sv[-1],se_sv[-1]):+.3f} p={p(b_sv[-1],se_sv[-1],df1):.3e}")
print(f"  tag SNP alone  beta={b_tag[-1]:+.4f} se={se_tag[-1]:.4f} "
      f"t={t(b_tag[-1],se_tag[-1]):+.3f} p={p(b_tag[-1],se_tag[-1],df2):.3e}")
print(f"  both: SV       beta={b_both[-2]:+.4f} se={se_both[-2]:.4f} "
      f"t={t(b_both[-2],se_both[-2]):+.3f} p={p(b_both[-2],se_both[-2],df3):.3e}")
print(f"  both: tag SNP  beta={b_both[-1]:+.4f} se={se_both[-1]:.4f} "
      f"t={t(b_both[-1],se_both[-1]):+.3f} p={p(b_both[-1],se_both[-1],df3):.3e}")
print(f"\n  SE inflation on the SV when the tag is added: "
      f"{se_both[-2]/se_sv[-1]:.2f}x")

summary = {
    "n": n, "panel": PANEL,
    "raw_max_abs_r_sv_tag": float(r_raw[best]),
    "n_snps_raw_r_ge_0.9": int(np.sum(r_raw >= 0.9)),
    "n_snps_raw_r_ge_0.8": int(np.sum(r_raw >= 0.8)),
    "n_snps_panel_polymorphic": int(Gk.shape[1]),
    "sv_alone": {"beta": float(b_sv[-1]), "p": float(p(b_sv[-1], se_sv[-1], df1))},
    "tag_alone": {"beta": float(b_tag[-1]), "p": float(p(b_tag[-1], se_tag[-1], df2))},
    "joint_sv": {"beta": float(b_both[-2]), "p": float(p(b_both[-2], se_both[-2], df3))},
    "joint_tag": {"beta": float(b_both[-1]), "p": float(p(b_both[-1], se_both[-1], df3))},
    "se_inflation_sv": float(se_both[-2] / se_sv[-1]),
}
(OUT / "tag_snp_head_to_head.json").write_text(json.dumps(summary, indent=2) + "\n")
print("\nwrote", OUT / "tag_snp_head_to_head.json")
