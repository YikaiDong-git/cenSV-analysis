#!/usr/bin/env python3
"""ADVERSARIAL CHECK on the chr8 haplogroup result."""
import json
import struct
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
CDR = ROOT / "analyses/cdr_qtl"
GRM_PREFIX = ROOT / "subproject_methylation/results/greml/pericen_1m_all"
OUT = CDR / "results/haplogroup_chr8"
CHROM, PHENO, SV_POS = "chr8", "cdr_pos_offset", 44300221
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
    d = np.column_stack([np.ones(len(y))] + [c if c.ndim == 2 else c[:, None] for c in cols])
    b, *_ = np.linalg.lstsq(d, y, rcond=None)
    r = y - d @ b
    df = max(1, len(y) - np.linalg.matrix_rank(d))
    v = np.sum(r ** 2) / df
    se = np.sqrt(np.maximum(0, np.diag(np.linalg.pinv(d.T @ d) * v)))
    return b, se, df, r


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
lead_i = sv_var.index[(sv_var["chrom"] == CHROM) & (sv_var["pos"] == SV_POS) &
                      (sv_var["cv_r2"] > 0.9)].tolist()[0]
lead = Gk[:, lead_i]

sd = Gk.std(axis=0); poly = sd > 0
r_lead = np.zeros(Gk.shape[1])
lc = lead - lead.mean()
for j in np.where(poly)[0]:
    r_lead[j] = abs(np.corrcoef(lc, Gk[:, j])[0, 1])

print("distribution of |r| with the lead across the 242 chr8 predictors:")
for lo, hi in ((0.0, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 0.8), (0.8, 1.01)):
    print(f"  {lo:.1f}-{hi:.1f}: {int(((r_lead >= lo) & (r_lead < hi) & poly).sum())}")

results = {}
print("\nD2  re-run at stricter exclusion thresholds")
print(f"{'thresh':>7} {'n_kept':>7} {'svPC1_varexp':>13} {'corr(svPC1,lead)':>17} "
      f"{'R2(lead~PCs)':>13} {'svPC1 t':>9} {'svPC1 P':>10} {'lead|PCs t':>11}")
for thr in (0.8, 0.5, 0.3, 0.2):
    block = (r_lead >= thr) | ~poly
    cols = np.where(~block)[0]
    if len(cols) < 12:
        print(f"{thr:>7.1f} {len(cols):>7} -- too few predictors, skipped")
        continue
    R = Gk[:, cols]
    Rz = (R - R.mean(axis=0)) / R.std(axis=0)
    U, S, _ = np.linalg.svd(Rz - Rz.mean(axis=0), full_matrices=False)
    ve = (S ** 2) / np.sum(S ** 2)
    svpc = U[:, :N_PC] * S[:N_PC]
    c1 = abs(np.corrcoef(svpc[:, 0], lead)[0, 1])
    d = np.column_stack([np.ones(n), svpc])
    bb, *_ = np.linalg.lstsq(d, lead, rcond=None)
    r2 = 1 - np.sum((lead - d @ bb) ** 2) / np.sum((lead - lead.mean()) ** 2)
    b, se, df, _ = fit(y_rint, [cov, svpc[:, 0]])
    t1 = b[-1] / se[-1]; p1 = 2 * stats.t.sf(abs(t1), df)
    b2, se2, df2, _ = fit(y_rint, [cov, svpc[:, :N_PC], lead])
    tl = b2[-1] / se2[-1]
    print(f"{thr:>7.1f} {len(cols):>7} {ve[0]*100:>12.1f}% {c1:>17.3f} {r2:>13.3f} "
          f"{t1:>9.2f} {p1:>10.3g} {tl:>11.2f}")
    results[str(thr)] = {"n_kept": int(len(cols)), "svPC1_varexp": float(ve[0]),
                         "corr_svPC1_lead": float(c1), "r2_lead_from_PCs": float(r2),
                         "svPC1_t": float(t1), "svPC1_p": float(p1),
                         "lead_given_PCs_t": float(tl)}

Gs = np.load(ROOT / "snp_extraction/midpoint1000" / f"{CHROM}_snp_genotypes.npz",
             allow_pickle=True)["genotypes"].T[[sidx[s] for s in keep], :].astype(float)
meta = pd.read_csv(ROOT / "snp_extraction/midpoint1000" / f"{CHROM}_snp_variants.tsv", sep="\t")
tag_j = int(meta.index[meta["pos"] == 46190487][0])
tag = Gs[:, tag_j]
block = (r_lead >= 0.8) | ~poly
R = Gk[:, np.where(~block)[0]]
Rz = (R - R.mean(axis=0)) / R.std(axis=0)
U, S, _ = np.linalg.svd(Rz - Rz.mean(axis=0), full_matrices=False)
svpc1 = (U[:, :1] * S[:1]).ravel()
print(f"\nD3  corr(svPC1, tag SNP 8:46,190,487) = {abs(np.corrcoef(svpc1, tag)[0,1]):.3f}")
print(f"    corr(svPC1, lead SV)               = {abs(np.corrcoef(svpc1, lead)[0,1]):.3f}")
b, se, df, _ = fit(y_rint, [cov, svpc1, tag])
print(f"    dip ~ cov + svPC1 + tagSNP: svPC1 t={b[-2]/se[-2]:+.2f} "
      f"P={2*stats.t.sf(abs(b[-2]/se[-2]), df):.3g} | tag t={b[-1]/se[-1]:+.2f} "
      f"P={2*stats.t.sf(abs(b[-1]/se[-1]), df):.3g}")

results["D3"] = {"corr_svPC1_tagSNP": float(abs(np.corrcoef(svpc1, tag)[0, 1])),
                 "corr_svPC1_lead": float(abs(np.corrcoef(svpc1, lead)[0, 1]))}
(OUT / "leakage_check.json").write_text(json.dumps(results, indent=2) + "\n")
print("\nwrote", OUT / "leakage_check.json")
