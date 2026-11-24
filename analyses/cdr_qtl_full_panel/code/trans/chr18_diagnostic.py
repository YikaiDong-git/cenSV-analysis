#!/usr/bin/env python
"""Attack the surviving trans association: chr8 INS-1868 -> chr18 CDR position, n = 54."""
import sys
import numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SUB = f"{P1}/subproject_methylation"
WU = f"{P1}/analyses/cdr_qtl_full_panel"
CANON = f"{P1}/analyses/cdr_qtl"
SV = "chr8_sv1_INS_non_best8_midpoint1000"
TARGET, PHENO = "chr18", "cdr_pos_offset_deepest"


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
BC = PC.groupby("sample")["basecaller"].first().to_dict()
Gm, gids = read_grm(f"{SUB}/results/greml/pericen_1m_all")
w, v = np.linalg.eigh(Gm)
PCMAP = {s: v[:, np.argsort(w)[::-1][:10]][i] for i, s in enumerate(gids)}

V = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t", engine="python")
Gz = np.load(f"{WU}/geno/full_genotypes.npz")
Gz = Gz[list(Gz.keys())[0]]
GS = open(f"{WU}/geno/full_samples.txt").read().split()
if Gz.shape[0] != len(GS):
    Gz = Gz.T
col = {s: i for i, s in enumerate(V.sv_locus)}[SV]
dos = pd.Series(Gz[:, col], index=GS).dropna()

y = PD_[PD_.chm13_chrom == TARGET].groupby("sample")[PHENO].mean().dropna()
have = set(y.index)
panel = [s for s in dos.index if s in PCMAP]
print(f"panel with genotype+PCs: {len(panel)}   with a {TARGET} {PHENO} call: "
      f"{len([s for s in panel if s in have])}")

avail = np.array([1.0 if s in have else 0.0 for s in panel])
d = dos.loc[panel].values
X = np.column_stack([np.ones(len(panel)), np.array([PCMAP[s] for s in panel])])
ar = avail - X @ np.linalg.lstsq(X, avail, rcond=None)[0]
dr = d - X @ np.linalg.lstsq(X, d, rcond=None)[0]
r = float((ar / ar.std()) @ (dr / dr.std()) / len(panel))
df = len(panel) - X.shape[1] - 1
t = r * np.sqrt(df) / np.sqrt(max(1 - r ** 2, 1e-30))
print(f"\n[1] COLLIDER CHECK -- does dosage predict having a {TARGET} call?")
print(f"    t = {t:+.3f}   p = {2*stats.t.sf(abs(t), df):.4f}   "
      f"mean dosage with call {d[avail==1].mean():.3f} vs without {d[avail==0].mean():.3f}")
print("    (a significant result here would mean the n=54 subset is selected on genotype)")

PED = f"{P1}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
ped = pd.read_csv(PED, sep=r"\s+")
sc = [c for c in ped.columns if c.lower() in ("sampleid", "individualid", "sample")][0]
spc = [c for c in ped.columns if "super" in c.lower()][0]
POP = ped.set_index(sc)[spc].to_dict()
a = pd.Series([POP.get(x, "NA") for x in panel if x in have]).value_counts()
b = pd.Series([POP.get(x, "NA") for x in panel if x not in have]).value_counts()
comp = pd.DataFrame({"with_call": a, "without_call": b}).fillna(0).astype(int)
print("\n[2] SUPERPOPULATION composition of the tested cell vs the rest")
print(comp.to_string())
print(f"    chi2 p = {stats.chi2_contingency(comp.values + 0.5).pvalue:.4f}")
carr = [x for x in panel if x in have and dos.loc[x] > 0]
noncarr = [x for x in panel if x in have and dos.loc[x] == 0]
print(f"    CARRIERS (n={len(carr)}) by superpop: "
      f"{dict(pd.Series([POP.get(x, 'NA') for x in carr]).value_counts())}")
print(f"    non-carriers (n={len(noncarr)}) by superpop: "
      f"{dict(pd.Series([POP.get(x, 'NA') for x in noncarr]).value_counts())}")
bcw = pd.Series([str(BC.get(s, "?")) for s in panel if s in have]).value_counts()
bco = pd.Series([str(BC.get(s, "?")) for s in panel if s not in have]).value_counts()
print(f"    basecaller with call {dict(bcw)}   without {dict(bco)}")

use = [s for s in panel if s in have]
x = dos.loc[use].values
yy = y.loc[use].values
print(f"\n[3] CARRIERS in the tested cell: n={len(use)}  "
      f"dosage counts {dict(pd.Series(x).value_counts().sort_index())}")


def fit(yv, xv, ss):
    Xc = np.column_stack([np.ones(len(ss)), np.array([PCMAP[s] for s in ss])])
    bc = np.array([1.0 if str(BC.get(s, "")) == "R1041" else 0.0 for s in ss])
    if bc.std() > 1e-10:
        Xc = np.column_stack([Xc, bc[:, None]])
    dfl = len(ss) - Xc.shape[1] - 1
    yr = yv - Xc @ np.linalg.lstsq(Xc, yv, rcond=None)[0]
    xr = xv - Xc @ np.linalg.lstsq(Xc, xv, rcond=None)[0]
    yr = (yr - yr.mean()) / yr.std(); xr = (xr - xr.mean()) / xr.std()
    rr = float(xr @ yr / len(ss))
    tt = rr * np.sqrt(dfl) / np.sqrt(max(1 - rr ** 2, 1e-30))
    return tt, 2 * stats.t.sf(abs(tt), dfl)


t0, p0 = fit(rint(yy), x, use)
lo, hi = np.quantile(yy, [0.05, 0.95])
t1, p1 = fit(rint(np.clip(yy, lo, hi)), x, use)
t2, p2 = fit(rint(yy), (x > 0).astype(float), use)
print(f"    RINT (as scanned)      t = {t0:+.3f}  p = {p0:.2e}")
print(f"    winsorised 5-95%       t = {t1:+.3f}  p = {p1:.2e}")
print(f"    binary carrier status  t = {t2:+.3f}  p = {p2:.2e}")
print("\ndone")
