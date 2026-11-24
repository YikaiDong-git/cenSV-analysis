#!/usr/bin/env python3
"""Boundary-ASYMMETRY test (closes humility item #4) for the chr8 headline cluster."""
import os, struct, numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
ROOT = CENSV_ROOT + ""
DISC = f"{ROOT}/analyses/cdr_qtl"
SUB = f"{ROOT}/subproject_methylation"
PHENO = f"{DISC}/output/per_hap_chrom_phenotypes.tsv"
GRM = f"{SUB}/results/greml/pericen_1m_all"
SVDROOT = f"{ROOT}/sv_extraction/decompose/mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
ANCHORS = [
    ("chr8", 44300221, "DEL", 113, "DEL-113 (headline)"),
    ("chr8", 44302147, "INS", 59, "INS-59"),
    ("chr8", 44460941, "DEL", 347, "DEL-347"),
    ("chr8", 44305108, "DEL", 338, "DEL-338"),
]

def read_grm_pcs(prefix, npc=10):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]; K[i, j] = v; K[j, i] = v
    ids = [l.split()[1] for l in open(f"{prefix}.grm.id")]
    ev, evec = np.linalg.eigh(K); o = np.argsort(ev)[::-1]
    P = pd.DataFrame(evec[:, o[:npc]], columns=[f"PC{i+1}" for i in range(npc)]); P["sample"] = ids
    return P.set_index("sample")

def rint(a): r = stats.rankdata(a); return stats.norm.ppf((r - 0.5) / len(r))

def ols(y, Xcols):
    X = np.column_stack([np.ones(len(y))] + Xcols)
    m = ~np.isnan(y) & ~np.isnan(X).any(axis=1)
    Xm, ym, n, k = X[m], y[m], int(m.sum()), X.shape[1]
    b = np.linalg.lstsq(Xm, ym, rcond=None)[0]; r = ym - Xm @ b
    s2 = (r ** 2).sum() / (n - k); se = np.sqrt(np.diag(s2 * np.linalg.inv(Xm.T @ Xm)))
    t = b / se; p = 2 * (1 - stats.t.cdf(np.abs(t), n - k))
    return float(t[1]), float(p[1]), n

ph = pd.read_csv(PHENO, sep="\t")
PCs = read_grm_pcs(GRM)
print(f"{'anchor':20} {'start~dos':>16} {'end~dos':>16} {'center~dos':>16} {'len~dos':>14} {'offset~dos(target)':>20}  verdict")
for chrom, pos, svt, svl, lab in ANCHORS:
    phc = ph[ph.chm13_chrom == chrom]
    ind = phc.groupby("sample").agg(off=("cdr_pos_offset","mean"), st=("active_start","mean"),
                                    en=("active_end","mean"), ln=("active_len","mean")).reset_index()
    ind["ctr"] = (ind["st"] + ind["en"]) / 2
    npz = np.load(f"{SVDROOT}/{chrom}/{chrom}_sv_genotypes.npz")
    G = npz[list(npz.keys())[0]]; samp = open(f"{SVDROOT}/{chrom}/{chrom}_sv_samples.txt").read().split()
    if G.shape[0] != len(samp) and G.shape[1] == len(samp): G = G.T
    svv = pd.read_csv(f"{SVDROOT}/{chrom}/{chrom}_sv_variants.tsv", sep="\t")
    cand = svv[(svv.pos == pos) & (svv.svtype.astype(str).str.startswith(svt)) & (svv.svlen.between(svl-3, svl+3))]
    dos = pd.Series(np.rint(G[:, cand.index[0]].astype(float)), index=samp).clip(0, 2)
    d = ind.merge(dos.rename("dos"), left_on="sample", right_index=True).merge(PCs, left_on="sample", right_index=True)
    pcm = [d[f"PC{i+1}"].values for i in range(10)]; dv = d["dos"].values
    ts = ols(rint(d["st"].values), [dv]+pcm); te = ols(rint(d["en"].values), [dv]+pcm)
    tc = ols(rint(d["ctr"].values), [dv]+pcm); tl = ols(rint(d["ln"].values), [dv]+pcm)
    to = ols(rint(d["off"].values), [dv]+pcm)
    edges_stable = ts[1] > 0.05 and te[1] > 0.05
    verdict = "REAL reposition (frame stable)" if edges_stable else "FRAME-SHIFT suspect (an edge moves)"
    f = lambda r: f"t={r[0]:+.2f} p={r[1]:.2f}"
    print(f"{lab:20} {f(ts):>16} {f(te):>16} {f(tc):>16} {f(tl):>14} {f(to):>20}  {verdict}")
