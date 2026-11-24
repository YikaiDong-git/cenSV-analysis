#!/usr/bin/env python
"""Second negative control for the local-methylation-window scan (met_08b_cdr)."""
import os
import struct

import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
CDR = f"{P1}/analyses/cdr_qtl_full_panel"
AUD = f"{P1}/analyses/cdr_calls"
GRM = f"{P1}/subproject_methylation/results/greml/pericen_1m_all"
MOUT = f"{CDR}/output/mqtl"
U = f"{P1}/analyses/meqtl_replication"
OUT = f"{U}/output/arms/window_availability"

MIN_IND = 50
NPERM = 10_000
CHUNK = 2_000
SEED = 20260901
RNG = np.random.default_rng(SEED)
SUF, RDIR = "_bam", "bam"
CHI2_MED = stats.chi2.ppf(0.5, 1)

os.makedirs(OUT, exist_ok=True)


def read_grm(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]
                K[i, j] = v
                K[j, i] = v
    return K, [l.split()[1] for l in open(f"{prefix}.grm.id")]


def fl_permutation(y, g, Z, nperm=NPERM):
    Zi = np.linalg.pinv(Z)
    ry = y - Z @ (Zi @ y)
    rg = g - Z @ (Zi @ g)
    denom = float(rg @ rg)
    if denom <= 1e-12:
        return np.nan, np.nan
    n, k = Z.shape
    dof = n - k - 1
    if dof < 1:
        return np.nan, np.nan
    S = float(ry @ ry)
    if S <= 1e-12:
        return np.nan, np.nan

    def t_of(c):
        c = np.asarray(c, float)
        s2 = (S - c * c / denom) / dof
        out = np.full(c.shape, np.nan)
        ok = s2 > 0
        out[ok] = (c[ok] / np.sqrt(denom)) / np.sqrt(s2[ok])
        return out

    t_obs = float(t_of(np.array([float(ry @ rg)]))[0])
    if not np.isfinite(t_obs):
        return np.nan, np.nan
    cnt, done = 0, 0
    base = np.tile(ry, (CHUNK, 1))
    while done < nperm:
        m = min(CHUNK, nperm - done)
        P = RNG.permuted(base[:m], axis=1)
        tp = t_of(P @ rg)
        cnt += int(np.sum(np.abs(tp) >= abs(t_obs) - 1e-12))
        done += m
    return t_obs, (cnt + 1) / (nperm + 1)


G = np.load(f"{CDR}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
vt = pd.read_csv(f"{CDR}/geno/full_variants.tsv", sep="\t", engine="python")
gsamp = np.array(open(f"{CDR}/geno/full_samples.txt").read().split())
gidx = {s: i for i, s in enumerate(gsamp)}
locus_row = {l: i for i, l in enumerate(vt.sv_locus.values)}

coh = pd.read_csv(f"{AUD}/output/cohort_corrected_final.tsv", sep="\t")
coh["is_r1041"] = coh.chemistry.astype(str).str.contains("1041").astype(float)
QI = coh.groupby("sample").agg(n_contigs=("n_contigs", "sum"),
                               is_r1041=("is_r1041", "max")).reset_index()
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCS = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]],
                   columns=[f"gPC{i+1}" for i in range(10)], index=gids)
COV = QI.set_index("sample").join(PCS, how="inner")
print(f"covariate table: {COV.shape[0]} individuals x {COV.shape[1]} covariates", flush=True)

rows, summary = [], []
for W in (10, 50, 100):
    Mth = pd.read_csv(f"{MOUT}/matrices/{RDIR}/matrix_meth_{W}kb{SUF}.tsv.gz",
                      sep="\t", index_col=0)
    print(f"\n=== +/-{W} kb: {Mth.shape[0]:,} SV loci x {Mth.shape[1]} individuals ===",
          flush=True)
    cols = [s for s in Mth.columns if s in gidx and s in COV.index]
    print(f"    individuals with genotype AND covariates: {len(cols)}", flush=True)
    Zc_all = COV.loc[cols].values.astype(float)
    n_deg = n_run = 0
    for sv in Mth.index:
        if sv not in locus_row:
            continue
        gvec = G[locus_row[sv], [gidx[s] for s in cols]]
        keep = np.isfinite(gvec)
        if keep.sum() < MIN_IND:
            continue
        gv = gvec[keep]
        if np.std(gv) < 1e-9:
            continue
        vals = Mth.loc[sv, cols].values.astype(float)[keep]
        y = np.isfinite(vals).astype(float)
        n = len(y)
        frac = float(y.mean())
        if y.std() < 1e-9:
            n_deg += 1
            rows.append(dict(sv_locus=sv, window_kb=W, n=n, frac_measured=frac,
                             t_avail=np.nan, p_avail=np.nan, degenerate=True))
            continue
        Z = np.column_stack([np.ones(n), Zc_all[keep]])
        t_a, p_a = fl_permutation(y, gv, Z)
        rows.append(dict(sv_locus=sv, window_kb=W, n=n, frac_measured=frac,
                         t_avail=t_a, p_avail=p_a, degenerate=False))
        n_run += 1
        if n_run % 500 == 0:
            print(f"    {n_run:,} tested", flush=True)

    D = pd.DataFrame([r for r in rows if r["window_kb"] == W and not r["degenerate"]])
    t = D.t_avail.values[np.isfinite(D.t_avail.values)]
    lam = float(np.median(t ** 2) / CHI2_MED) if len(t) else np.nan
    sig = D.p_avail <= 0.05
    try:
        SC = pd.read_csv(f"{MOUT}/scans/{RDIR}/scan_{W}kb{SUF}.tsv", sep="\t")
        cov_frac = float((SC.p_cov <= 0.05).mean())
        cov_n = len(SC)
    except FileNotFoundError:
        cov_frac, cov_n = np.nan, 0
    summary.append(dict(
        window_kb=W, n_loci_with_missingness=int(len(D)), n_loci_fully_measured=n_deg,
        median_frac_measured=float(D.frac_measured.median()) if len(D) else np.nan,
        avail_p05_n=int(sig.sum()),
        avail_p05_frac=float(sig.mean()) if len(D) else np.nan,
        avail_p05_frac_negative=float((D.loc[sig, "t_avail"] < 0).mean()) if int(sig.sum()) else np.nan,
        lambda_gc_avail=lam,
        depth_control_p05_frac=cov_frac, depth_control_n=cov_n))
    print(f"    tested {len(D):,} loci with missingness "
          f"({n_deg:,} fully measured, untestable)", flush=True)
    print(f"    dosage predicts MEASUREMENT PRESENCE (p<=0.05): "
          f"{int(sig.sum()):,} ({100*sig.mean():.1f}%)   lambda_GC = {lam:.3f}", flush=True)
    print(f"    depth control on the same windows (p<=0.05): {100*cov_frac:.1f}%", flush=True)

R = pd.DataFrame(rows).merge(
    vt[["sv_locus", "chrom", "pos", "svtype", "svlen", "af"]], on="sv_locus", how="left")
R.to_csv(f"{OUT}/w1_availability_control.tsv", sep="\t", index=False, float_format="%.6g")
S = pd.DataFrame(summary)
S.to_csv(f"{OUT}/w2_summary.tsv", sep="\t", index=False, float_format="%.6g")
print("\n=== SUMMARY ===", flush=True)
print(S.to_string(index=False), flush=True)
print(f"\nwrote {OUT}/w1_availability_control.tsv and w2_summary.tsv", flush=True)
