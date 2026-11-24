#!/usr/bin/env python
"""TRANS FOLLOW-UP -- try to kill the trans signal before believing it."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SUB = f"{P1}/subproject_methylation"
WU = f"{P1}/analyses/cdr_qtl_full_panel"
CANON = f"{P1}/analyses/cdr_qtl"
OUT = f"{WU}/output/trans"
sys.path.insert(0, f"{WU}/code")
from phenotype_family import PRIMARY

NBG = int(sys.argv[1]) if len(sys.argv) > 1 else 300
SEED = 20260805
os.makedirs(OUT, exist_ok=True)


def rint(a):
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def read_grm(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    vals = np.fromfile(f"{prefix}.grm.bin", dtype=np.float32)
    G = np.zeros((n, n)); k = 0
    for i in range(n):
        for j in range(i + 1):
            G[i, j] = G[j, i] = vals[k]; k += 1
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
vi = {s: i for i, s in enumerate(V.sv_locus)}

RES = pd.read_csv(f"{OUT}/focused_trans_results.tsv", sep="\t")
HITS = RES[(RES.relation == "trans") & (RES.bh_q_trans_escalated < 0.05)].copy()
print(f"trans hits at BH q<0.05: {len(HITS)}  targets {sorted(set(HITS.target_chrom))}  "
      f"SVs {HITS.sv_locus.nunique()}")


def series(ph, chrom):
    src = PD_ if ph == "cdr_pos_offset_deepest" else PC
    s = src[src.chm13_chrom == chrom]
    return s.groupby("sample")[ph].mean().dropna() if ph in src.columns else pd.Series(dtype=float)


def alen(chrom):
    return PC[PC.chm13_chrom == chrom].groupby("sample")["active_len"].mean().dropna()


def fit(samples, y, x, extra=None):
    bc = np.array([1.0 if str(BC.get(s, "")) == "R1041" else 0.0 for s in samples])
    cols = [np.ones(len(y)), np.array([PCMAP[s] for s in samples])]
    if bc.std() > 1e-10:
        cols.append(bc[:, None])
    if extra is not None:
        cols.append(rint(extra)[:, None])
    X = np.column_stack(cols)
    df = len(y) - X.shape[1] - 1
    yr = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    xr = x - X @ np.linalg.lstsq(X, x, rcond=None)[0]
    if yr.std() < 1e-12 or xr.std() < 1e-12:
        return np.nan, np.nan, len(y)
    yr = (yr - yr.mean()) / yr.std(); xr = (xr - xr.mean()) / xr.std()
    r = float(xr @ yr / len(yr))
    t = r * np.sqrt(df) / np.sqrt(max(1 - r ** 2, 1e-30))
    return t, 2 * stats.t.sf(abs(t), df), len(y)


rows = []
for h in HITS.itertuples():
    ytr = series(h.phenotype, h.target_chrom)
    ycis = series(h.phenotype, h.sv_chrom)
    al = alen(h.target_chrom)
    col = vi[h.sv_locus]
    base = [s for s in ytr.index if s in PCMAP and s in GS]
    d = {s: Gz[GS.index(s), col] for s in base}
    base = [s for s in base if not np.isnan(d[s])]
    x0 = np.array([d[s] for s in base])
    t0, p0, n0 = fit(base, rint(ytr.loc[base].values), x0)

    sh = [s for s in base if s in ycis.index]
    t1, p1, n1 = fit(sh, rint(ytr.loc[sh].values), np.array([d[s] for s in sh]),
                     extra=ycis.loc[sh].values)
    sl = [s for s in base if s in al.index]
    t2, p2, n2 = fit(sl, rint(ytr.loc[sl].values), np.array([d[s] for s in sl]),
                     extra=al.loc[sl].values)
    both = [s for s in ytr.index if s in ycis.index]
    rr = stats.spearmanr(ytr.loc[both], ycis.loc[both]) if len(both) > 10 else None
    rows.append(dict(sv_locus=h.sv_locus, phenotype=h.phenotype, target=h.target_chrom,
                     n=n0, t_raw=t0, p_raw=p0,
                     n_cisadj=n1, t_adj_cis=t1, p_adj_cis=p1,
                     n_lenadj=n2, t_adj_alen=t2, p_adj_alen=p2,
                     rho_cis_trans_pheno=(rr.statistic if rr else np.nan),
                     p_rho=(rr.pvalue if rr else np.nan)))
F = pd.DataFrame(rows)
F.to_csv(f"{OUT}/trans_followup.tsv", sep="\t", index=False, float_format="%.6g")
print("\n=== H_SHARED / H_LENGTH: does the trans signal survive adjustment? ===")
print(F.to_string(index=False))

print(f"\n=== H_CHROM: {NBG} AF-matched random chr8 SVs through the SAME trans cells ===")
rng = np.random.default_rng(SEED)
c8 = V[(V.chrom == "chr8") & (~V.sv_locus.isin(HITS.sv_locus))].copy()
targets = sorted(set(zip(HITS.phenotype, HITS.target_chrom)))
bg = []
for ph, tc in targets:
    ytr = series(ph, tc)
    obs = HITS[(HITS.phenotype == ph) & (HITS.target_chrom == tc)]
    for afq in sorted(set(np.round(obs.af, 3))):
        pool = c8[(c8.af - afq).abs() <= 0.05]
        if len(pool) < 20:
            pool = c8.iloc[(c8.af - afq).abs().argsort()[:200]]
        pick = pool.sample(min(NBG, len(pool)), random_state=int(afq * 1e6) % 2**31)
        ts = []
        for r in pick.itertuples():
            col = vi[r.sv_locus]
            base = [s for s in ytr.index if s in PCMAP and s in GS]
            d = {s: Gz[GS.index(s), col] for s in base}
            base = [s for s in base if not np.isnan(d[s])]
            if len(base) < 40:
                continue
            t, p, n = fit(base, rint(ytr.loc[base].values), np.array([d[s] for s in base]))
            if not np.isnan(t):
                ts.append(abs(t))
        ts = np.array(ts)
        o = float(obs[np.isclose(obs.af.round(3), afq)].abs_t.max())
        bg.append(dict(phenotype=ph, target=tc, af=afq, n_bg=len(ts), obs_abs_t=o,
                       bg_median=np.median(ts) if len(ts) else np.nan,
                       bg_p95=np.quantile(ts, .95) if len(ts) else np.nan,
                       bg_max=ts.max() if len(ts) else np.nan,
                       emp_p=float((ts >= o).sum() + 1) / (len(ts) + 1) if len(ts) else np.nan))
B = pd.DataFrame(bg)
B.to_csv(f"{OUT}/trans_background.tsv", sep="\t", index=False, float_format="%.6g")
print(B.to_string(index=False))
print("\nemp_p = fraction of AF-matched chr8 SVs reaching the observed |t| in the SAME cell.")
print("done")
