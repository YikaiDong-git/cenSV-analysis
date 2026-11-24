#!/usr/bin/env python
"""Missing gate: test active_len ~ dosage directly, for every surviving lead."""
import numpy as np, pandas as pd, struct, os
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
CAN = CENSV_ROOT + "/analyses/cdr_phenotypes/output/compat/compat_root_span"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
LENGTH_NORMALISED = {"cdr_pos_offset", "cdr_pos_offset_deepest", "cdr_fraction", "cdrs_per_mb",
                     "cdr_pos_var", "cdr_pos_iqr"}


def read_grm(p):
    n = sum(1 for _ in open(f"{p}.grm.id")); K = np.zeros((n, n))
    with open(f"{p}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]; K[i, j] = v; K[j, i] = v
    return K, [l.split()[1] for l in open(f"{p}.grm.id")]


def rint(a): return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def slope(y, x, extra):
    X = np.column_stack([np.ones(len(y)), np.asarray(x, float), np.asarray(extra, float)])
    m = np.isfinite(X).all(1) & np.isfinite(y)
    X, y = X[m], np.asarray(y, float)[m]
    if len(y) <= X.shape[1] + 1 or X[:, 1].std() == 0:
        return None
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X) * (r ** 2).sum() / (len(y) - X.shape[1])))
    t = b[1] / se[1]
    return dict(t=float(t), p=float(2 * (1 - stats.t.cdf(abs(t), len(y) - X.shape[1]))))


ph = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
phd = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t")
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], index=gids)
G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
samp = np.array(open(f"{WU}/geno/full_samples.txt").read().split())
loc2row = {l: i for i, l in enumerate(vt.sv_locus)}
gate = pd.read_csv(f"{WU}/output/gates/full_hits_gate_1e7.tsv", sep="\t")
surv = gate[gate.overall == "SURVIVES"]

print(f"{'chrom':6s} {'SV':22s} {'phenotype':24s} {'t(len)':>8s} {'p(len)':>10s} "
      f"{'t(start)':>9s} {'t(end)':>8s}  verdict")
rows = []
for r in surv.itertuples():
    src = phd if r.phenotype == "cdr_pos_offset_deepest" else ph
    s = src[src.chm13_chrom == r.chrom]
    agg = s.groupby("sample").agg(y=(r.phenotype, "mean"))
    b = ph[ph.chm13_chrom == r.chrom].groupby("sample").agg(
        aa_s=("active_start", "mean"), aa_e=("active_end", "mean"), aa_l=("active_len", "mean"))
    bc = ph[ph.chm13_chrom == r.chrom].groupby("sample")["basecaller"].agg(
        lambda v: 1.0 if v.astype(str).str.contains("1041").any() else 0.0)
    agg = agg.join(b).join(bc.rename("bc"))
    dos = pd.Series(G[loc2row[r.sv_locus], :], index=samp).replace(-1, np.nan)
    keep = sorted(set(agg.dropna(subset=["y"]).index) & set(dos.dropna().index) & set(PCs.index))
    d = dos.loc[keep].values
    COV = np.column_stack([PCs.loc[keep].values, agg.loc[keep, "bc"].values[:, None]])
    L = slope(rint(agg.loc[keep, "aa_l"].values), d, COV)
    S = slope(rint(agg.loc[keep, "aa_s"].values), d, COV)
    E = slope(rint(agg.loc[keep, "aa_e"].values), d, COV)
    lennorm = r.phenotype in LENGTH_NORMALISED
    bad = lennorm and L is not None and L["p"] < 0.05
    verdict = "LENGTH MOVES -> mechanical" if bad else ("len stable" if lennorm else "n/a (not normalised)")
    print(f"{r.chrom:6s} {r.svtype+'-'+str(int(r.svlen))+'@'+str(int(r.pos)):22s} {r.phenotype:24s} "
          f"{L['t']:+8.2f} {L['p']:10.2e} {S['t']:+9.2f} {E['t']:+8.2f}  {verdict}")
    rows.append(dict(chrom=r.chrom, sv=f"{r.svtype}-{int(r.svlen)}@{int(r.pos)}", af=r.af,
                     phenotype=r.phenotype, p_scan=r.p_best, t_active_len=L["t"], p_active_len=L["p"],
                     t_active_start=S["t"], t_active_end=E["t"],
                     length_normalised=lennorm, new_verdict=("DEMOTED_length_moves" if bad else "SURVIVES")))

out = pd.DataFrame(rows)
p = CENSV_ROOT + "/analyses/cdr_phenotypes/output/per_frame/frames/gates_span/active_len_gate.tsv"
out.to_csv(p, sep="\t", index=False, float_format="%.6g")
print("\n=== revised survivor count ===")
print(out.groupby(["chrom", "new_verdict"]).size().to_string())
print(f"\nwrote {p}")
