#!/usr/bin/env python
"""Does the chr9 cdr_fraction signal have a mechanical account -- and does it survive dropping carriers?"""
import numpy as np, pandas as pd, struct, os
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
CAN = f"{ROOT}/analyses/cdr_qtl"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
CHROM = "chr9"
PHENOS = ["cdr_fraction", "total_cdr_bp", "active_len", "n_cdrs", "median_cdr_width", "cdrs_per_mb"]


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
    return dict(t=float(t), p=float(2 * (1 - stats.t.cdf(abs(t), len(y) - X.shape[1]))), n=int(len(y)))


ph = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
p9 = ph[ph.chm13_chrom == CHROM]
per_ind = p9.groupby("sample")[PHENOS].mean()
bc = p9.groupby("sample")["basecaller"].agg(lambda v: 1.0 if v.astype(str).str.contains("1041").any() else 0.0)
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], index=gids)
G = np.load(f"{WU}/geno/perchrom/{CHROM}/{CHROM}_sv_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(f"{WU}/geno/perchrom/{CHROM}/{CHROM}_sv_variants.tsv", sep="\t")
svs = open(f"{WU}/geno/perchrom/{CHROM}/{CHROM}_sv_samples.txt").read().split()
svi = {s: i for i, s in enumerate(svs)}
gate = pd.read_csv(f"{WU}/output/gates/full_hits_gate_1e7.tsv", sep="\t")
leads = gate[(gate.overall == "SURVIVES") & (gate.chrom == CHROM)]

rows = []
for L in leads.itertuples():
    idx = vt[vt.sv_locus == L.sv_locus].index
    if not len(idx):
        print(f"{L.sv_locus} not found"); continue
    idx = idx[0]
    dos = pd.Series(G[idx, :], index=svs).replace(-1, np.nan)
    keep = sorted(set(per_ind.dropna(subset=["cdr_fraction"]).index) & set(dos.dropna().index)
                  & set(PCs.index) & set(bc.index))
    d = dos.loc[keep].values
    COV = np.column_stack([PCs.loc[keep].values, bc.loc[keep].values[:, None]])
    ncar = int((d > 0).sum())
    print("=" * 94)
    print(f"=== {CHROM} {L.svtype}-{L.svlen} @{int(L.pos):,}  af={L.af:.4f}  n={len(keep)}  "
          f"carriers={ncar}  (scan p={L.p_best:.1e}) ===")
    print(f"  {'phenotype':18s} {'t':>7s} {'p':>10s}   direction   role")
    role = {"cdr_fraction": "THE HIT", "total_cdr_bp": "NUMERATOR", "active_len": "DENOMINATOR",
            "n_cdrs": "numerator component", "median_cdr_width": "numerator component",
            "cdrs_per_mb": "related (count/len)"}
    res = {}
    for phn in PHENOS:
        r = slope(rint(per_ind.loc[keep, phn].values), d, COV)
        if r is None:
            print(f"  {phn:18s}    n/a"); continue
        res[phn] = r
        print(f"  {phn:18s} {r['t']:+7.2f} {r['p']:10.2e}   {'up' if r['t']>0 else 'down':9s}   {role[phn]}")
        rows.append(dict(sv=f"{L.svtype}-{L.svlen}@{int(L.pos)}", af=float(L.af), carriers=ncar,
                         phenotype=phn, t=r["t"], p=r["p"], n=r["n"]))
    if "cdr_fraction" in res and "total_cdr_bp" in res and "active_len" in res:
        sf, sn, sd = np.sign(res["cdr_fraction"]["t"]), np.sign(res["total_cdr_bp"]["t"]), np.sign(res["active_len"]["t"])
        num_ok = (sf == sn) and res["total_cdr_bp"]["p"] < 0.05
        den_ok = (sf == -sd) and res["active_len"]["p"] < 0.05
        print(f"  --> numerator explains it: {'YES' if num_ok else 'no'}"
              f"   |  denominator explains it: {'YES' if den_ok else 'no'}"
              f"   |  {'MECHANICALLY ACCOUNTED' if (num_ok or den_ok) else 'NO ACCOUNT -- treat with caution'}")
    y = rint(per_ind.loc[keep, "cdr_fraction"].values)
    base = slope(y, d, COV)
    ts = []
    carriers = [i for i, v in enumerate(d) if v > 0]
    for i in carriers:
        m = np.ones(len(d), bool); m[i] = False
        r = slope(rint(per_ind.loc[[k for k, f in zip(keep, m) if f], "cdr_fraction"].values), d[m], COV[m])
        if r:
            ts.append(r["t"])
    if ts:
        ts = np.array(ts)
        print(f"  leave-one-CARRIER-out (n={len(ts)}): full t={base['t']:+.2f}  "
              f"range [{ts.min():+.2f}, {ts.max():+.2f}]  weakest |t|={np.abs(ts).min():.2f}  "
              f"still |t|>3 in {100*np.mean(np.abs(ts)>3):.0f}% of drops")

out = f"{WU}/output/chr9/chr9_numerator_check.tsv"
pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.6g")
print(f"\nwrote {out}")
