#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Array-length gate on the loci the Gao unit newly makes study-wide significant."""
import struct
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
U = ROOT + "/analyses/cdr_definition"
WU = ROOT + "/analyses/cdr_qtl_full_panel"
GRM = ROOT + "/subproject_methylation/results/greml/pericen_1m_all"
ALPHA = 0.05 / (11129 * 8)
LENGTH_NORMALISED = {"cdr_pos_offset", "cdr_pos_offset_deepest", "cdr_fraction", "cdrs_per_mb",
                     "cdr_pos_var", "cdr_pos_iqr"}


def read_grm(p):
    n = sum(1 for _ in open(p + ".grm.id"))
    K = np.zeros((n, n))
    with open(p + ".grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]
                K[i, j] = v
                K[j, i] = v
    return K, [l.split()[1] for l in open(p + ".grm.id")]


def slope(y, x, extra):
    X = np.column_stack([np.ones(len(y)), np.asarray(x, float), np.asarray(extra, float)])
    y = np.asarray(y, float)
    m = np.isfinite(X).all(1) & np.isfinite(y)
    X, y = X[m], y[m]
    if len(y) <= X.shape[1] + 1 or X[:, 1].std() == 0:
        return None
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X) * (r ** 2).sum() / (len(y) - X.shape[1])))
    t = float(b[1] / se[1])
    return dict(t=t, p=float(2 * (1 - stats.t.cdf(abs(t), len(y) - X.shape[1]))), n=len(y))


ph = pd.read_csv(U + "/output/compat/compat_dip_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], index=gids)
G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
samp = np.array(open(WU + "/geno/full_samples.txt").read().split())
loc2row = {l: i for i, l in enumerate(vt.sv_locus)}

M = pd.read_csv(U + "/output/gw/escalated_1e7_both_arms.tsv", sep="\t")
gcol = "geno_group_dip" if "geno_group_dip" in M.columns else None
new = M[(M.p_emp_perm_gao < ALPHA) & (M.p_emp_perm_dip >= ALPHA)].copy()
print("loci significant in the Gao arm ONLY: %d rows, %d genotype groups"
      % (len(new), new.groupby(["chrom", "phenotype"] + ([gcol] if gcol else [])).ngroups))
print()

rows = []
seen = set()
for _, r in new.sort_values(["chrom", "phenotype", "pos"]).iterrows():
    key = (r.chrom, r.phenotype, r[gcol] if gcol else r.pos)
    if key in seen:
        continue
    seen.add(key)
    if r.sv_locus not in loc2row:
        print("  [skip] %s not in the full genotype table" % r.sv_locus)
        continue
    dose = pd.Series(G[loc2row[r.sv_locus]], index=samp)
    c = ph[ph.chm13_chrom == r.chrom]
    per = c.groupby("sample").agg(active_len=("active_len", "mean"),
                                  active_start=("active_start", "mean"),
                                  active_end=("active_end", "mean"))
    keep = sorted(set(per.index) & set(PCs.index) & set(dose.dropna().index))
    x = dose.loc[keep].to_numpy()
    P = PCs.loc[keep].to_numpy()
    res = {k: slope(per.loc[keep, k].to_numpy(), x, P)
           for k in ("active_len", "active_start", "active_end")}
    if res["active_len"] is None:
        print("  [skip] %s: degenerate design" % r.sv_locus)
        continue
    disq = res["active_len"]["p"] < 0.05
    rows.append(dict(chrom=r.chrom, pos=int(r.pos), phenotype=r.phenotype,
                     sv_locus=r.sv_locus,
                     geno_group=(r[gcol] if gcol else np.nan),
                     t_gao=r.t_obs_gao, p_gao=r.p_emp_perm_gao,
                     t_dip=r.t_obs_dip, p_dip=r.p_emp_perm_dip,
                     n=res["active_len"]["n"],
                     t_len=res["active_len"]["t"], p_len=res["active_len"]["p"],
                     t_start=res["active_start"]["t"], t_end=res["active_end"]["t"],
                     length_normalised=r.phenotype in LENGTH_NORMALISED,
                     verdict="DISQUALIFIED" if (disq and r.phenotype in LENGTH_NORMALISED)
                             else "SURVIVES"))

D = pd.DataFrame(rows)
print("%-7s %-14s %-34s %8s %9s %9s %8s %8s  %s"
      % ("chrom", "phenotype", "sv_locus", "t_gao", "t(len)", "p(len)", "t(start)", "t(end)",
         "verdict"))
print("-" * 118)
for _, r in D.iterrows():
    print("%-7s %-14s %-34s %+8.3f %+9.3f %9.3g %+8.3f %+8.3f  %s"
          % (r.chrom, r.phenotype, r.sv_locus[:34], r.t_gao, r.t_len, r.p_len,
             r.t_start, r.t_end, r.verdict))

D.to_csv(U + "/output/gw/new_candidates_length_gate.tsv", sep="\t", index=False,
         float_format="%.6g")
print()
print("SURVIVES %d / DISQUALIFIED %d  (of %d distinct candidate cells)"
      % (int((D.verdict == "SURVIVES").sum()), int((D.verdict == "DISQUALIFIED").sum()), len(D)))
print("[written] %s/output/gw/new_candidates_length_gate.tsv" % U)

lead = "chr8_sv30_DEL_mon_best8_midpoint1000"
dose = pd.Series(G[loc2row[lead]], index=samp)
c = ph[ph.chm13_chrom == "chr8"]
per = c.groupby("sample").agg(active_len=("active_len", "mean"))
keep = sorted(set(per.index) & set(PCs.index) & set(dose.dropna().index))
s = slope(per.loc[keep, "active_len"].to_numpy(), dose.loc[keep].to_numpy(),
          PCs.loc[keep].to_numpy())
print()
print("CONTROL  the manuscript's chr8 lead against the same gate: t(len) = %+.3f, p = %.3f -> %s"
      % (s["t"], s["p"], "SURVIVES" if s["p"] >= 0.05 else "DISQUALIFIED"))
assert s["p"] >= 0.05, \
    ("the chr8 lead FAILS the length gate in this implementation (p = %.4f); the gate is "
     "mis-specified and the verdicts above cannot be read" % s["p"])
print("         gate reproduces the published verdict for the lead -- the table above is readable")
