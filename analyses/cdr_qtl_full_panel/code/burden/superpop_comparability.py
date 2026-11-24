#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
CAN = ROOT + "/analyses/cdr_qtl"
WU = ROOT + "/analyses/cdr_qtl_full_panel"
RNG = np.random.default_rng(20260807)
NPERM = 1000

PRIMARY = ["cdr_pos_offset", "cdr_pos_offset_deepest", "cdr_pos_iqr", "n_cdrs",
           "cdrs_per_mb", "cdr_fraction", "median_cdr_width", "mean_5mC"]
EXTRA = ["cdr_pos_var", "mean_5mC_covwt", "active_len", "total_cdr_bp"]
NAME = {"cdr_pos_offset": "dip position, median", "cdr_pos_offset_deepest": "dip position, deepest",
        "cdr_pos_iqr": "dip spread, IQR", "n_cdrs": "dip count", "cdrs_per_mb": "dip density",
        "cdr_fraction": "hypomethylated fraction", "median_cdr_width": "dip width",
        "mean_5mC": "array methylation level", "cdr_pos_var": "dip spread, variance",
        "mean_5mC_covwt": "methylation, cov-weighted", "active_len": "ARRAY LENGTH",
        "total_cdr_bp": "total dip bp"}


def rint(a):
    a = np.asarray(a, float)
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


corr = pd.read_csv(CAN + "/output/per_hap_chrom_phenotypes.tsv", sep="\t",
                   engine="python")
deep = pd.read_csv(CAN + "/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t",
                   engine="python")
k = ["sample", "hap_local", "chm13_chrom"]
allph = corr.merge(deep[k + ["cdr_pos_offset_deepest"]], on=k, how="left")
allph["total_cdr_bp"] = allph.cdr_fraction * allph.active_len
PHEN = PRIMARY + EXTRA
P = (allph.groupby(["sample", "chm13_chrom"])[PHEN].mean().reset_index()
     .rename(columns={"chm13_chrom": "chrom"}))

SP = pd.read_csv(WU + "/output/burden/superpop_burden_individual.tsv", sep="\t")
D = P.merge(SP[["sample", "superpop"]], on="sample", how="inner")
print("rows (individual x chromosome): %d  [%d individuals x %d chromosomes]"
      % (len(D), D["sample"].nunique(), D.chrom.nunique()))
print("individuals by superpopulation: %s"
      % D.drop_duplicates("sample").superpop.value_counts().to_dict())

for c in PHEN:
    D[c + "_z"] = D.groupby("chrom")[c].transform(
        lambda s: rint(s) if s.notna().sum() == len(s) and s.std() > 0 else np.nan)

POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
DUM = pd.get_dummies(D.superpop, prefix="sp", drop_first=True).astype(float)
D = pd.concat([D, DUM], axis=1)
SPC = list(DUM.columns)


def partial_r2(d, y, block, extra=()):
    cols = [y] + list(block) + list(extra)
    dd = d[cols + ["sample"]].replace([np.inf, -np.inf], np.nan).dropna()
    yv = dd[y].values.astype(float)
    X0 = sm.add_constant(dd[list(extra)].values.astype(float)) if extra else \
        np.ones((len(dd), 1))
    X1 = sm.add_constant(np.column_stack([dd[list(block) + list(extra)].values.astype(float)]))
    m0, m1 = sm.OLS(yv, X0).fit(), sm.OLS(yv, X1).fit()
    ss0 = float(np.sum(m0.resid ** 2))
    ss1 = float(np.sum(m1.resid ** 2))
    q = len(block)
    dfd = len(dd) - X1.shape[1]
    F = ((ss0 - ss1) / q) / (ss1 / dfd)
    return (ss0 - ss1) / ss0, F, len(dd), dd


def perm_p(dd, y, block, extra, F_obs):
    ids = dd["sample"].unique()
    lab = dd.drop_duplicates("sample").set_index("sample")[list(block)]
    cnt = 0
    for _ in range(NPERM):
        perm = dict(zip(ids, RNG.permutation(ids)))
        rep = lab.loc[[perm[s] for s in dd["sample"]]].values
        d2 = dd.copy()
        d2[list(block)] = rep
        _, F, _, _ = partial_r2(d2, y, block, extra)
        cnt += (F >= F_obs)
    return (cnt + 1) / (NPERM + 1)


def icc(dd, y):
    g = dd.groupby("sample")[y]
    mb = g.mean()
    n = g.size()
    k = float(n.mean())
    msb = float(((mb - dd[y].mean()) ** 2 * n).sum() / (len(mb) - 1))
    msw = float(g.apply(lambda s: ((s - s.mean()) ** 2).sum()).sum() / (len(dd) - len(mb)))
    v = (msb - msw) / k
    return max(0.0, v) / (max(0.0, v) + msw)


print("\n" + "=" * 106)
print("SUPERPOPULATION EFFECT, before and after adjusting for active-array length")
print("=" * 106)
print("%-26s %9s %8s %10s %9s %8s %10s  %8s" %
      ("phenotype", "R2_pop", "F", "perm P", "R2_pop|len", "F", "perm P", "reduction"))
print("-" * 106)
rows = []
for y in PHEN:
    yz = y + "_z"
    r2a, Fa, na, dda = partial_r2(D, yz, SPC)
    pa = perm_p(dda, yz, SPC, (), Fa)
    if y == "active_len":
        print("%-26s %9.4f %8.2f %10.4f %9s %8s %10s  %8s"
              % (NAME[y], r2a, Fa, pa, "-", "-", "-", "(covariate)"))
        rows.append(dict(phenotype=y, r2_pop=r2a, F=Fa, p_perm=pa,
                         r2_pop_adj=np.nan, F_adj=np.nan, p_perm_adj=np.nan, n=na))
        continue
    r2b, Fb, nb, ddb = partial_r2(D, yz, SPC, extra=("active_len_z",))
    pb = perm_p(ddb, yz, SPC, ("active_len_z",), Fb)
    print("%-26s %9.4f %8.2f %10.4f %9.4f %8.2f %10.4f  %7.0f%%"
          % (NAME[y], r2a, Fa, pa, r2b, Fb, pb, 100 * (1 - r2b / max(r2a, 1e-12))))
    rows.append(dict(phenotype=y, r2_pop=r2a, F=Fa, p_perm=pa,
                     r2_pop_adj=r2b, F_adj=Fb, p_perm_adj=pb, n=nb))

print("\n" + "=" * 106)
print("BETWEEN-INDIVIDUAL SHARE (intraclass correlation), before and after")
print("=" * 106)
print("%-26s %14s %16s" % ("phenotype", "ICC", "ICC | array length"))
print("-" * 60)
for y in PHEN:
    yz = y + "_z"
    if y == "active_len":
        dd = D[[yz, "sample"]].replace([np.inf, -np.inf], np.nan).dropna()
        a = icc(dd, yz)
        print("%-26s %14.4f %16s" % (NAME[y], a, "(covariate)"))
        for r in rows:
            if r["phenotype"] == y:
                r["icc"], r["icc_adj"] = a, np.nan
        continue
    dd = D[[yz, "active_len_z", "sample"]].replace([np.inf, -np.inf], np.nan).dropna()
    a = icc(dd, yz)
    m = sm.OLS(dd[yz].values.astype(float),
               sm.add_constant(dd["active_len_z"].values.astype(float))).fit()
    dd = dd.assign(res=m.resid)
    b = icc(dd, "res")
    print("%-26s %14.4f %16.4f" % (NAME[y], a, b))
    for r in rows:
        if r["phenotype"] == y:
            r["icc"], r["icc_adj"] = a, b

out = WU + "/output/burden/superpop_phenotype_comparability.tsv"
pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.6g")
print("\nwrote %s" % out)

print("\n--- which phenotypes are length-normalised by construction ---")
for y in PHEN:
    r = np.corrcoef(D[y + "_z"].dropna(),
                    D.loc[D[y + "_z"].notna(), "active_len_z"])[0, 1] \
        if D[y + "_z"].notna().sum() > 10 else np.nan
    print("  %-26s r with array length = %+.3f" % (NAME[y], r))
