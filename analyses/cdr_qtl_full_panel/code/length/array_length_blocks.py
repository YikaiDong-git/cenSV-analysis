#!/usr/bin/env python
import os
import runpy
import sys
import warnings

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

warnings.filterwarnings("ignore")
ROOT = CENSV_ROOT + ""
WU = ROOT + "/analyses/cdr_qtl_full_panel"
OUT = WU + "/output/length"
SRC = WU + "/code/length/array_length_budget.py"

CACHE = OUT + "/model_table.tsv"
if os.path.exists(CACHE) and "--rebuild" not in sys.argv:
    print("reusing the cached modelling table (pass --rebuild to regenerate)")
    D = pd.read_csv(CACHE, sep="\t")
    G = runpy.run_path(SRC, run_name="__adj_r2_only__") if False else None
    def adj_r2(y, X):
        X = (np.column_stack([np.ones(len(y))] + [np.asarray(c, float) for c in X])
             if len(X) else np.ones((len(y), 1)))
        keep = [0] + [j for j in range(1, X.shape[1]) if np.std(X[:, j]) > 1e-12]
        X = X[:, keep]
        n, p = X.shape
        if n <= p + 1:
            return np.nan
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        ss_res = float(np.sum((y - X @ beta) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        if ss_tot <= 0:
            return np.nan
        return 1 - (1 - (1 - ss_res / ss_tot)) * (n - 1) / (n - p)
else:
    print("running the original to obtain its modelling table ...")
    _out = sys.stdout
    sys.stdout = open("/dev/null", "w")
    G = runpy.run_path(SRC)
    sys.stdout = _out
    D, adj_r2 = G["D"], G["adj_r2"]
    D.to_csv(CACHE, sep="\t", index=False, float_format="%.8g")
print("modelling table: %d haplotype-chromosomes, %d chromosomes"
      % (len(D), D.chrom.nunique()))

BLOCKS = {"technical": ["coverage", "is_r1041", "n_contigs", "n_intervals"],
          "ancestry": ["gPC%d" % (i+1) for i in range(10)],
          "SVs": ["burden", "netbp"]}
NAMES = ["technical", "ancestry", "SVs", "lineage"]


def cols_of(g, name):
    if name == "lineage":
        d = pd.get_dummies(g.haplogroup.fillna("NA"), drop_first=True)
        return [d[k].values.astype(float) for k in d.columns]
    return [g[k].fillna(g[k].median()).values.astype(float) for k in BLOCKS[name] if k in g]


def cv_r2(y, X, k=5, seed=0):
    if not X:
        return np.nan
    A = np.column_stack([np.ones(len(y))] + X)
    keep = [0] + [j for j in range(1, A.shape[1]) if np.std(A[:, j]) > 1e-12]
    A = A[:, keep]
    n = len(y)
    if n < 3 * A.shape[1]:
        return np.nan
    rng = np.random.default_rng(seed)
    fold = rng.permutation(n) % k
    pred = np.full(n, np.nan)
    for f in range(k):
        tr, te = fold != f, fold == f
        if tr.sum() <= A.shape[1] + 1:
            return np.nan
        beta, *_ = np.linalg.lstsq(A[tr], y[tr], rcond=None)
        pred[te] = A[te] @ beta
    return float(1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2))


cv, com, lin = [], [], []
for c, g in D.groupby("chrom"):
    g = g.dropna(subset=["alen_z"])
    if len(g) < 60:
        continue
    y = g.alen_z.values.astype(float)
    blocks = {nm: cols_of(g, nm) for nm in NAMES}
    allc = [x for nm in NAMES for x in blocks[nm]]

    r = dict(chrom=c, n=len(g), n_haplogroups=int(g.haplogroup.nunique()))
    for nm in NAMES:
        r["cv_" + nm] = cv_r2(y, blocks[nm])
        r["ncol_" + nm] = len(blocks[nm])
    r["cv_all"] = cv_r2(y, allc)
    cv.append(r)

    r2_all = adj_r2(y, allc)
    u = {}
    for nm in NAMES:
        rest = [x for m in NAMES if m != nm for x in blocks[m]]
        u[nm] = r2_all - adj_r2(y, rest)
    row = dict(chrom=c, n=len(g), R2_all=r2_all,
               **{"uniq_" + k2: v for k2, v in u.items()})
    row["shared"] = r2_all - sum(u.values())
    com.append(row)

    r2_a, r2_l = adj_r2(y, blocks["ancestry"]), adj_r2(y, blocks["lineage"])
    r2_al = adj_r2(y, blocks["ancestry"] + blocks["lineage"])
    pc_by_hg = []
    if len(blocks["lineage"]):
        hgd = pd.get_dummies(g.haplogroup.fillna("NA"), drop_first=True)
        for p in BLOCKS["ancestry"]:
            ok = g[p].notna().values
            if ok.sum() >= 40:
                pc_by_hg.append(adj_r2(g[p].values.astype(float)[ok],
                                       [hgd[k].values.astype(float)[ok] for k in hgd.columns]))
    lin.append(dict(chrom=c, n=len(g), n_haplogroups=int(g.haplogroup.nunique()),
                    R2_ancestry=r2_a, R2_lineage=r2_l, R2_both=r2_al,
                    shared_anc_lin=(r2_a + r2_l - r2_al),
                    hg_predicts_gPC_mean=float(np.nanmean(pc_by_hg)) if pc_by_hg else np.nan,
                    hg_predicts_gPC_max=float(np.nanmax(pc_by_hg)) if pc_by_hg else np.nan))

CV, CM, LN = pd.DataFrame(cv), pd.DataFrame(com), pd.DataFrame(lin)
for d, f in ((CV, "variance_budget_cv"), (CM, "variance_commonality"),
             (LN, "lineage_ancestry")):
    d.to_csv("%s/%s.tsv" % (OUT, f), sep="\t", index=False, float_format="%.6g")

print()
print("=" * 104)
print("1. CROSS-VALIDATED R2 -- the same question with the df penalty measured, not assumed")
print(CV[["chrom", "n"] + ["cv_" + n for n in NAMES] + ["cv_all"]].round(3).to_string(index=False))
print()
print("  median CV R2:  " + "   ".join("%s %+.3f" % (n, CV["cv_" + n].median()) for n in NAMES)
      + "   ALL %+.3f" % CV.cv_all.median())
print("  block largest by ADJUSTED R2 (the current panel) vs by CV R2:")
adj = pd.read_csv(OUT + "/variance_budget.tsv", sep="\t")
w_adj = adj.set_index("chrom")[["R2_" + n for n in NAMES]].idxmax(axis=1).str.replace("R2_", "")
w_cv = CV.set_index("chrom")[["cv_" + n for n in NAMES]].idxmax(axis=1).str.replace("cv_", "")
print("    adjusted: %s" % w_adj.value_counts().to_dict())
print("    CV      : %s" % w_cv.value_counts().to_dict())
print("    they disagree on %d of %d chromosomes"
      % (int((w_adj.reindex(w_cv.index) != w_cv).sum()), len(w_cv)))
print("  chromosomes where NO block predicts out of sample (all CV R2 < 0): %d of %d"
      % (int((CV[["cv_" + n for n in NAMES]] < 0).all(axis=1).sum()), len(CV)))

print()
print("=" * 104)
print("2. IS 'HAPLOGROUP LINEAGE' JUST ANCESTRY?")
print(LN.round(3).to_string(index=False))
print()
print("  shared between the ancestry and lineage blocks, adjusted R2 units:")
print("    median %+.4f   over median lineage-alone %+.4f  ->  %.0f%% of the lineage block's"
      % (LN.shared_anc_lin.median(), LN.R2_lineage.median(),
         100*LN.shared_anc_lin.median()/max(1e-9, abs(LN.R2_lineage.median()))))
print("    marginal R2 is also explainable by ancestry")
print("  haplogroup dummies predicting the ancestry PCs: mean adj R2 median %+.4f, max %+.4f"
      % (LN.hg_predicts_gPC_mean.median(), LN.hg_predicts_gPC_max.max()))

print()
print("=" * 104)
print("3. AN ADDITIVE DECOMPOSITION -- unique contributions plus one shared term")
print(CM.round(3).to_string(index=False))
chk = (CM[["uniq_" + n for n in NAMES]].sum(axis=1) + CM.shared - CM.R2_all).abs().max()
print("\n  max |sum of parts - joint fit| = %.2e  (additive by construction)" % chk)
print("  median unique:  " + "   ".join("%s %+.4f" % (n, CM["uniq_" + n].median())
                                        for n in NAMES))
print("  median shared:  %+.4f" % CM.shared.median())
print("  chromosomes where the shared term exceeds every unique term: %d of %d"
      % (int((CM.shared > CM[["uniq_" + n for n in NAMES]].max(axis=1)).sum()), len(CM)))
print("\nwrote %s/{variance_budget_cv,variance_commonality,lineage_ancestry}.tsv" % OUT)
