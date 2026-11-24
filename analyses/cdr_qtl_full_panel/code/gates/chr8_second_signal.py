#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Is there a SECOND, separable signal at the chromosome-8 dip locus?"""
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = ROOT + "/analyses/cdr_qtl_full_panel"
CAN = ROOT + "/analyses/cdr_qtl"
LO, HI = 44_200_000, 44_500_000
PHENO = "cdr_pos_offset"

vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
samples = np.array(open(WU + "/geno/full_samples.txt").read().split())
assert G.shape == (len(vt), len(samples)), (G.shape, len(vt), len(samples))

w = vt[(vt.chrom == "chr8") & (vt.pos >= LO) & (vt.pos <= HI)].copy()
w = w.sort_values("pos").reset_index()
print("catalogued chr8 variants in %s-%s: %d  (%d distinct genotype groups)"
      % (f"{LO:,}", f"{HI:,}", len(w), w.geno_group.nunique()))

ph = pd.read_csv(CAN + "/output/per_hap_chrom_phenotypes.tsv", sep="\t",
                 engine="python")
p8 = ph[ph.chm13_chrom == "chr8"].groupby("sample")[PHENO].mean().dropna()
keep = [i for i, s in enumerate(samples) if s in p8.index]
sk = samples[keep]
y = p8.loc[sk].values.astype(float)
Gw = G[w["index"].values][:, keep]
print("individuals with chr8 dip position: %d" % len(y))

w["key"] = ["|".join("%.0f" % v if v == v else "n" for v in row) for row in Gw]
u = w.drop_duplicates("key").reset_index(drop=True)
Gu = Gw[[list(w.key).index(k) for k in u.key]]
u["label"] = ["%s-%d@%d" % (r.svtype, abs(int(r.svlen)), int(r.pos)) for r in u.itertuples()]
print("distinct genotype vectors: %d" % len(u))


def assoc(idx, cond=()):
    cols = [Gu[idx]] + [Gu[c] for c in cond]
    X = np.column_stack(cols)
    m = np.isfinite(X).all(axis=1) & np.isfinite(y)
    if m.sum() < 30 or np.nanstd(X[m, 0]) == 0:
        return np.nan, np.nan, int(m.sum()), np.nan
    fit = sm.OLS(y[m], sm.add_constant(X[m])).fit()
    if len(cond):
        r2 = sm.OLS(X[m, 0], sm.add_constant(X[m, 1:])).fit().rsquared
        vif = 1 / max(1e-12, 1 - r2)
    else:
        vif = 1.0
    return fit.tvalues[1], fit.pvalues[1], int(m.sum()), vif


print("\n=== 1. marginal association with dip position ===")
res = []
for i, r in u.iterrows():
    t, p, n, _ = assoc(i)
    res.append(dict(i=i, label=r.label, af=r.af, cv_r2=r.cv_r2, t=t, p=p, n=n))
R = pd.DataFrame(res).sort_values("p")
print(R.head(12).to_string(index=False, float_format=lambda v: "%.4g" % v))

lead = int(R.iloc[0].i)
print("\nlead by marginal p: %s  (t = %.2f, p = %.3g, af = %.3f)"
      % (u.label[lead], R.iloc[0].t, R.iloc[0].p, u.af[lead]))

print("\n=== 2. LD with the lead ===")
ld = []
for i in range(len(u)):
    m = np.isfinite(Gu[i]) & np.isfinite(Gu[lead])
    r = np.corrcoef(Gu[i][m], Gu[lead][m])[0, 1] if m.sum() > 30 and \
        np.std(Gu[i][m]) > 0 else np.nan
    ld.append(r)
u["r_lead"] = ld

print("\n=== 3. each variant conditioned on the lead ===")
rows = []
for i, r in u.iterrows():
    if i == lead:
        continue
    tm, pm, n, _ = assoc(i)
    tc, pc, nc, vif = assoc(i, cond=(lead,))
    rows.append(dict(label=r.label, af=r.af, cv_r2=r.cv_r2, r_lead=r.r_lead,
                     t_marg=tm, p_marg=pm, t_cond=tc, p_cond=pc, vif=vif, n=nc))
C = pd.DataFrame(rows).sort_values("p_cond")
print(C.head(14).to_string(index=False, float_format=lambda v: "%.4g" % v))

print("\n=== 4. forward selection (add while p_cond < 0.05) ===")
sel = [lead]
while True:
    best, bp = None, 1.0
    for i in range(len(u)):
        if i in sel:
            continue
        t, p, n, vif = assoc(i, cond=tuple(sel))
        if p == p and p < bp:
            best, bp, bt, bvif = i, p, t, vif
    if best is None or bp >= 0.05:
        break
    sel.append(best)
    print("  + %-22s p = %.3g  t = %+.2f  VIF = %.2f  |r| with lead = %.2f"
          % (u.label[best], bp, bt, bvif, abs(u.r_lead[best])))
print("independent signals supported by the window: %d  (%s)"
      % (len(sel), ", ".join(u.label[i] for i in sel)))

print("\n=== 5. joint model of the selected set ===")
X = np.column_stack([Gu[i] for i in sel])
m = np.isfinite(X).all(axis=1) & np.isfinite(y)
fit = sm.OLS(y[m], sm.add_constant(X[m])).fit()
for k, i in enumerate(sel):
    r2 = sm.OLS(X[m, k], sm.add_constant(np.delete(X[m], k, axis=1))).fit().rsquared \
        if len(sel) > 1 else 0.0
    print("  %-22s beta = %+.4f  t = %+.2f  p = %.3g  VIF = %.2f"
          % (u.label[i], fit.params[k + 1], fit.tvalues[k + 1], fit.pvalues[k + 1],
             1 / max(1e-12, 1 - r2)))
print("  joint R2 = %.4f on n = %d" % (fit.rsquared, int(m.sum())))

print("\n=== 6. is the second signal real, or the best of 70 screened vectors? ===")
NPERM = 1000
RNG = np.random.default_rng(20260807)
obs_best = float(C.p_cond.min())
cand = [i for i in range(len(u)) if i != lead]
Xl = Gu[lead]
cnt = 0
for _ in range(NPERM):
    yp = RNG.permutation(y)
    best = 1.0
    for i in cand:
        X = np.column_stack([Gu[i], Xl])
        m = np.isfinite(X).all(axis=1)
        if m.sum() < 30 or np.nanstd(X[m, 0]) == 0:
            continue
        f = sm.OLS(yp[m], sm.add_constant(X[m])).fit()
        if f.pvalues[1] < best:
            best = float(f.pvalues[1])
    cnt += (best <= obs_best)
pcal = (cnt + 1) / (NPERM + 1)
print("  best conditional p observed        : %.4g  (%s)" % (obs_best, C.iloc[0].label))
print("  candidates screened                : %d" % len(cand))
print("  Bonferroni threshold over them     : %.2g" % (0.05 / len(cand)))
print("  permutation-calibrated p           : %.4f  (%d of %d nulls at least this extreme)"
      % (pcal, cnt, NPERM))
print("  VERDICT: %s" % ("a second signal is supported" if pcal < 0.05 else
                         "NOT supported once the screen is accounted for"))

out = WU + "/output/chr8_multisignal/second_signal_forward_selection.tsv"
C.to_csv(out, sep="\t", index=False, float_format="%.6g")
print("\nwrote %s" % out)
