#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Two multiplicity framings for \"is there a second signal at the chr8 dip locus\"."""
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
NPERM = 2000
RNG = np.random.default_rng(20260807)

vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
samples = np.array(open(WU + "/geno/full_samples.txt").read().split())

w = vt[(vt.chrom == "chr8") & (vt.pos >= LO) & (vt.pos <= HI)].copy().sort_values("pos")
w = w.reset_index()
ph = pd.read_csv(CAN + "/output/per_hap_chrom_phenotypes.tsv", sep="\t",
                 engine="python")
p8 = ph[ph.chm13_chrom == "chr8"].groupby("sample")[PHENO].mean().dropna()
keep = [i for i, s in enumerate(samples) if s in p8.index]
y = p8.loc[samples[keep]].values.astype(float)
Gw = G[w["index"].values][:, keep]

w["key"] = ["|".join("%.0f" % v if v == v else "n" for v in r) for r in Gw]
u = w.drop_duplicates("key").reset_index(drop=True)
Gu = Gw[[list(w.key).index(k) for k in u.key]]
u["label"] = ["%s-%d@%d" % (r.svtype, abs(int(r.svlen)), int(r.pos)) for r in u.itertuples()]
print("distinct genotype vectors in the window: %d   n = %d individuals" % (len(u), len(y)))

gate = pd.read_csv(WU + "/output/gates/active_len_gate.tsv", sep="\t")
g8 = gate[gate.chrom == "chr8"].drop_duplicates(subset=["sv", "af"])
print("\npre-specified chr8 leads from the scan (gates/active_len_gate.tsv):")
for _, r in g8.iterrows():
    print("   %-24s af %.3f" % (r.sv, r.af))
PRE_SV = sorted({s.split("@")[0] + "@" + s.split("@")[1] for s in g8.sv})
pre_idx = [i for i, r in u.iterrows() if r.label in PRE_SV]
print("matched to %d genotype vectors in the window: %s"
      % (len(pre_idx), [u.label[i] for i in pre_idx]))
assert len(pre_idx) >= 4, "pre-specified leads did not match the window table"


def cond_p(i, lead, yv):
    X = np.column_stack([Gu[i], Gu[lead]])
    m = np.isfinite(X).all(axis=1) & np.isfinite(yv)
    if m.sum() < 30 or np.nanstd(X[m, 0]) == 0:
        return np.nan
    return float(sm.OLS(yv[m], sm.add_constant(X[m])).fit().pvalues[1])


def marg_p(i, yv):
    m = np.isfinite(Gu[i]) & np.isfinite(yv)
    return float(sm.OLS(yv[m], sm.add_constant(Gu[i][m])).fit().pvalues[1])


lead = int(np.argmin([marg_p(i, y) for i in range(len(u))]))
print("\ntop lead by marginal p: %s (af %.3f)" % (u.label[lead], u.af[lead]))

ANGLES = {
    "A  pre-specified leads only": [i for i in pre_idx if i != lead],
    "B  open search, whole window": [i for i in range(len(u)) if i != lead],
}

print("\n" + "=" * 92)
print("%-30s %8s %-24s %10s %12s" % ("framing", "screened", "best candidate",
                                     "best p", "calibrated p"))
print("=" * 92)
res = []
for name, cand in ANGLES.items():
    ps = {i: cond_p(i, lead, y) for i in cand}
    ps = {i: p for i, p in ps.items() if p == p}
    best_i = min(ps, key=ps.get)
    obs = ps[best_i]
    cnt = 0
    for _ in range(NPERM):
        yp = RNG.permutation(y)
        b = 1.0
        for i in cand:
            p = cond_p(i, lead, yp)
            if p == p and p < b:
                b = p
        cnt += (b <= obs)
    cal = (cnt + 1) / (NPERM + 1)
    print("%-30s %8d %-24s %10.4g %12.4f" % (name, len(cand), u.label[best_i], obs, cal))
    res.append(dict(framing=name, n_screened=len(cand), best=u.label[best_i],
                    p_best=obs, p_calibrated=cal))

print("\nper pre-specified lead, conditional on the top lead:")
print("%-24s %6s %8s %10s %10s" % ("lead", "af", "|r|lead", "p_marginal", "p_conditional"))
for i in pre_idx:
    if i == lead:
        continue
    m = np.isfinite(Gu[i]) & np.isfinite(Gu[lead])
    r = np.corrcoef(Gu[i][m], Gu[lead][m])[0, 1]
    print("%-24s %6.3f %8.2f %10.3g %10.4f"
          % (u.label[i], u.af[i], abs(r), marg_p(i, y), cond_p(i, lead, y)))

out = WU + "/output/chr8_multisignal/second_signal_two_framings.tsv"
pd.DataFrame(res).to_csv(out, sep="\t", index=False, float_format="%.6g")
print("\nwrote %s" % out)
