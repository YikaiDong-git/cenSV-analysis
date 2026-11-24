#!/usr/bin/env python
# -*- coding: utf-8 -*-
import importlib.util
import sys
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
U = ROOT + "/analyses/cdr_definition"
WU = ROOT + "/analyses/cdr_qtl_full_panel"
SCAN = (ROOT + "/analyses/cdr_qtl/code/scan/"
        "b6_cdr_pos_qtl_perchrom.py")
GRM = ROOT + "/subproject_methylation/results/greml/pericen_1m_all"
CHROM = "chr8"
LEAD_POS = 44300221
T_PUBLISHED = -6.173024844715221
T_GAO = -6.172461
GAP = 150_000
NREP = 1000
SEED = 20260831

spec = importlib.util.spec_from_file_location("scn", SCAN)
scn = importlib.util.module_from_spec(spec)
sys.modules["scn"] = scn
spec.loader.exec_module(scn)


def rint(a):
    r = stats.rankdata(a)
    return stats.norm.ppf((r - 0.5) / len(r))


G = np.load(WU + "/geno/perchrom/%s/%s_sv_genotypes.npz" % (CHROM, CHROM))
sv = G[list(G.keys())[0]]
sv_samples = open(WU + "/geno/perchrom/%s/%s_sv_samples.txt" % (CHROM, CHROM)).read().split()
sv_var = pd.read_csv(WU + "/geno/perchrom/%s/%s_sv_variants.tsv" % (CHROM, CHROM), sep="\t")
if sv.shape[0] != len(sv_samples):
    sv = sv.T
K, grm_samples = scn.read_gcta_grm_bin(GRM)
ev, evec = np.linalg.eigh(K)
PCs = evec[:, np.argsort(ev)[::-1][:10]]
pc = pd.DataFrame(PCs, columns=["PC%d" % (i + 1) for i in range(10)])
pc["sample"] = grm_samples

ph = pd.read_csv(U + "/output/compat/compat_dip_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
pc8 = ph[ph.chm13_chrom == CHROM]
bc_map = pc8.groupby("sample")["basecaller"].first().to_dict()

calls = pd.read_csv(U + "/output/phenotypes/per_call_positions_dip.tsv", sep="\t")
calls = calls[calls.chm13_chrom == CHROM].copy()
calls = calls.sort_values(["sample", "hap_local", "asm_contig", "call_start"])

gid, out = [], []
for (s, h, c), g in calls.groupby(["sample", "hap_local", "asm_contig"], sort=False):
    k, prev_end = 0, None
    for _, r in g.iterrows():
        if prev_end is not None and r.call_start - prev_end > GAP:
            k += 1
        prev_end = max(prev_end, r.call_end) if prev_end is not None else r.call_end
        gid.append(k)
    out.append(g)
calls = pd.concat(out)
calls["cdr_group"] = gid
ngrp = calls.groupby(["sample", "hap_local", "asm_contig"]).cdr_group.nunique().sum()
print("chr8 dips: %d in %d merged CDRs over %d haplotypes"
      % (len(calls), ngrp, calls.groupby(["sample", "hap_local"]).ngroups))

lead_j = sv_var.index[(sv_var.chrom == CHROM) & (sv_var.pos == LEAD_POS)
                      & (sv_var.sv_locus == "chr8_sv30_DEL_mon_best8_midpoint1000")]
assert len(lead_j) == 1, "lead SV not uniquely found: %d matches" % len(lead_j)
lead_j = int(lead_j[0])


def t_for(pos_by_hap):
    per_ind = pos_by_hap.groupby(level=0).mean()
    keep = sorted(set(bc_map) & set(per_ind.index) & set(grm_samples) & set(sv_samples))
    y = np.array([per_ind[s] for s in keep], dtype=float)
    m = ~np.isnan(y)
    keep = [k for i, k in enumerate(keep) if m[i]]
    y = y[m]
    pcm = pc.set_index("sample")
    P = pcm.loc[keep, ["PC%d" % (i + 1) for i in range(10)]].to_numpy(dtype=float)
    b = np.array([1.0 if bc_map.get(s) == "R1041" else 0.0 for s in keep])
    idx = {s: i for i, s in enumerate(sv_samples)}
    x = sv[[idx[s] for s in keep], lead_j].astype(float)
    X = np.column_stack([np.ones(len(y)), P] + ([b[:, None]] if b.std() > 1e-10 else []))
    df = len(y) - X.shape[1] - 1
    yr = rint(y)
    yr = yr - X @ np.linalg.lstsq(X, yr, rcond=None)[0]
    xr = x - X @ np.linalg.lstsq(X, x, rcond=None)[0]
    if xr.std() < 1e-10:
        return np.nan, len(y)
    yr = (yr - yr.mean()) / yr.std()
    xr = (xr - xr.mean()) / xr.std()
    r = float(xr @ yr) / len(yr)
    return r * np.sqrt(df) / np.sqrt(max(1 - r * r, 1e-30)), len(y)


full = calls.groupby(["sample", "hap_local"]).pos_splice.median()
t_full, n_full = t_for(full)
print()
print("GATE  all dips retained: t = %.12f   (published %.12f)  n = %d"
      % (t_full, T_PUBLISHED, n_full))
assert abs(t_full - T_PUBLISHED) < 1e-6, \
    ("the reimplementation does not reproduce the published t (diff %.3g); the jackknife "
     "distribution below would be shifted by the same amount and cannot be read"
     % abs(t_full - T_PUBLISHED))
print("      reproduced to %.2e -- the model below is the published model" % abs(t_full - T_PUBLISHED))

rng = np.random.default_rng(SEED)
key = calls[["sample", "hap_local", "asm_contig", "cdr_group"]].astype(str).agg("|".join, axis=1)
calls = calls.assign(_k=key.to_numpy())
grp_index = {k: np.asarray(v) for k, v in calls.groupby("_k").indices.items()}
pos = calls.pos_splice.to_numpy()
sh = calls[["sample", "hap_local"]].agg(tuple, axis=1).to_numpy()

ts = []
for it in range(NREP):
    pick = np.array([g[rng.integers(len(g))] for g in grp_index.values()])
    s = pd.Series(pos[pick], index=pd.MultiIndex.from_tuples(sh[pick],
                                                             names=["sample", "hap_local"]))
    s = s.groupby(level=[0, 1]).median()
    t, _ = t_for(s)
    ts.append(t)
ts = np.array([t for t in ts if not np.isnan(t)])

print()
print("JACKKNIFE  one random dip kept per merged CDR, %d replicates" % len(ts))
print("   published (all dips)     t = %+.4f" % T_PUBLISHED)
print("   Gao unit (merged span)   t = %+.4f" % T_GAO)
print("   jackknife  mean          t = %+.4f   sd = %.4f" % (ts.mean(), ts.std()))
print("   jackknife  median        t = %+.4f" % np.median(ts))
print("   jackknife  2.5 / 97.5%%   t = %+.4f / %+.4f"
      % (np.percentile(ts, 2.5), np.percentile(ts, 97.5)))
print("   jackknife  min / max     t = %+.4f / %+.4f" % (ts.min(), ts.max()))
print("   replicates with the same SIGN as published : %.1f%%"
      % (100.0 * np.mean(np.sign(ts) == np.sign(T_PUBLISHED))))
print("   replicates with |t| >= 5                    : %.1f%%" % (100.0 * np.mean(np.abs(ts) >= 5)))
print("   replicates with |t| >= 4                    : %.1f%%" % (100.0 * np.mean(np.abs(ts) >= 4)))
np.save(U + "/output/h2h/jackknife_chr8_t.npy", ts)
print()
print("[written] %s/output/h2h/jackknife_chr8_t.npy" % U)
