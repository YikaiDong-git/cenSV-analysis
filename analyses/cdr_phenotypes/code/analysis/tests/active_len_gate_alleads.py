#!/usr/bin/env python
"""Array-length test for EVERY gate-considered lead, under both length definitions."""
import os
import struct
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = ROOT + "/analyses/cdr_qtl_full_panel"
UNIT = ROOT + "/analyses/cdr_phenotypes"
GRM = ROOT + "/subproject_methylation/results/greml/pericen_1m_all"
V2 = UNIT + "/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv"
OUT = UNIT + "/output/per_frame/gate_rederive/len_gate_alleads"


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


def rint(a):
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


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
    return dict(t=float(t), p=float(2 * (1 - stats.t.cdf(abs(t), len(y) - X.shape[1]))),
                n=int(len(y)))


os.makedirs(OUT, exist_ok=True)

ph = pd.read_csv(V2, sep="\t")
assert "active_len_span" in ph.columns and "active_len_sum" in ph.columns, ph.columns.tolist()
print("corrected phenotype table: %d rows, %d samples, %d chromosomes"
      % (len(ph), ph["sample"].nunique(), ph.chm13_chrom.nunique()))

K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], index=gids)

G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
samp = np.array(open(WU + "/geno/full_samples.txt").read().split())
loc2row = {l: i for i, l in enumerate(vt.sv_locus)}

gate = pd.read_csv(WU + "/output/gates/full_hits_gate_1e7.tsv", sep="\t")
print("gate rows considered: %d  (SURVIVES %d, other %d)"
      % (len(gate), int((gate.overall == "SURVIVES").sum()),
         int((gate.overall != "SURVIVES").sum())))
print("verdict tally: %s" % gate.overall.value_counts().to_dict())

rows = []
for r in gate.itertuples():
    if r.sv_locus not in loc2row:
        print("  skip %s: absent from the variant table" % r.sv_locus)
        continue
    s = ph[ph.chm13_chrom == r.chrom]
    agg = s.groupby("sample").agg(
        span=("active_len_span", "mean"), summ=("active_len_sum", "mean"))
    bc = s.groupby("sample")["basecaller"].agg(
        lambda v: 1.0 if v.astype(str).str.contains("1041").any() else 0.0)
    agg = agg.join(bc.rename("bc"))
    dos = pd.Series(G[loc2row[r.sv_locus], :], index=samp).replace(-1, np.nan)
    keep = sorted(set(agg.dropna().index) & set(dos.dropna().index) & set(PCs.index))
    if len(keep) < 30:
        continue
    d = dos.loc[keep].values
    COV = np.column_stack([PCs.loc[keep].values, agg.loc[keep, "bc"].values[:, None]])
    a = slope(rint(agg.loc[keep, "span"].values), d, COV)
    b = slope(rint(agg.loc[keep, "summ"].values), d, COV)
    if a is None or b is None:
        continue
    rows.append(dict(chrom=r.chrom, sv_locus=r.sv_locus, phenotype=r.phenotype,
                     overall=r.overall, af=getattr(r, "af", np.nan),
                     t_span=a["t"], p_span=a["p"], n_span=a["n"],
                     t_sum=b["t"], p_sum=b["p"], n_sum=b["n"],
                     dt=b["t"] - a["t"]))

R = pd.DataFrame(rows).drop_duplicates(subset=["chrom", "sv_locus", "phenotype"])
R = R.sort_values(["chrom", "t_sum"])
R.to_csv(OUT + "/active_len_gate_alleads.tsv", sep="\t", index=False, float_format="%.6g")

print("")
print("%-6s %-42s %-24s %-22s %9s %10s %9s %10s %8s"
      % ("chrom", "sv_locus", "phenotype", "gate verdict",
         "t_SPAN", "P_SPAN", "t_SUM", "P_SUM", "dt"))
for r in R.itertuples():
    print("%-6s %-42s %-24s %-22s %9.3f %10.2e %9.3f %10.2e %+8.3f"
          % (r.chrom, r.sv_locus[:42], r.phenotype[:24], r.overall,
             r.t_span, r.p_span, r.t_sum, r.p_sum, r.dt))

print("")
print("=== per-chromosome strongest |t| under the CORRECTED (sum) length ===")
for c, g in R.groupby("chrom"):
    i = g.t_sum.abs().idxmax()
    print("  %-6s t_SUM = %+.3f  P = %.3g  n = %d   (SPAN was t = %+.3f, P = %.3g)"
          % (c, g.loc[i, "t_sum"], g.loc[i, "p_sum"], g.loc[i, "n_sum"],
             g.loc[i, "t_span"], g.loc[i, "p_span"]))
print("")
print("[written] " + OUT + "/active_len_gate_alleads.tsv")
