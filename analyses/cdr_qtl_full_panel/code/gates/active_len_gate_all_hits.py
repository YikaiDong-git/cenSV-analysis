#!/usr/bin/env python
import os
import struct

import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = ROOT + "/analyses/cdr_qtl_full_panel"
CAN = ROOT + "/analyses/cdr_qtl"
GRM = ROOT + "/subproject_methylation/results/greml/pericen_1m_all"
OUT = WU + "/output/gates/active_len_gate_all_hits.tsv"

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


ph = pd.read_csv(CAN + "/output/per_hap_chrom_phenotypes.tsv", sep="\t")
phd = pd.read_csv(CAN + "/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t")
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], index=gids)
G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
samp = np.array(open(WU + "/geno/full_samples.txt").read().split())
loc2row = {l: i for i, l in enumerate(vt.sv_locus)}

gate = pd.read_csv(WU + "/output/gates/full_hits_gate_1e7.tsv", sep="\t")
key = ["chrom", "sv_locus", "phenotype"]
assert not gate.duplicated(key).any(), "a hit appears twice on %s" % key
print("testing all %d hits on %d chromosomes (the shipped gate tested %d, on 2)"
      % (len(gate), gate.chrom.nunique(), 14))

rows = []
for r in gate.itertuples():
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
    moves = L is not None and L["p"] < 0.05
    rows.append(dict(
        chrom=r.chrom, sv="%s-%d@%d" % (r.svtype, int(r.svlen), int(r.pos)), af=r.af,
        sv_locus=r.sv_locus, phenotype=r.phenotype, p_scan=r.p_best,
        t_active_len=L["t"], p_active_len=L["p"], n=L["n"],
        t_active_start=S["t"], p_active_start=S["p"],
        t_active_end=E["t"], p_active_end=E["p"],
        p_geometry=min(L["p"], S["p"], E["p"]),
        geometry_moves=min(L["p"], S["p"], E["p"]) < 0.05,
        which_moves=",".join(k for k, v in (("len", L["p"]), ("start", S["p"]),
                                            ("end", E["p"])) if v < 0.05) or "none",
        length_normalised=lennorm,
        pass_headline=bool(r.pass_headline),
        prior_verdict=r.overall, boundary_verdict=r.boundary_verdict,
        ancestry_verdict=r.ancestry_verdict,
        length_verdict=("LENGTH_MOVES" if (lennorm and moves) else
                        ("length_stable" if lennorm else
                         ("LENGTH_MOVES_not_normalised" if moves else "n/a not normalised")))))

out = pd.DataFrame(rows)

old = pd.read_csv(WU + "/output/gates/active_len_gate.tsv", sep="\t")
mk = ["chrom", "sv", "phenotype", "af"]
old["af"] = old.af.round(6)
out["af6"] = out.af.round(6)
m = old.merge(out.rename(columns={"af6": "af_k"}), left_on=mk,
              right_on=["chrom", "sv", "phenotype", "af_k"], suffixes=("_old", "_new"))
assert len(m) == len(old), ("matched %d of the shipped gate's %d rows; the new run does not "
                            "cover what the old one did" % (len(m), len(old)))
worst = float((m.t_active_len_old - m.t_active_len_new).abs().max())
assert worst < 5e-3, ("re-running the shipped gate's own rows gives a different t (worst "
                      "%.4g): the estimator here is not the one whose numbers are in the "
                      "manuscript, so its new rows cannot be trusted either" % worst)
print("reproduction: %d shipped rows re-derived, worst |dt| = %.2e" % (len(m), worst))
ALPHA = 0.05 / (11129 * 8)
assert ((out.p_scan <= ALPHA) == out.pass_headline).all(), \
    "pass_headline disagrees with p_scan <= alpha, so the flag cannot be used as the definition"
out["tier"] = np.where(out.pass_headline, "significant", "suggestive")
assert (out.p_scan <= 2 * ALPHA).all(), \
    ("%d leads sit beyond 2*alpha, so 'suggestive' would be a name for a ragged remainder "
     "rather than a tier" % int((out.p_scan > 2 * ALPHA).sum()))

out = out.drop(columns=["af6"])

out.to_csv(OUT, sep="\t", index=False, float_format="%.6g")
print("\nwrote %s" % OUT)
_written = pd.read_csv(OUT, sep="\t")
assert "tier" in _written.columns, \
    ("tier was computed but is not in %s -- something below the write is assigning columns "
     "that never reach the file" % OUT)
assert set(_written.tier) == {"significant", "suggestive"}, sorted(set(_written.tier))
print("  verified: the written table carries tier (%d significant, %d suggestive)"
      % (int((_written.tier == "significant").sum()), int((_written.tier == "suggestive").sum())))

print("\n%-6s %-22s %-24s %8s %11s  %s"
      % ("chrom", "SV", "phenotype", "t(len)", "p(len)", "length verdict"))
for r in out.sort_values("t_active_len").itertuples():
    print("%-6s %-22s %-24s %+8.2f %11.2e  %s"
          % (r.chrom, r.sv, r.phenotype, r.t_active_len, r.p_active_len, r.length_verdict))

print("\n=== by chromosome: does the lead's dosage predict array length? ===")
g = out.assign(moves=out.p_active_len < 0.05).groupby("chrom").agg(
    hits=("sv", "size"), moves=("moves", "sum"), worst_t=("t_active_len",
                                                          lambda s: s.abs().max()))
print(g.to_string())
print("\nhits whose lead predicts array length at p<0.05: %d of %d, on %d chromosomes"
      % (int((out.p_active_len < 0.05).sum()), len(out),
         out[out.p_active_len < 0.05].chrom.nunique()))



hl = out[out.pass_headline]
print("\ntiers: significant p <= %.3g | suggestive %.3g < p <= %.3g (Bonferroni 0.05 and 0.10)"
      % (ALPHA, ALPHA, 2 * ALPHA))
print("study-wide significant: %d of %d hits, %d chromosomes"
      % (len(hl), len(out), hl.chrom.nunique()))
print("  lead predicts array length: %d of %d, on %d chromosomes"
      % (int((hl.p_active_len < 0.05).sum()), len(hl),
         hl[hl.p_active_len < 0.05].chrom.nunique()))
print(hl.assign(moves=hl.p_active_len < 0.05).groupby("chrom")
        .agg(hits=("sv", "size"), moves=("moves", "sum")).to_string())
print("\n=== array GEOMETRY, not only length ===")
print("A length-normalised phenotype is built from active_start and active_len, so a lead whose")
print("dosage predicts either one is confounded. Counted per coordinate, at P < 0.05 each:")
for nm, col in (("length", "p_active_len"), ("start", "p_active_start"),
                ("end", "p_active_end")):
    print("  %-7s %2d of %d leads" % (nm, int((out[col] < 0.05).sum()), len(out)))
print("  ANY     %2d of %d leads, on %d chromosomes"
      % (int(out.geometry_moves.sum()), len(out), out[out.geometry_moves].chrom.nunique()))
print("\n  moves on start or end but NOT length: %d"
      % int((out.geometry_moves & (out.p_active_len >= 0.05)).sum()))
for r in out[out.geometry_moves & (out.p_active_len >= 0.05)].itertuples():
    print("    %-6s %-24s %-24s %-11s which=%s"
          % (r.chrom, r.sv, r.phenotype, r.tier, r.which_moves))
print()
for nm, sub in (("significant", out[out.pass_headline]), ("suggestive", out[~out.pass_headline])):
    print("%s: length %d of %d | geometry %d of %d, on %d chromosomes"
          % (nm, int((sub.p_active_len < 0.05).sum()), len(sub),
             int(sub.geometry_moves.sum()), len(sub),
             sub[sub.geometry_moves].chrom.nunique()))
    print(sub.groupby("chrom").agg(leads=("sv", "size"),
                                   length=("p_active_len", lambda s: int((s < 0.05).sum())),
                                   geometry=("geometry_moves", "sum")).to_string())
    print()

sg = out[~out.pass_headline]
print("\nsuggestive: %d of %d hits, %d chromosomes" % (len(sg), len(out), sg.chrom.nunique()))
print("  lead predicts array length: %d of %d, on %d chromosomes"
      % (int((sg.p_active_len < 0.05).sum()), len(sg),
         sg[sg.p_active_len < 0.05].chrom.nunique()))
print(sg.assign(moves=sg.p_active_len < 0.05).groupby("chrom")
        .agg(hits=("sv", "size"), moves=("moves", "sum")).to_string())
