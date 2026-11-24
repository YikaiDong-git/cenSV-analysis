#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Are the BLOCK-only study-wide hits an artefact of gap structure, or is SPLICE losing signal?"""
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = ROOT + "/analyses/cdr_qtl_full_panel"
UNIT = ROOT + "/analyses/cdr_phenotypes"

HITS = [("chr1", "chr1_sv1595_DEL_non_best8_midpoint1000", "cdr_pos_iqr"),
        ("chr4", "chr4_sv912_DEL_mon_best8_fksize100k", "cdr_pos_offset")]

ph = pd.read_csv(UNIT + "/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv", sep="\t")
print("phenotype rows: %d" % len(ph))

print("")
print("=" * 92)
print("T1  does each frame's dip position track gap structure?")
print("=" * 92)
print("%-7s %6s %7s %10s %10s %10s %10s" %
      ("chrom", "n", "frag%", "r_blk_gap", "r_spl_gap", "r_blk_nblk", "r_spl_nblk"))
rows = []
for c, g in ph.groupby("chm13_chrom"):
    if c == "chrY" or len(g) < 40:
        continue
    frag = 100.0 * (g.n_blocks > 1).mean()
    out = dict(chrom=c, n=len(g), frag=frag)
    for tag, col in (("blk", "cdr_pos_offset_block"), ("spl", "cdr_pos_offset_splice")):
        for gtag, gcol in (("gap", "gap_bp"), ("nblk", "n_blocks")):
            d = g[[col, gcol]].dropna()
            if len(d) < 30 or d[gcol].nunique() < 2:
                out["r_%s_%s" % (tag, gtag)] = np.nan
                continue
            out["r_%s_%s" % (tag, gtag)] = stats.spearmanr(d[col], d[gcol]).correlation
    rows.append(out)
    print("%-7s %6d %7.1f %10s %10s %10s %10s" %
          (c, out["n"], frag,
           "%+.3f" % out["r_blk_gap"] if np.isfinite(out.get("r_blk_gap", np.nan)) else "-",
           "%+.3f" % out["r_spl_gap"] if np.isfinite(out.get("r_spl_gap", np.nan)) else "-",
           "%+.3f" % out["r_blk_nblk"] if np.isfinite(out.get("r_blk_nblk", np.nan)) else "-",
           "%+.3f" % out["r_spl_nblk"] if np.isfinite(out.get("r_spl_nblk", np.nan)) else "-"))

T = pd.DataFrame(rows)
frg = T[T.frag > 10]
print("")
print("Across the %d chromosomes with >10%% fragmentation:" % len(frg))
for tag, name in (("blk", "BLOCK"), ("spl", "SPLICE")):
    v = frg["r_%s_gap" % tag].dropna()
    w = frg["r_%s_nblk" % tag].dropna()
    print("  %-6s median |r| vs gap_bp = %.3f   vs n_blocks = %.3f"
          % (name, np.median(np.abs(v)) if len(v) else np.nan,
             np.median(np.abs(w)) if len(w) else np.nan))
vb = np.abs(frg["r_blk_gap"].dropna())
vs = np.abs(frg["r_spl_gap"].dropna())
if len(vb) > 2 and len(vs) > 2:
    st, p = stats.wilcoxon(vb, vs) if len(vb) == len(vs) else stats.mannwhitneyu(vb, vs)
    print("  BLOCK vs SPLICE |r| with gap_bp: P = %.4f" % p)
    if p < 0.05 and np.median(vb) > np.median(vs) * 1.5:
        print("  -> BLOCK tracks gap structure more than SPLICE, as the stability argument "
              "predicts")
    elif p < 0.05 and np.median(vs) > np.median(vb) * 1.5:
        print("  -> SPLICE tracks gap structure MORE than BLOCK -- opposite to expectation")
    else:
        print("  -> INCONCLUSIVE: medians %.3f vs %.3f at P = %.2f. T1 provides NO evidence "
              "either way; the frame choice rests on the a priori stability argument and on "
              "T2, not on this." % (np.median(vb), np.median(vs), p))

print("")
print("=" * 92)
print("T2  does the block-only association survive on SINGLE-BLOCK haplotypes?")
print("=" * 92)
G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
samp = np.array(open(WU + "/geno/full_samples.txt").read().split())
loc2row = {l: i for i, l in enumerate(vt.sv_locus)}


def rint(a):
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def assoc(y, x):
    m = np.isfinite(y) & np.isfinite(x)
    if m.sum() < 25 or np.std(x[m]) == 0:
        return None
    r, p = stats.pearsonr(rint(y[m]), x[m])
    n = int(m.sum())
    t = r * np.sqrt((n - 2) / max(1e-12, 1 - r ** 2))
    return dict(t=float(t), p=float(p), n=n, r=float(r))


for chrom, locus, pheno in HITS:
    print("")
    print("--- %s  %s  %s ---" % (chrom, locus, pheno))
    sub = ph[ph.chm13_chrom == chrom]
    n_single = int((sub.n_blocks == 1).sum())
    print("    haplotypes: %d total, %d single-block (%.1f%%)"
          % (len(sub), n_single, 100.0 * n_single / len(sub)))
    if locus not in loc2row:
        print("    locus absent from the variant table -- skipped")
        continue
    dos = pd.Series(G[loc2row[locus], :], index=samp).replace(-1, np.nan)

    for label, mask in (("ALL haplotypes", sub.index),
                        ("SINGLE-BLOCK only", sub.index[sub.n_blocks == 1])):
        s = sub.loc[mask]
        if not len(s):
            print("    %-20s : no haplotypes" % label)
            continue
        for frame in ("block", "splice"):
            col = "%s_%s" % (pheno, frame) if pheno != "cdr_pos_iqr" else "cdr_pos_iqr_%s" % frame
            if col not in s.columns:
                continue
            agg = s.groupby("sample")[col].mean()
            keep = sorted(set(agg.dropna().index) & set(dos.dropna().index))
            if len(keep) < 25:
                print("    %-20s %-7s: n=%d too small" % (label, frame, len(keep)))
                continue
            res = assoc(agg.loc[keep].values, dos.loc[keep].values)
            if res is None:
                print("    %-20s %-7s: not estimable" % (label, frame))
                continue
            print("    %-20s %-7s: t = %+7.3f   P = %9.2e   n = %d"
                  % (label, frame, res["t"], res["p"], res["n"]))

print("")
print("=" * 92)
print("READING THE RESULT")
print("  T2 collapse under single-block restriction => the block-only hits lived in the gapped")
print("     haplotypes, i.e. they are gap artefacts and SPLICE is right to drop them.")
print("  T2 persistence => SPLICE is losing real signal and the canonical-frame choice must be")
print("     revisited before the paper claims four chromosomes rather than six.")
