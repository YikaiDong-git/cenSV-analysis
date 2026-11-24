#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Per-CALL position tables for both arms, and the Figure 5C dosage medians recomputed."""
import glob
import importlib.util
import os
import sys
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
U = ROOT + "/analyses/cdr_definition"
WU = ROOT + "/analyses/cdr_qtl_full_panel"
MOD = U + "/code/pipeline/03_recompute_phenotypes.py"
LEAD = 44300221
SHIPPED_MEDIANS = [0.82, 0.71, 0.55]

spec = importlib.util.spec_from_file_location("rp", MOD)
rp = importlib.util.module_from_spec(spec)
sys.modules["rp"] = rp
spec.loader.exec_module(rp)

ARMS = {"dip": rp.CORR, "gao": U + "/output/callsets/gao_cdrs"}
PHENO = {"dip": U + "/output/phenotypes/per_hap_chrom_phenotypes_reproduced.tsv",
         "gao": U + "/output/phenotypes/per_hap_chrom_phenotypes_gao.tsv"}

gmap = rp.build_globe_map()

bins_by = {}
for tagdir in sorted(glob.glob(rp.CORR + "/*")):
    if not os.path.isdir(tagdir):
        continue
    for ci, lst in rp.load_bed_by_contig(tagdir + "/binned_freq.bed", 4).items():
        bins_by.setdefault(ci, []).extend(lst)

tables = {}
for arm, root in ARMS.items():
    calls_by = {}
    for tagdir in sorted(glob.glob(rp.CORR + "/*")):
        if not os.path.isdir(tagdir):
            continue
        bed = os.path.join(root, os.path.basename(tagdir), "cdrs.bed")
        for ci, lst in rp.load_bed_by_contig(bed).items():
            calls_by.setdefault(ci, []).extend(lst)

    rows = []
    for (sample, hap, chrom), by_contig in sorted(gmap.items()):
        if len(by_contig) > 1:
            continue
        contig, raw = next(iter(by_contig.items()))
        if contig not in bins_by:
            continue
        blocks = rp.merge(raw)
        aa_s, aa_e = blocks[0][0], blocks[-1][1]
        span, S = aa_e - aa_s, sum(b - a for a, b in blocks)
        if S <= 0:
            continue
        for (cs, ce, _x) in calls_by.get(contig, []):
            for a, b in blocks:
                if ce > a and cs < b:
                    lo, hi = max(cs, a), min(ce, b)
                    m = (lo + hi) / 2.0
                    r_span, r_block, r_splice = rp.frames(m, blocks, aa_s, aa_e, S)
                    rows.append(dict(sample=sample, hap_local=hap, chm13_chrom=chrom,
                                     asm_contig=contig, call_start=int(lo), call_end=int(hi),
                                     call_width=int(hi - lo), pos_span=r_span,
                                     pos_block=r_block, pos_splice=r_splice,
                                     active_len_sum=S, active_len_span=span,
                                     n_blocks=len(blocks)))
                    break
    D = pd.DataFrame(rows)

    ph = pd.read_csv(PHENO[arm], sep="\t")
    agg = (D.groupby(["sample", "hap_local", "chm13_chrom"])
            .agg(n=("pos_splice", "size"), med=("pos_splice", "median")).reset_index())
    m = ph.merge(agg, on=["sample", "hap_local", "chm13_chrom"], how="inner")
    bad = m[m.n != m.n_cdrs]
    assert bad.empty, "%s: per-call count disagrees with n_cdrs on %d rows" % (arm, len(bad))
    worst = float(np.nanmax((m.med - m.cdr_pos_offset_splice).abs()))
    assert worst < 1e-5, \
        "%s: median of per-call positions disagrees with cdr_pos_offset_splice (max %.3g)" \
        % (arm, worst)

    out = "%s/output/phenotypes/per_call_positions_%s.tsv" % (U, arm)
    D.to_csv(out, sep="\t", index=False, float_format="%.6g")
    tables[arm] = D
    print("[%s] %d calls over %d haplotype-chromosomes; both checks passed -> %s"
          % (arm, len(D), m.shape[0], os.path.basename(out)))

vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
loc = vt[(vt.chrom == "chr8") & (vt.pos == LEAD)].sv_locus.iloc[0]
G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
sm = open(WU + "/geno/full_samples.txt").read().split()
k = {l: j for j, l in enumerate(vt.sv_locus)}[loc]
dose = pd.Series(G[k], index=sm, name="dose")
print()
print("lead locus %s   dosage classes: %s"
      % (loc, sorted(pd.Series(dose).dropna().unique())))

res = {}
for arm, D in tables.items():
    d = D[D.chm13_chrom == "chr8"].copy()
    d["rel"] = (d.pos_splice + 1.0) / 2.0
    assert d.rel.between(-1e-9, 1 + 1e-9).all(), "%s: relative position outside [0, 1]" % arm
    d = d.join(dose, on="sample")
    d = d[d.dose.notna()]
    res[arm] = d

DOSE = sorted(res["dip"].dose.unique())
print()
print("Figure 5C -- median relative position of a call, by lead-deletion dosage")
print("%-6s %10s %12s %10s %12s %10s" %
      ("dose", "n_dip", "median_dip", "n_gao", "median_gao", "shipped"))
print("-" * 66)
med_dip, med_gao = [], []
for i, dv in enumerate(DOSE):
    a = res["dip"][res["dip"].dose == dv]
    b = res["gao"][res["gao"].dose == dv]
    ma, mb = float(a.rel.median()), float(b.rel.median())
    med_dip.append(ma)
    med_gao.append(mb)
    sh = SHIPPED_MEDIANS[i] if i < len(SHIPPED_MEDIANS) else float("nan")
    print("%-6.1f %10d %12.4f %10d %12.4f %10.2f" % (dv, len(a), ma, len(b), mb, sh))

print()
d = max(abs(round(x, 2) - s) for x, s in zip(med_dip, SHIPPED_MEDIANS))
assert d < 1e-9, ("the dip arm does not reproduce the published 0.82 / 0.71 / 0.55; got %s. "
                  "Nothing in the Gao column can be interpreted until this matches."
                  % [round(x, 2) for x in med_dip])
print("[gate] the dip arm reproduces the published medians %s exactly"
      % [round(x, 2) for x in med_dip])
print("       under the Gao CDR unit they are %s" % [round(x, 2) for x in med_gao])
sp_dip = med_dip[0] - med_dip[-1]
sp_gao = med_gao[0] - med_gao[-1]
print("       monotone decrease across classes: dip %+.3f   gao %+.3f   (%.0f%% retained)"
      % (sp_dip, sp_gao, 100.0 * sp_gao / sp_dip if sp_dip else float("nan")))

from scipy.stats import kruskal
for arm in ("dip", "gao"):
    g = [res[arm][res[arm].dose == dv].rel.dropna() for dv in DOSE]
    h, p = kruskal(*g)
    print("       Kruskal-Wallis over dosage classes, %s arm: H = %.2f, P = %.3g" % (arm, h, p))
