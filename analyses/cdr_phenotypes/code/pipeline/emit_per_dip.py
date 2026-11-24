#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Emit a PER-DIP table: one row per called CDR, with its position along the array."""
import glob
import importlib.util
import os
import sys
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
UNIT = ROOT + "/analyses/cdr_phenotypes"
MOD = UNIT + "/code/pipeline/compute_phenotypes.py"
PH = UNIT + "/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv"
OUT = UNIT + "/output/phenotypes/per_dip_positions.tsv"

spec = importlib.util.spec_from_file_location("rp2", MOD)
rp2 = importlib.util.module_from_spec(spec)
sys.modules["rp2"] = rp2
spec.loader.exec_module(rp2)
print("imported recompute_phenotypes_v2 (main() not run)")

gmap = rp2.build_globe_map()
print("globe map entries: %d" % len(gmap))

cdrs_by, bins_by = {}, {}
tags = 0
for tagdir in sorted(glob.glob(rp2.CORR + "/*")):
    if not os.path.isdir(tagdir):
        continue
    tags += 1
    for ci, lst in rp2.load_bed_by_contig(tagdir + "/cdrs.bed").items():
        cdrs_by.setdefault(ci, []).extend(lst)
    for ci, lst in rp2.load_bed_by_contig(tagdir + "/binned_freq.bed", 4).items():
        bins_by.setdefault(ci, []).extend(lst)
print("tags %d; contigs with dips %d, with windows %d" % (tags, len(cdrs_by), len(bins_by)))

rows = []
skipped_multicontig = skipped_nodata = 0
for (sample, hap, chrom), by_contig in sorted(gmap.items()):
    if len(by_contig) > 1:
        skipped_multicontig += 1
        continue
    contig, raw = next(iter(by_contig.items()))
    if contig not in bins_by:
        skipped_nodata += 1
        continue
    blocks = rp2.merge(raw)
    aa_s, aa_e = blocks[0][0], blocks[-1][1]
    span, S = aa_e - aa_s, sum(b - a for a, b in blocks)
    if S <= 0:
        continue
    for (cs, ce, _x) in cdrs_by.get(contig, []):
        for a, b in blocks:
            if ce > a and cs < b:
                lo, hi = max(cs, a), min(ce, b)
                m = (lo + hi) / 2.0
                r_span, r_block, r_splice = rp2.frames(m, blocks, aa_s, aa_e, S)
                rows.append(dict(sample=sample, hap_local=hap, chm13_chrom=chrom,
                                 asm_contig=contig, dip_start=int(lo), dip_end=int(hi),
                                 dip_width=int(hi - lo),
                                 pos_span=r_span, pos_block=r_block, pos_splice=r_splice,
                                 active_len_sum=S, active_len_span=span,
                                 n_blocks=len(blocks)))
                break

D = pd.DataFrame(rows)
print("")
print("per-dip rows: %d  over %d haplotype-chromosomes  (%d samples, %d chromosomes)"
      % (len(D), D.groupby(["sample", "hap_local", "chm13_chrom"]).ngroups,
         D["sample"].nunique(), D.chm13_chrom.nunique()))
print("skipped: %d multi-contig, %d without window data" % (skipped_multicontig, skipped_nodata))

ph = pd.read_csv(PH, sep="\t")
agg = (D.groupby(["sample", "hap_local", "chm13_chrom"])
        .agg(n=("pos_splice", "size"), med=("pos_splice", "median")).reset_index())
m = ph.merge(agg, on=["sample", "hap_local", "chm13_chrom"], how="inner")
print("matched haplotype-chromosomes: %d of %d in the phenotype table" % (len(m), len(ph)))

bad_n = m[m.n != m.n_cdrs]
assert bad_n.empty, ("dip COUNT disagrees with n_cdrs for %d haplotype-chromosomes:\n%s"
                     % (len(bad_n), bad_n[["sample", "chm13_chrom", "n", "n_cdrs"]].head()))
print("[check] per-dip counts match n_cdrs for all %d haplotype-chromosomes" % len(m))

d = (m.med - m.cdr_pos_offset_splice).abs()
worst = float(np.nanmax(d))
TOL = 1e-5
assert worst < TOL, ("median of per-dip splice positions disagrees with cdr_pos_offset_splice "
                     "beyond storage precision; max |diff| = %.3g (tolerance %.0e)" % (worst, TOL))
print("[check] median of per-dip positions reproduces cdr_pos_offset_splice to storage precision "
      "(max |diff| = %.2g, table stored at 6 significant figures)" % worst)

D.to_csv(OUT, sep="\t", index=False, float_format="%.6g")
print("")
print("dips per chromosome:")
for c, n in D.chm13_chrom.value_counts().sort_index().items():
    print("   %-7s %6d" % (c, n))
print("[written] %s" % OUT)
