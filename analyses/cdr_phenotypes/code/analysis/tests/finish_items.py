#!/usr/bin/env python
"""Two loose ends for the manuscript diff."""
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = ROOT + "/analyses/cdr_qtl_full_panel"
UNIT = ROOT + "/analyses/cdr_phenotypes"

vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t",
                 usecols=["sv_locus", "chrom", "pos", "end", "svtype", "svlen", "af"])
R = pd.read_csv(UNIT + "/output/per_frame/gate_rederive/len_gate_alleads/active_len_gate_alleads.tsv",
                sep="\t")
R = R.merge(vt, on="sv_locus", how="left", suffixes=("", "_v"))

print("=== chr9 gate leads, by strength of the array-length effect ===")
c9 = R[R.chrom == "chr9"].sort_values("t_sum")
for r in c9.itertuples():
    print("  %-40s %s-%d@%d  AF %.3f   t_SUM %+.3f  P %.3g  (SPAN %+.3f)"
          % (r.sv_locus[:40], r.svtype, r.svlen, r.pos, r.af, r.t_sum, r.p_sum, r.t_span))

print("")
print("=== chr8 gate leads ===")
c8 = R[R.chrom == "chr8"].sort_values("t_sum")
for r in c8.itertuples():
    print("  %-40s %s-%d@%d  AF %.3f   %-24s t_SUM %+.3f  P %.3g"
          % (r.sv_locus[:40], r.svtype, r.svlen, r.pos, r.af, r.phenotype,
             r.t_sum, r.p_sum))

print("")
print("=== chr5 gate lead ===")
for r in R[R.chrom == "chr5"].itertuples():
    print("  %-40s %s-%d@%d  AF %.3f   t_SUM %+.3f  P %.3g  (SPAN %+.3f, P %.3g)"
          % (r.sv_locus[:40], r.svtype, r.svlen, r.pos, r.af,
             r.t_sum, r.p_sum, r.t_span, r.p_span))

print("")
print("=" * 90)
print("ITEM 26  array length vs assembly coverage")
print("=" * 90)
ph = pd.read_csv(UNIT + "/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv", sep="\t")
cov_src = WU + "/output/burden/individual_burden_and_quality.tsv"
cov = pd.read_csv(cov_src, sep="\t")
covcol = [c for c in cov.columns if "cov" in c.lower()]
print("coverage table: %s" % cov_src)
print("  columns: %s" % ", ".join(cov.columns))
assert covcol, "no coverage column found"
cc = covcol[0]
print("  using coverage column: %s" % cc)

per_sample = ph.groupby("sample").agg(
    len_sum=("active_len_sum", "mean"), len_span=("active_len_span", "mean")).reset_index()
m = per_sample.merge(cov[["sample", cc]], on="sample", how="inner").dropna()
print("  samples with both length and coverage: %d" % len(m))

for label, col in (("CORRECTED (merged intervals)", "len_sum"),
                   ("SHIPPED (convex hull)", "len_span")):
    r, p = stats.pearsonr(m[col], m[cc])
    q = pd.qcut(m[cc], 5, labels=False, duplicates="drop")
    med = m.groupby(q)[col].median() / 1e6
    print("")
    print("  %s" % label)
    print("     r(length, coverage) = %+.3f   P = %.3g   n = %d" % (r, p, len(m)))
    print("     median length by coverage quintile (Mb): %s"
          % ", ".join("%.2f" % v for v in med.values))
