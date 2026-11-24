#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compare the two CDR-unit arms of the chr8 scan."""
import os
import sys
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition")
S = U + "/output/scan/scan_h2h"
LEAD_POS = 44300221
IVAL = (44250000, 44460000)
CHROM = sys.argv[1] if len(sys.argv) > 1 else "chr8"


def load(arm, pheno):
    p = "%s/%s/perchrom/%s_%s_perm.tsv" % (S, arm, CHROM, pheno)
    if not os.path.exists(p):
        return None
    return pd.read_csv(p, sep="\t")


def line(s=""):
    print(s)


a, b = load("dip", "mean_5mC"), load("gao", "mean_5mC")
line("=" * 78)
line("1. NEGATIVE CONTROL  mean_5mC  (window-derived; cannot depend on the call set)")
if a is None or b is None:
    line("   MISSING -- control not run, so nothing below is protected against contamination")
else:
    m = a.merge(b, on=["chrom", "pos", "sv_locus"], suffixes=("_dip", "_gao"))
    d = (m.t_obs_dip - m.t_obs_gao).abs().max()
    line("   %d loci compared; max |t_dip - t_gao| = %.3g" % (len(m), d))
    if d == 0:
        line("   CONTROL PASSED -- the arms differ only where the call set enters")
    else:
        line("   CONTROL FAILED -- the merge leaked into a phenotype it cannot affect")

for pheno in ("cdr_pos_offset", "cdr_pos_offset_deepest"):
    a, b = load("dip", pheno), load("gao", pheno)
    line()
    line("=" * 78)
    line("PHENOTYPE  %s" % pheno)
    if a is None or b is None:
        line("   MISSING (dip=%s gao=%s)" % (a is not None, b is not None))
        continue
    m = a.merge(b, on=["chrom", "pos", "sv_locus"], suffixes=("_dip", "_gao"))
    line("   %d loci in both arms (dip %d, gao %d)" % (len(m), len(a), len(b)))

    line()
    line("   2. the manuscript's lead, chr8:%d DEL 113 bp" % LEAD_POS)
    L = m[m.pos == LEAD_POS]
    if L.empty:
        line("      NOT PRESENT in the merged table")
    else:
        for _, r in L.iterrows():
            line("      %-42s  af=%.3f  n=%d" % (r.sv_locus, r.get("af_dip", np.nan),
                                                 r.get("n_indiv_dip", -1)))
            line("        dip arm : t = %+8.4f   r = %+.4f   p_emp = %.3g"
                 % (r.t_obs_dip, r.r_obs_dip, r.p_emp_perm_dip))
            line("        gao arm : t = %+8.4f   r = %+.4f   p_emp = %.3g"
                 % (r.t_obs_gao, r.r_obs_gao, r.p_emp_perm_gao))
            keep = abs(r.t_obs_gao) / abs(r.t_obs_dip) if r.t_obs_dip else np.nan
            line("        |t| retained under the Gao unit: %.1f%%   sign %s"
                 % (100 * keep,
                    "PRESERVED" if np.sign(r.t_obs_dip) == np.sign(r.t_obs_gao) else "FLIPPED"))

    line()
    line("   3. do the arms agree on which loci rank high?")
    line("      Spearman r(t_dip, t_gao) over all %d loci = %.4f"
         % (len(m), m.t_obs_dip.corr(m.t_obs_gao, method="spearman")))
    line("      Pearson  r(t_dip, t_gao)                 = %.4f"
         % m.t_obs_dip.corr(m.t_obs_gao))
    for k in (10, 50):
        td = set(m.nlargest(k, "abs_t_obs_dip").sv_locus)
        tg = set(m.nlargest(k, "abs_t_obs_gao").sv_locus)
        line("      top-%d overlap: %d/%d" % (k, len(td & tg), k))

    line()
    line("   4. the study-wide interval chr8:44.25-44.46 Mb")
    I = m[(m.pos >= IVAL[0]) & (m.pos <= IVAL[1])]
    line("      %d loci in interval" % len(I))
    line("      max |t| : dip %.4f   gao %.4f"
         % (I.abs_t_obs_dip.max(), I.abs_t_obs_gao.max()))
    line("      loci with |t| > 5 : dip %d   gao %d"
         % (int((I.abs_t_obs_dip > 5).sum()), int((I.abs_t_obs_gao > 5).sum())))

    line()
    line("   top 8 by |t| in the GAO arm:")
    line("      %-44s %9s %9s" % ("sv_locus", "t_gao", "t_dip"))
    for _, r in m.nlargest(8, "abs_t_obs_gao").iterrows():
        line("      %-44s %+9.4f %+9.4f" % (r.sv_locus[:44], r.t_obs_gao, r.t_obs_dip))

    out = "%s/h2h_%s_%s.tsv" % (U + "/output", CHROM, pheno)
    m.to_csv(out, sep="\t", index=False, float_format="%.6g")
    line("   [written] %s" % out)
