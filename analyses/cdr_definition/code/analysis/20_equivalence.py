#!/usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition")
S = U + "/output/scan/scan_gw_1e5"
POS = ("cdr_pos_offset", "cdr_pos_offset_deepest")

rows = []
for p in sorted(glob.glob(S + "/dip/perchrom/*_perm.tsv")):
    base = os.path.basename(p)
    q = S + "/gao/perchrom/" + base
    if not os.path.exists(q):
        continue
    chrom = base.split("_")[0]
    pheno = base[len(chrom) + 1:-len("_perm.tsv")]
    a = pd.read_csv(p, sep="\t")
    b = pd.read_csv(q, sep="\t")
    m = a.merge(b, on=["chrom", "pos", "sv_locus"], suffixes=("_d", "_g"))
    m["phenotype"] = pheno
    rows.append(m[["chrom", "pos", "sv_locus", "phenotype", "t_obs_d", "t_obs_g",
                   "abs_t_obs_d", "abs_t_obs_g"]])
M = pd.concat(rows, ignore_index=True)
M["d"] = M.t_obs_g - M.t_obs_d
M["ad"] = M.d.abs()

print("loci-by-phenotype cells compared: %d" % len(M))
print()
print("EQUIVALENCE MARGIN on t, position phenotypes only (%s):" % ", ".join(POS))
P = M[M.phenotype.isin(POS)]
print("   n = %d" % len(P))
for q in (0.5, 0.9, 0.99, 0.999, 1.0):
    print("   %6.1f%% of loci move by <= %.4f" % (100 * q, P.ad.quantile(q)))
print("   mean signed shift  %+.5f   (a systematic shrink would show as a shift toward zero)"
      % P.d.mean())
print("   mean |t| dip %.4f   gao %.4f   ratio %.4f"
      % (P.abs_t_obs_d.mean(), P.abs_t_obs_g.mean(),
         P.abs_t_obs_g.mean() / P.abs_t_obs_d.mean()))

print()
print("Restricted to loci that could matter (|t| >= 4 in either arm):")
B = P[(P.abs_t_obs_d >= 4) | (P.abs_t_obs_g >= 4)]
print("   n = %d" % len(B))
if len(B):
    print("   max |t_gao - t_dip| = %.4f   (worst case, %s %s at %d)"
          % (B.ad.max(), B.loc[B.ad.idxmax(), "chrom"], B.loc[B.ad.idxmax(), "phenotype"],
             B.loc[B.ad.idxmax(), "pos"]))
    print("   median |t_gao - t_dip| = %.4f" % B.ad.median())
    print("   |t| retained, median %.1f%%   min %.1f%%"
          % (100 * (B.abs_t_obs_g / B.abs_t_obs_d).median(),
             100 * (B.abs_t_obs_g / B.abs_t_obs_d).min()))
    print("   sign preserved on %d of %d (%.1f%%)"
          % (int((np.sign(B.t_obs_d) == np.sign(B.t_obs_g)).sum()), len(B),
             100 * (np.sign(B.t_obs_d) == np.sign(B.t_obs_g)).mean()))
    print()
    print("   the 8 loci that move most:")
    print("   %-8s %-24s %10s %9s %9s %8s" %
          ("chrom", "phenotype", "pos", "t_dip", "t_gao", "|delta|"))
    for _, r in B.nlargest(8, "ad").iterrows():
        print("   %-8s %-24s %10d %+9.4f %+9.4f %8.4f"
              % (r.chrom, r.phenotype, r.pos, r.t_obs_d, r.t_obs_g, r.ad))

print()
print("For contrast, the phenotypes the unit is NOT expected to preserve:")
print("   %-24s %8s %12s %12s" % ("phenotype", "n", "median |d|", "max |d|"))
for ph, g in M[~M.phenotype.isin(POS)].groupby("phenotype"):
    print("   %-24s %8d %12.4f %12.4f" % (ph, len(g), g.ad.median(), g.ad.max()))

M.to_csv(U + "/output/gw/equivalence_margins.tsv", sep="\t", index=False, float_format="%.6g")
print()
print("[written] %s/output/gw/equivalence_margins.tsv" % U)
