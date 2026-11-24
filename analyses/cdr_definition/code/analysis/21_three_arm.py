#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""chr8 under all three call sets, side by side."""
import os
import pandas as pd
import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition")
S = U + "/output/scan/scan_h2h"
ARMS = ["dip", "gao", "extgao"]
LEAD = 44300221


def load(arm, pheno):
    p = "%s/%s/perchrom/chr8_%s_perm.tsv" % (S, arm, pheno)
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


print("NEGATIVE CONTROL  mean_5mC")
ref = load("dip", "mean_5mC")
for arm in ARMS[1:]:
    o = load(arm, "mean_5mC")
    if o is None:
        print("   %-8s MISSING" % arm)
        continue
    m = ref.merge(o, on=["chrom", "pos", "sv_locus"], suffixes=("_a", "_b"))
    d = (m.t_obs_a - m.t_obs_b).abs().max()
    print("   %-8s vs dip: %d loci, max |dt| = %.3g   %s"
          % (arm, len(m), d, "PASS" if d == 0 else "FAIL"))

for pheno in ("cdr_pos_offset", "cdr_pos_offset_deepest"):
    print()
    print("=" * 74)
    print("PHENOTYPE  %s" % pheno)
    tabs = {a: load(a, pheno) for a in ARMS}
    if any(v is None for v in tabs.values()):
        print("   missing: %s" % [a for a, v in tabs.items() if v is None])
        continue
    base = tabs["dip"][["chrom", "pos", "sv_locus", "t_obs", "abs_t_obs"]].rename(
        columns={"t_obs": "t_dip", "abs_t_obs": "at_dip"})
    for a in ARMS[1:]:
        base = base.merge(tabs[a][["chrom", "pos", "sv_locus", "t_obs", "abs_t_obs"]].rename(
            columns={"t_obs": "t_" + a, "abs_t_obs": "at_" + a}),
            on=["chrom", "pos", "sv_locus"])
    print("   %d loci in all three arms" % len(base))

    L = base[base.pos == LEAD].drop_duplicates("sv_locus")
    print()
    print("   lead 113-bp DEL at chr8:%d" % LEAD)
    for _, r in L.iterrows():
        print("      %-44s dip %+8.4f   gao %+8.4f   extgao %+8.4f"
              % (r.sv_locus[:44], r.t_dip, r.t_gao, r.t_extgao))

    print()
    print("   agreement with the dip arm over all %d loci:" % len(base))
    for a in ARMS[1:]:
        sp = base.t_dip.corr(base["t_" + a], method="spearman")
        pe = base.t_dip.corr(base["t_" + a])
        mx = (base.t_dip - base["t_" + a]).abs().max()
        top = len(set(base.nlargest(50, "at_dip").sv_locus)
                  & set(base.nlargest(50, "at_" + a).sv_locus))
        print("      %-8s Spearman %.4f   Pearson %.4f   max |dt| %.4f   top-50 shared %d/50"
              % (a, sp, pe, mx, top))

    print()
    print("   gao vs extgao (does correcting extent ON TOP of merging change anything?)")
    sp = base.t_gao.corr(base.t_extgao, method="spearman")
    mx = (base.t_gao - base.t_extgao).abs().max()
    print("      Spearman %.4f   max |dt| %.4f" % (sp, mx))

    out = "%s/output/h2h/three_arm_chr8_%s.tsv" % (U, pheno)
    base.to_csv(out, sep="\t", index=False, float_format="%.6g")
    print("   [written] %s" % os.path.basename(out))
