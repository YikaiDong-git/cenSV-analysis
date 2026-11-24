#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition")
OUT = U + "/output"
V1D = (CENSV_ROOT + "/analyses/"
       "cdr_qtl/output")
FRAME = "splice"
PASSTHRU = ["sample", "hap_local", "basecaller", "chm13_chrom", "asm_contig",
            "active_start", "active_end", "n_cdrs", "total_cdr_bp", "median_cdr_width",
            "mean_5mC", "mean_5mC_covwt", "n_windows", "n_intervals", "n_blocks", "gap_bp"]

v2 = pd.read_csv(OUT + "/per_hap_chrom_phenotypes_extgao.tsv", sep="\t", engine="python")
d = v2[PASSTHRU].copy()
d["active_len"] = v2["active_len_sum"]
d["active_len_span"] = v2["active_len_span"]
d["cdr_fraction"] = v2["cdr_fraction_sum"]
d["cdrs_per_mb"] = v2["cdrs_per_mb_sum"]
d["cdr_pos_offset"] = v2["cdr_pos_offset_" + FRAME]
d["cdr_pos_iqr"] = v2["cdr_pos_iqr_" + FRAME]
d["cdr_pos_var"] = v2["cdr_pos_var_" + FRAME]
dd = d.copy()
dd["cdr_pos_offset_deepest"] = v2["cdr_pos_offset_deepest_" + FRAME]
dd["deepest_min5mC"] = v2["deepest_min5mC"]
dd["deepest_cdr_mid"] = float("nan")

dest = "%s/compat/compat_extgao_%s" % (OUT, FRAME)
os.makedirs(dest, exist_ok=True)
d.to_csv(dest + "/per_hap_chrom_phenotypes.tsv", sep="\t", index=False,
         float_format="%.6g")
dd.to_csv(dest + "/per_hap_chrom_phenotypes_deepest.tsv", sep="\t", index=False,
          float_format="%.6g")
print("[written] compat_extgao_%s  %d rows" % (FRAME, len(d)))

for name in ("corrected", "deepest"):
    want = set(pd.read_csv(V1D + "/per_hap_chrom_phenotypes_%s.tsv" % name, sep="\t",
                           engine="python", nrows=1).columns)
    got = set(pd.read_csv("%s/per_hap_chrom_phenotypes_%s.tsv" % (dest, name),
                          sep="\t", nrows=1).columns)
    assert not want - got, "missing v1 columns for _%s: %s" % (name, sorted(want - got))
    print("[gate] all %d v1 columns present for _%s" % (len(want), name))

ref = pd.read_csv(OUT + "/compat/compat_dip_splice/per_hap_chrom_phenotypes.tsv", sep="\t")
KEY = ["sample", "hap_local", "chm13_chrom"]
m = ref[KEY + ["mean_5mC", "active_len", "n_windows", "cdr_pos_offset"]].merge(
    d[KEY + ["mean_5mC", "active_len", "n_windows", "cdr_pos_offset"]], on=KEY,
    suffixes=("_dip", "_ext"))
assert len(m) == len(ref) == len(d), "row sets differ (%d/%d/%d)" % (len(m), len(ref), len(d))
for c in ("mean_5mC", "active_len", "n_windows"):
    mx = (m[c + "_dip"] - m[c + "_ext"]).abs().max()
    assert mx == 0, "%s differs by up to %g -- extension leaked into a window quantity" % (c, mx)
print("[gate] mean_5mC, active_len, n_windows identical to the dip arm (%d rows)" % len(m))

ok = m.cdr_pos_offset_dip.notna() & m.cdr_pos_offset_ext.notna()
print()
print("cdr_pos_offset, extent+merge arm vs dip arm: r = %.4f over %d rows"
      % (m.loc[ok, "cdr_pos_offset_dip"].corr(m.loc[ok, "cdr_pos_offset_ext"]), int(ok.sum())))
nd = int((m.cdr_pos_offset_dip != m.cdr_pos_offset_ext).sum())
assert nd > 0.5 * len(m), \
    "the extgao arm is not distinct from the dip arm (%d of %d rows differ)" % (nd, len(m))
print("differs on %d of %d rows (%.1f%%) -- the arm is distinct data" % (nd, len(m),
                                                                        100.0 * nd / len(m)))
