#!/usr/bin/env python
import os
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

RD = (CENSV_ROOT + "/discovery"
      "/cdr_phenotypes")
OUT = RD + "/output"
FRAMES = ["span", "block", "splice"]

v2 = pd.read_csv(OUT + "/per_hap_chrom_phenotypes_frames.tsv", sep="\t", engine="python")

PASSTHRU = ["sample", "hap_local", "basecaller", "chm13_chrom", "asm_contig",
            "active_start", "active_end", "n_cdrs", "total_cdr_bp", "median_cdr_width",
            "mean_5mC", "mean_5mC_covwt", "n_windows", "n_intervals", "n_blocks", "gap_bp"]

for frame in FRAMES:
    d = v2[PASSTHRU].copy()
    d["active_len"] = v2["active_len_sum"]
    d["active_len_span"] = v2["active_len_span"]
    d["cdr_fraction"] = v2["cdr_fraction_sum"]
    d["cdrs_per_mb"] = v2["cdrs_per_mb_sum"]
    d["cdr_pos_offset"] = v2["cdr_pos_offset_" + frame]
    d["cdr_pos_iqr"] = v2["cdr_pos_iqr_" + frame]
    d["cdr_pos_var"] = v2["cdr_pos_var_" + frame]

    dd = d.copy()
    dd["cdr_pos_offset_deepest"] = v2["cdr_pos_offset_deepest_" + frame]
    dd["deepest_min5mC"] = v2["deepest_min5mC"]
    dd["deepest_cdr_mid"] = float("nan")

    dest = "%s/compat_%s" % (OUT, frame)
    os.makedirs(dest, exist_ok=True)
    d.to_csv(dest + "/per_hap_chrom_phenotypes.tsv", sep="\t",
             index=False, float_format="%.6g")
    dd.to_csv(dest + "/per_hap_chrom_phenotypes_deepest.tsv", sep="\t",
              index=False, float_format="%.6g")
    print("[written] compat_%s  (%d rows, %d/%d cols, active_len max %.2f Mb)"
          % (frame, len(d), len(d.columns), len(dd.columns), d.active_len.max() / 1e6))

V1D = (CENSV_ROOT + "/discovery"
       "/cdr_qtl/output")
for name in ("corrected", "deepest"):
    want = set(pd.read_csv(V1D + "/per_hap_chrom_phenotypes_%s.tsv" % name, sep="\t",
                           engine="python", nrows=1).columns)
    got = set(pd.read_csv(OUT + "/compat_block/per_hap_chrom_phenotypes_%s.tsv" % name,
                          sep="\t", nrows=1).columns)
    missing = want - got
    assert not missing, "compat _%s is missing v1 columns: %s" % (name, sorted(missing))
    print("[gate] compat _%s carries all %d v1 columns" % (name, len(want)))

a = pd.read_csv(OUT + "/compat_span/per_hap_chrom_phenotypes.tsv", sep="\t")
b = pd.read_csv(OUT + "/compat_block/per_hap_chrom_phenotypes.tsv", sep="\t")
assert (a.active_len <= a.active_len_span + 1).all(), "corrected length exceeds the hull"
sb = a.n_blocks == 1
assert (a.loc[sb, "cdr_pos_offset"] - b.loc[sb, "cdr_pos_offset"]).abs().max() < 1e-9, \
    "frames disagree on single-block arrays -- the shim has reordered rows"
print("[gate] shim invariants PASSED (%d single-block rows checked)" % int(sb.sum()))
