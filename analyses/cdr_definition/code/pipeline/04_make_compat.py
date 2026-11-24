#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

U = (CENSV_ROOT + "/analyses/"
     "cdr_definition")
OUT = U + "/output"
V1D = (CENSV_ROOT + "/analyses/"
       "cdr_qtl/output")
SHIPPED = (CENSV_ROOT + "/analyses/"
           "cdr_phenotypes/output/compat/compat_splice")

FRAME = "splice"
PASSTHRU = ["sample", "hap_local", "basecaller", "chm13_chrom", "asm_contig",
            "active_start", "active_end", "n_cdrs", "total_cdr_bp", "median_cdr_width",
            "mean_5mC", "mean_5mC_covwt", "n_windows", "n_intervals", "n_blocks", "gap_bp"]

ARMS = {"dip": "per_hap_chrom_phenotypes_reproduced.tsv",
        "gao": "per_hap_chrom_phenotypes_gao.tsv"}

built = {}
for arm, fn in ARMS.items():
    src = os.path.join(OUT, fn)
    assert os.path.exists(src), "missing phenotype table for arm %s: %s" % (arm, src)
    v2 = pd.read_csv(src, sep="\t", engine="python")

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

    dest = "%s/compat/compat_%s_%s" % (OUT, arm, FRAME)
    os.makedirs(dest, exist_ok=True)
    d.to_csv(dest + "/per_hap_chrom_phenotypes.tsv", sep="\t",
             index=False, float_format="%.6g")
    dd.to_csv(dest + "/per_hap_chrom_phenotypes_deepest.tsv", sep="\t",
              index=False, float_format="%.6g")
    built[arm] = d
    print("[written] compat_%s_%s  %d rows" % (arm, FRAME, len(d)))

for name in ("corrected", "deepest"):
    a = os.path.join(SHIPPED, "per_hap_chrom_phenotypes_%s.tsv" % name)
    b = "%s/compat/compat_dip_%s/per_hap_chrom_phenotypes_%s.tsv" % (OUT, FRAME, name)
    assert os.path.exists(a), "shipped shim not found: %s" % a
    if open(a, "rb").read() != open(b, "rb").read():
        A = pd.read_csv(a, sep="\t", engine="python")
        B = pd.read_csv(b, sep="\t", engine="python")
        print("SHAPE shipped=%s mine=%s" % (A.shape, B.shape))
        print("cols only in shipped:", sorted(set(A.columns) - set(B.columns)))
        print("cols only in mine   :", sorted(set(B.columns) - set(A.columns)))
        raise SystemExit("GATE FAIL: the dip arm does not reproduce the shipped compat_%s" % name)
    print("[gate] dip arm reproduces the shipped compat_splice/_%s byte-for-byte" % name)

for name in ("corrected", "deepest"):
    want = set(pd.read_csv(V1D + "/per_hap_chrom_phenotypes_%s.tsv" % name, sep="\t",
                           engine="python", nrows=1).columns)
    for arm in ARMS:
        got = set(pd.read_csv("%s/compat/compat_%s_%s/per_hap_chrom_phenotypes_%s.tsv"
                              % (OUT, arm, FRAME, name), sep="\t", nrows=1).columns)
        missing = want - got
        assert not missing, "compat_%s _%s missing v1 columns: %s" % (arm, name, sorted(missing))
    print("[gate] both arms carry all %d v1 columns for _%s" % (len(want), name))

key = ["sample", "hap_local", "chm13_chrom"]
m = built["dip"][key + ["mean_5mC", "mean_5mC_covwt", "active_len", "n_windows"]].merge(
    built["gao"][key + ["mean_5mC", "mean_5mC_covwt", "active_len", "n_windows"]],
    on=key, suffixes=("_dip", "_gao"))
assert len(m) == len(built["dip"]) == len(built["gao"]), \
    "the two arms do not cover the same (sample, hap, chrom) set: %d vs %d vs %d" \
    % (len(m), len(built["dip"]), len(built["gao"]))
for c in ("mean_5mC", "mean_5mC_covwt", "active_len", "n_windows"):
    d = (m[c + "_dip"] - m[c + "_gao"]).abs().max()
    assert d == 0, "%s differs between arms by up to %g -- the merge leaked into it" % (c, d)
print("[gate] mean_5mC, mean_5mC_covwt, active_len, n_windows IDENTICAL across arms (%d rows)"
      % len(m))

j = built["dip"][key + ["n_cdrs", "cdr_pos_offset", "cdr_fraction", "median_cdr_width"]].merge(
    built["gao"][key + ["n_cdrs", "cdr_pos_offset", "cdr_fraction", "median_cdr_width"]],
    on=key, suffixes=("_dip", "_gao"))
print()
print("rows: %d   (%d with a position in both arms)"
      % (len(j), int((j.cdr_pos_offset_dip.notna() & j.cdr_pos_offset_gao.notna()).sum())))
for c in ("n_cdrs", "cdr_pos_offset", "cdr_fraction", "median_cdr_width"):
    a, b = j[c + "_dip"], j[c + "_gao"]
    ok = a.notna() & b.notna()
    r = a[ok].corr(b[ok]) if ok.sum() > 2 else float("nan")
    print("  %-18s dip median %10.4g   gao median %10.4g   r = %.4f"
          % (c, a[ok].median(), b[ok].median(), r))
