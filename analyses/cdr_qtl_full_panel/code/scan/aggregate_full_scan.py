#!/usr/bin/env python
import pandas as pd, numpy as np, os, sys, glob
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

_D = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_D, os.path.dirname(_D)]
from phenotype_family import (PRIMARY, SUPPLEMENTARY, ALL_SCANNED, NPH, denominators, alpha,
                              headline)

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")
SCAN = sys.argv[1] if len(sys.argv) > 1 else f"{WU}/output/scan_1e5"
TAG = sys.argv[2] if len(sys.argv) > 2 else "1e5"
if not os.path.isabs(SCAN):
    SCAN = f"{WU}/{SCAN}"

files = sorted(glob.glob(f"{SCAN}/perchrom/*_perm.tsv"))
print(f"aggregating {len(files)} scan files from {SCAN}")
df = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
print(f"total tests: {len(df)}   chroms={df.chrom.nunique()}   "
      f"phenotypes={df.phenotype.nunique()}")
unknown = set(df.phenotype) - set(ALL_SCANNED)
assert not unknown, f"scan holds phenotypes absent from phenotype_family.py: {unknown}"

vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
FRAMINGS = denominators(vt)
hl_lab, hl_n, hl_a = headline(vt)
print(f"\n=== multiple-testing framings: alpha = 0.05 / (n x {NPH} primary phenotypes) ===")
for nm, n in FRAMINGS:
    print(f"  {nm:34s} n={n:6d}   alpha={alpha(n):.3e}"
          + ("   <-- HEADLINE" if nm == hl_lab else ""))

print(f"\n=== the {len(SUPPLEMENTARY)} supplementary twins are NOT in the denominator ===")
loosest = max(alpha(n) for _, n in FRAMINGS)
for tw, (pr, _why) in SUPPLEMENTARY.items():
    s = df[df.phenotype == tw]
    if not len(s):
        print(f"  {tw:16s} not present in this scan dir")
        continue
    bp = s.p_emp_perm.min()
    print(f"  {tw:16s} best p = {bp:.3g}   vs loosest alpha {loosest:.3g}   -> "
          f"{'NO HIT at any threshold considered' if bp > loosest else '*** HIT -- REVISIT THE SPLIT ***'}")

floor = df.p_emp_perm.min()
n_floor = int((df.p_emp_perm <= floor * 1.0001).sum())
print(f"\npermutation floor in this run: p={floor:.3e}   tests at the floor: {n_floor}")

cells = (df[df.p_emp_perm <= floor * 1.0001]
         .groupby(["chrom", "phenotype"]).size().reset_index(name="n_floor")
         .sort_values("n_floor", ascending=False))
esc = f"{WU}/code/escalate_from_{TAG}.txt"
with open(esc, "w") as fh:
    for r in cells.itertuples():
        fh.write(f"{r.chrom}\t{r.phenotype}\n")
print(f"cells holding a floor-saturated locus: {len(cells)} -> {esc}")

hits = df[df.p_emp_perm <= loosest].copy()
for nm, n in FRAMINGS:
    key = ("row" if nm.startswith("row") else "uniq") + ("_ADMIT" if "ADMIT" in nm else "_FULL")
    hits[f"pass_{key}"] = hits.p_emp_perm <= alpha(n)
hits["primary"] = hits.phenotype.isin(PRIMARY)
hits = hits.sort_values(["p_emp_perm", "chrom", "pos"])
out_h = f"{WU}/output/raw_hits_{TAG}.tsv"
tmp = f"{out_h}.tmp.{os.getpid()}"
hits.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
os.replace(tmp, out_h)

print(f"\n=== RAW hits, gate-coverage set (loosest alpha {loosest:.3e}): "
      f"{len(hits)} locus-phenotype tests ===")
if len(hits):
    print(f"  distinct loci {hits.sv_locus.nunique()}   genotype-groups {hits.geno_group.nunique()}")
    print(f"  by chrom:  {hits.groupby('chrom').size().to_dict()}")
    print(f"  by status: {hits.groupby('status').size().to_dict()}")
    hh = hits[hits.p_emp_perm <= hl_a]
    print(f"\n=== at the HEADLINE alpha ({hl_a:.3e}, {hl_lab}): {len(hh)} tests, "
          f"{hh.sv_locus.nunique()} loci, {hh.geno_group.nunique()} genotype-groups ===")
    print(f"  by chrom: {hh.groupby('chrom').size().to_dict()}")
    print("  genotype-groups per phenotype:")
    for ph in ALL_SCANNED:
        s = hh[hh.phenotype == ph]
        tagp = "" if ph in PRIMARY else "  [supplementary]"
        print(f"    {ph:24s} {s.geno_group.nunique():3d}  {sorted(set(s.chrom))}{tagp}")
print(f"wrote {out_h}")

out_all = f"{WU}/output/all_tests_{TAG}.tsv.gz"
tmp = f"{out_all}.tmp.{os.getpid()}"
df.sort_values(["chrom", "phenotype", "p_emp_perm"]).to_csv(
    tmp, sep="\t", index=False, float_format="%.6g", compression="gzip")
os.replace(tmp, out_all)
print(f"wrote {out_all} ({len(df)} rows)")
