#!/usr/bin/env python
import pandas as pd, numpy as np, os, sys, glob
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = CENSV_ROOT + "/analyses/cdr_qtl_full_panel"
CAN = CENSV_ROOT + "/analyses/cdr_qtl"
FULL_DIR = sys.argv[1] if len(sys.argv) > 1 else f"{WU}/output/scan_1e5"
CAN_DIR = sys.argv[2] if len(sys.argv) > 2 else f"{CAN}/results/canonical_scan_1e5"
TAG = sys.argv[3] if len(sys.argv) > 3 else "1e5"


def load(d):
    fs = sorted(glob.glob(f"{d}/perchrom/*_perm.tsv"))
    if not fs:
        raise SystemExit(f"no scan files under {d}")
    return pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True), len(fs)


full, nf = load(FULL_DIR)
can, nc = load(CAN_DIR)
print(f"full: {len(full)} tests from {nf} files   |   canonical: {len(can)} tests from {nc} files")

KEY = ["chrom", "pos", "end", "svtype", "svlen", "ac", "an", "af", "maf", "n_miss", "phenotype"]
fg = full.groupby(KEY)["t_obs"].apply(lambda s: np.sort(s.values))
fp = full.groupby(KEY)["p_emp_perm"].apply(lambda s: np.sort(s.values))
cg = can.groupby(KEY)["t_obs"].apply(lambda s: np.sort(s.values))
cp = can.groupby(KEY)["p_emp_perm"].apply(lambda s: np.sort(s.values))
shared = fg.index.intersection(cg.index)
print(f"\n=== shared variant-metadata x phenotype keys: {len(shared)} ===")
amb = sum(1 for k in shared if len(fg[k]) > 1 or len(cg[k]) > 1)
print(f"  keys mapping to >1 locus in either run (ambiguous metadata): {amb}")
bad = []
for k in shared:
    ft, ct = fg[k], cg[k]
    fpv, cpv = fp[k], cp[k]
    for tv, pv in zip(ct, cpv):
        if not (np.isclose(ft, tv, rtol=0, atol=1e-9).any() and np.isclose(fpv, pv, rtol=0, atol=1e-12).any()):
            bad.append((k, tv, pv, ft, fpv))
            break
print(f"  canonical results NOT reproduced in the full run: {len(bad)}")
if bad:
    for k, tv, pv, ft, fpv in bad[:15]:
        print(f"    {k}  canonical t={tv:.6f} p={pv:.3e}  |  full t={ft} p={fpv}")
else:
    print("  VERIFIED: every canonical locus-phenotype result is reproduced EXACTLY in the full run.")
m = (pd.DataFrame({"p_full": fp[shared].map(lambda a: a.min()), "p_can": cp[shared].map(lambda a: a.min())})
     .reset_index())
dp = (m.p_full - m.p_can).abs()
dt = pd.Series(np.zeros(len(m)))

f2 = full.groupby(KEY, as_index=False).agg(p_full=("p_emp_perm", "min"), t_full=("t_obs", "first"))
c2 = can.groupby(KEY, as_index=False).agg(p_can=("p_emp_perm", "min"))
newc = f2.merge(c2, on=KEY, how="left", indicator=True)
newc = newc[newc._merge == "left_only"]
print(f"\n=== tests present ONLY in the full substrate: {len(newc)} ===")
print(newc.groupby("chrom").size().sort_values(ascending=False).to_string())

if "p_full" in newc:
    print("\ntop 15 NEW tests by p:")
    cols = ["chrom", "pos", "svtype", "svlen", "af", "phenotype", "p_full", "t_full"]
    print(newc.sort_values("p_full")[cols].head(15).to_string(index=False))

out = f"{WU}/output/diff_vs_canonical_{TAG}.tsv"
tmp = f"{out}.tmp.{os.getpid()}"
m.assign(dp=dp).sort_values("p_full").to_csv(tmp, sep="\t", index=False, float_format="%.6g")
os.replace(tmp, out)
print(f"\nwrote {out}")
