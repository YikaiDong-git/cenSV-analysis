"""Aggregate the re-scan, and check it against the shipped scan where it must agree."""
import glob, os
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

RD = (CENSV_ROOT + "/discovery"
      "/cdr_phenotypes")
WU = (CENSV_ROOT + "/discovery"
      "/cdr_qtl_full_panel")
FLOOR = 1e-5
CLEAN = ["chr2", "chr5", "chr6", "chr8", "chr11", "chr12", "chr18", "chr20", "chr21", "chr22", "chrX"]

def load(d, frame=None):
    out = []
    for f in sorted(glob.glob(d + "/perchrom/*_perm.tsv")):
        try:
            t = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if len(t):
            t["frame"] = frame if frame else "shipped"
            out.append(t)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()

ship = load(WU + "/output/scan/scan_1e5")
mine = pd.concat([load(RD + "/output/scan/scan_1e5/" + f, f) for f in ("span", "block", "splice")],
                 ignore_index=True)
print("shipped rows %s | re-scan rows %s" % (format(len(ship), ","), format(len(mine), ",")))

key = ["chrom", "sv_locus", "phenotype"]
POS = ["cdr_pos_offset", "cdr_pos_offset_deepest", "cdr_pos_iqr", "cdr_pos_var"]
sp = mine[mine.frame == "span"]
m = ship.merge(sp, on=key, suffixes=("_ship", "_new"))
ctl = m[m.chrom.isin(CLEAN) & m.phenotype.isin(POS)]
print("\n=== POSITIVE CONTROL: span frame vs shipped, position phenotypes, clean chromosomes ===")
if len(ctl):
    d = (ctl.t_obs_ship - ctl.t_obs_new).abs()
    print("n compared        : %s" % format(len(ctl), ","))
    print("max |dt|          : %.3e" % d.max())
    print("median |dt|       : %.3e" % d.median())
    print("r(t_ship, t_new)  : %.8f" % np.corrcoef(ctl.t_obs_ship, ctl.t_obs_new)[0, 1])
    print("VERDICT           : %s" % ("REPRODUCED" if d.max() < 1e-6 else
          ("near-identical" if d.max() < 0.05 else "*** DIVERGED ***")))
else:
    print("no overlapping cells to compare")

moved = m[m.chrom.isin(["chr3", "chr4", "chr10"]) & m.phenotype.isin(["cdr_fraction", "cdrs_per_mb"])]
if len(moved):
    d2 = (moved.t_obs_ship - moved.t_obs_new).abs()
    print("\ncontrast -- chr3/4/10 density (where correction SHOULD bite): "
          "n=%d, max |dt| = %.3f, median %.3f" % (len(moved), d2.max(), d2.median()))

print("\n=== CELLS AT THE 1e-5 PERMUTATION FLOOR (escalation candidates) ===")
fl = mine[mine.p_emp_perm <= FLOOR]
if len(fl):
    g = (fl.groupby(["frame", "chrom", "phenotype"])
           .agg(n_sv=("sv_locus", "nunique"), max_abs_t=("abs_t_obs", "max"))
           .reset_index().sort_values("max_abs_t", ascending=False))
    print(g.to_string(index=False))
    print("\nchromosomes reaching the floor, by frame:")
    for fr in ("span", "block", "splice"):
        s = sorted(fl[fl.frame == fr].chrom.unique())
        print("  %-7s %s" % (fr, s if s else "(none)"))
else:
    print("none")

shipfl = ship[ship.p_emp_perm <= FLOOR]
print("\nshipped scan, same criterion: %s" % sorted(shipfl.chrom.unique()))

os.makedirs(RD + "/output/agg", exist_ok=True)
mine.to_csv(RD + "/output/agg/rescan_all_frames.tsv.gz", sep="\t", index=False,
            compression="gzip")
if len(fl):
    fl.to_csv(RD + "/output/agg/floor_cells.tsv", sep="\t", index=False, float_format="%.6g")
print("\n[written] output/agg/rescan_all_frames.tsv.gz, floor_cells.tsv")
