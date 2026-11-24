#!/usr/bin/env python
"""How many of the 11,129 are low-frequency in OUR OWN cohort (svfreq / af < 0.05)?"""
import pandas as pd, numpy as np, os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
EN = f"{ROOT}/results/elastic_net/decompose/cosine_nj"

vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
N_IND = 188
print(f"=== af of the 11,129, computed on the {N_IND} HPRC individuals (376 alleles) ===")
print(f"  af   : min={vt.af.min():.4f}  1%={vt.af.quantile(.01):.4f}  median={vt.af.median():.4f}  max={vt.af.max():.4f}")
print(f"  maf  : min={vt.maf.min():.4f}  median={vt.maf.median():.4f}")
print(f"  smallest non-zero af corresponds to {vt.af.min()*2*N_IND:.0f} alt alleles")
for thr in [0.01, 0.02, 0.05, 0.10]:
    n = int((vt.maf < thr).sum())
    print(f"  maf < {thr:.2f} : {n:6d} / {len(vt)} rows ({100*n/len(vt):5.1f}%)   "
          f"geno_groups {vt[vt.maf < thr].geno_group.nunique():5d} / {vt.geno_group.nunique()}")

print("\n=== by ADMIT status ===")
for st, g in vt.groupby("status"):
    print(f"  {st:24s} n={len(g):5d}  maf median={g.maf.median():.4f}  "
          f"maf<0.05: {100*(g.maf<0.05).mean():5.1f}%   maf<0.01: {100*(g.maf<0.01).mean():5.1f}%")

print("\n=== what an af >= 0.05 filter would leave ===")
keep = vt[vt.maf >= 0.05]
print(f"  rows {len(keep):,} / {len(vt):,}   geno_groups {keep.geno_group.nunique():,} / {vt.geno_group.nunique():,}")
print(f"  chroms retained: {keep.chrom.nunique()} / {vt.chrom.nunique()}")
lost = sorted(set(vt.chrom) - set(keep.chrom))
print(f"  chroms lost entirely: {lost if lost else 'none'}")

print("\n=== do the surviving gated leads pass a 0.05 floor? ===")
gp = f"{WU}/output/gates/full_hits_gate_1e7.tsv"
if os.path.exists(gp):
    g = pd.read_csv(gp, sep="\t")
    s = g[g.overall == "SURVIVES"].drop_duplicates(["chrom", "pos", "svtype", "svlen", "af"])
    for r in s.sort_values(["chrom", "af"]).itertuples():
        mafv = min(r.af, 1 - r.af)
        print(f"  {r.chrom:5s} {r.svtype}-{r.svlen}@{r.pos:<10d} af={r.af:.4f} maf={mafv:.4f}  "
              f"{'PASSES 0.05' if mafv >= 0.05 else '*** BELOW 0.05 ***'}")

print("\n=== was a frequency floor ever applied upstream? (universe from sv_metrics) ===")
tot, below = 0, 0
mins = []
for ch in sorted([f"chr{i}" for i in range(1, 23)] + ["chrX"], key=lambda c: (len(c), c)):
    for s in ["midpoint1000", "midpoint500", "edge1000"]:
        p = f"{EN}/{s}_{SUF}/{ch}/sv_metrics.tsv"
        if os.path.exists(p):
            m = pd.read_csv(p, sep="\t", usecols=["svfreq"])
            v = pd.to_numeric(m.svfreq, errors="coerce").dropna()
            tot += len(v); below += int((v < 0.05).sum()); mins.append(v.min())
            break
print(f"  sv_metrics rows: {tot:,}   svfreq < 0.05: {below:,} ({100*below/max(tot,1):.1f}%)   "
      f"global min svfreq = {min(mins) if mins else float('nan'):.4f}")
print("  (a hard floor would show as a sharp minimum at the threshold; a smooth tail down to ~1/376 means none)")
