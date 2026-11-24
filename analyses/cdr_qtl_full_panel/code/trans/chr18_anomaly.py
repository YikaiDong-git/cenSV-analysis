#!/usr/bin/env python
"""Diagnose chr18's two anomalies, which FINDINGS 8h/8j record as probably ONE anomaly, undiagnosed."""
import glob, os
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
CANON = f"{P1}/analyses/cdr_qtl"
AUD = f"{P1}/analyses/cdr_calls"
WU = f"{P1}/analyses/cdr_qtl_full_panel"
OUT = f"{WU}/output/trans"
os.makedirs(OUT, exist_ok=True)

files = sorted(glob.glob(f"{WU}/output/mqtl/per_hap_bam/*.tsv.gz"))
print(f"per-haplotype lifts: {len(files)}")
acc = []
for i, f in enumerate(files):
    D = pd.read_csv(f, sep="\t")
    D = D[D.window_bp == 100000]
    g = D.groupby("chrom").agg(n_sv=("sv_locus", "size"),
                               blocks=("n_blocks", lambda s: (s > 0).sum()),
                               bins=("n_bins", lambda s: (s > 0).sum()),
                               meth=("meth", lambda s: s.notna().sum())).reset_index()
    g["tag"] = os.path.basename(f).replace(".tsv.gz", "")
    acc.append(g)
    if (i + 1) % 80 == 0:
        print(f"  {i+1}/{len(files)}")
A = pd.concat(acc, ignore_index=True)
A.to_csv(f"{OUT}/chr18_anomaly_by_hap.tsv", sep="\t", index=False)

S = A.groupby("chrom").agg(haps=("tag", "nunique"), sv_cells=("n_sv", "sum"),
                           with_blocks=("blocks", "sum"), with_bins=("bins", "sum"),
                           with_meth=("meth", "sum")).reset_index()
S["pct_blocks"] = 100 * S.with_blocks / S.sv_cells
S["pct_bins"] = 100 * S.with_bins / S.sv_cells
S["pct_meth"] = 100 * S.with_meth / S.sv_cells
S["pct_aligned_but_no_meth"] = 100 * (S.with_blocks - S.with_bins) / S.with_blocks.replace(0, np.nan)
S = S.sort_values("pct_meth")
S.to_csv(f"{OUT}/chr18_anomaly_by_chrom.tsv", sep="\t", index=False, float_format="%.4g")
print("\nPER CHROMOSOME (window = +/-100 kb, pooled over haplotypes)")
print(S[["chrom", "haps", "sv_cells", "pct_blocks", "pct_bins", "pct_meth",
         "pct_aligned_but_no_meth"]].to_string(index=False))

PC = pd.read_csv(f"{CANON}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
PD_ = pd.read_csv(f"{CANON}/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t")
rows = []
for c in sorted(set(PC.chm13_chrom.dropna())):
    a = PC[PC.chm13_chrom == c]
    d = PD_[PD_.chm13_chrom == c]
    rows.append(dict(chrom=c, hap_rows=len(a),
                     indiv_any_row=a["sample"].nunique(),
                     indiv_active_len=a.dropna(subset=["active_len"])["sample"].nunique(),
                     indiv_n_cdrs_gt0=a[a.n_cdrs > 0]["sample"].nunique(),
                     indiv_cdr_pos=a.dropna(subset=["cdr_pos_offset"])["sample"].nunique(),
                     indiv_deepest=(d.dropna(subset=["cdr_pos_offset_deepest"])["sample"].nunique()
                                    if "cdr_pos_offset_deepest" in d.columns else 0),
                     indiv_mean5mC=a.dropna(subset=["mean_5mC"])["sample"].nunique(),
                     median_active_len=a.active_len.median(),
                     median_n_cdrs=a.n_cdrs.median()))
P = pd.DataFrame(rows).sort_values("indiv_deepest")
P.to_csv(f"{OUT}/chr18_phenotype_availability.tsv", sep="\t", index=False, float_format="%.6g")
print("\nPHENOTYPE AVAILABILITY per chromosome (individuals, not haplotypes)")
print(P.to_string(index=False))

print("\nchr18 vs the panel:")
c18 = P[P.chrom == "chr18"]
if len(c18):
    r = c18.iloc[0]
    med = P.median(numeric_only=True)
    print(f"  individuals with ANY chr18 row      {r.indiv_any_row:.0f}  (panel median {med.indiv_any_row:.0f})")
    print(f"  with active_len defined             {r.indiv_active_len:.0f}")
    print(f"  with n_cdrs > 0                     {r.indiv_n_cdrs_gt0:.0f}")
    print(f"  with cdr_pos_offset_deepest         {r.indiv_deepest:.0f}  (panel median {med.indiv_deepest:.0f})")
    print(f"  with mean_5mC                       {r.indiv_mean5mC:.0f}")
    print(f"  median active_len                   {r.median_active_len:,.0f} bp "
          f"(panel median {med.median_active_len:,.0f})")
    print(f"  median n_cdrs                       {r.median_n_cdrs:.1f} (panel median {med.median_n_cdrs:.1f})")
    print("\n  READING: if 'ANY row' is near the panel median but 'n_cdrs > 0' collapses, the array is")
    print("  assembled and measured but CARRIES NO DETECTABLE DIP -- a biological/detection property,")
    print("  not a pipeline failure, and the 54 individuals are those where a dip WAS called.")
print("\ndone")
