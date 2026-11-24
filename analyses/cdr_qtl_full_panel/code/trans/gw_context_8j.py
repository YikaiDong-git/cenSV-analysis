#!/usr/bin/env python
"""Place 8j's surviving chr18 association in the genome-wide frame."""
import glob
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")
GW = f"{WU}/output/trans/gw"

gmax = np.max(np.vstack([np.load(f) for f in sorted(glob.glob(f"{GW}/*/*_maxT.npy"))]), axis=0)
NP = len(gmax)
R = pd.read_csv(f"{WU}/output/trans/gw_results.tsv", sep="\t")
F = pd.read_csv(f"{WU}/output/trans/focused_trans_results.tsv", sep="\t")
hits = F[(F.relation == "trans") & (F.bh_q_trans_escalated < 0.05)]

print(f"global max-T null over {NP:,} shared permutations, 192 cells, 2,045,280 trans tests")
for q in (0.05, 0.25, 0.5, 0.75, 0.95, 0.99):
    print(f"  quantile {q:.2f}: {np.quantile(gmax, q):.3f}")
print(f"  max {gmax.max():.3f}")

print("\n8j's BH-significant rows, re-scored against the GENOME-WIDE family:")
rows = []
for h in hits.itertuples():
    m = R[(R.sv_locus == h.sv_locus) & (R.phenotype == h.phenotype)
          & (R.target_chrom == h.target_chrom)]
    fw = float((gmax >= abs(h.t)).sum() + 1) / (NP + 1)
    pct = float((gmax < abs(h.t)).mean())
    rows.append(dict(sv_locus=h.sv_locus, phenotype=h.phenotype, target=h.target_chrom,
                     n=h.n, t=h.t,
                     focused_fwer_960=h.fwer_global, focused_bh=h.bh_q_trans_escalated,
                     gw_fwer=fw, pct_of_gw_null_below=100 * pct,
                     in_gw_candidates=bool(len(m))))
D = pd.DataFrame(rows).sort_values("gw_fwer")
print(D.to_string(index=False))
D.to_csv(f"{WU}/output/trans/gw_context_for_8j.tsv", sep="\t", index=False, float_format="%.6g")

top = D.iloc[0]
print(f"\nThe chr18 association (|t| = {abs(top.t):.3f}) sits below "
      f"{100 - top.pct_of_gw_null_below:.1f}% of the genome-wide null's maxima "
      f"(median null max = {np.median(gmax):.3f}).")
print("READING: significant within a pre-specified 5-SV family, ordinary against the full catalogue.")
print("\ndone")
