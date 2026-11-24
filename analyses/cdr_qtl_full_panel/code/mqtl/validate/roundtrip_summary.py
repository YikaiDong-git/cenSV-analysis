#!/usr/bin/env python
"""Summarise the CDR coordinate round-trip across haplotypes."""
import glob, os
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")
D = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{WU}/output/mqtl/validation/roundtrip/*.tsv"))],
              ignore_index=True)
print(f"haplotypes {D.tag.nunique()}   CDR rows {len(D)}")
print("\nlift outcome:")
print(D.reason.value_counts().head(8).to_string())

ok = D[D.reason == "ok"].copy()
print(f"\nlifted through the BAM: {len(ok)}/{len(D)} = {100*len(ok)/len(D):.1f}%")
print(f"chromosome agreement with the chain: {100*ok.chrom_match.mean():.3f}%")

m = ok[ok.chrom_match].copy()
print(f"\non {len(m)} CDRs landing on the SAME chromosome:")
for lab, g in [("ALL", m)] + [(f"label {'AGREES' if k else 'INVERTED'}", gg)
                              for k, gg in m.groupby("label_agrees")]:
    print(f"  {lab:18s} n={len(g):6d}  haps={g.tag.nunique():3d}  "
          f"median d_start {g.d_start.median():+,.0f}  "
          f"|d|<=1kb {100*(g.d_start.abs()<=1000).mean():5.1f}%  "
          f"|d|<=10kb {100*(g.d_start.abs()<=10000).mean():5.1f}%  "
          f"p95|d| {g.d_start.abs().quantile(0.95):,.0f}")

# A per-haplotype median offset is the signature of a WRONGLY PAIRED haplotype: arithmetic errors
ph = m.groupby(["tag", "label_agrees"]).agg(n=("d_start", "size"),
                                            med=("d_start", "median"),
                                            frac1kb=("d_start", lambda s: (s.abs() <= 1000).mean())
                                            ).reset_index().sort_values("frac1kb")
print(f"\nper-haplotype, WORST 10 by fraction within 1 kb:")
print(ph.head(10).to_string(index=False))
print(f"\nper-haplotype frac-within-1kb: median {ph.frac1kb.median():.4f}  min {ph.frac1kb.min():.4f}")
print(f"per-haplotype median offset: median {ph.med.median():+,.0f}  "
      f"max|median| {ph.med.abs().max():,.0f}")

m["wr"] = m.bam_width / m.chain_width.replace(0, np.nan)
print(f"\nwidth ratio bam/chain: median {m.wr.median():.4f}  "
      f"IQR {m.wr.quantile(.25):.4f}-{m.wr.quantile(.75):.4f}")
os.makedirs(f"{WU}/output/mqtl/validation", exist_ok=True)
ph.to_csv(f"{WU}/output/mqtl/validation/roundtrip_by_hap.tsv", sep="\t", index=False,
          float_format="%.6g")
print("\ndone")
