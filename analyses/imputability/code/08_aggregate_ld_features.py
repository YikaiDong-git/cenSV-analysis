#!/usr/bin/env python3
"""Aggregate per-chrom LD-feature chunks into a single TSV + best-of-8 collapse."""
import sys
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"


def main():
    chunks = sorted(glob.glob(str(OUT / "ld_chunks" / "ld_features_per_sv_strategy_chr*.tsv")))
    if not chunks:
        print("[FAIL] no chunk files; run 07 per chrom first", flush=True)
        sys.exit(1)
    print(f"[INFO] aggregating {len(chunks)} chunks", flush=True)
    frames = [pd.read_csv(c, sep="\t") for c in chunks]
    full = pd.concat(frames, ignore_index=True)
    out_long = OUT / "ld_features_per_sv_strategy.tsv"
    full.to_csv(out_long, sep="\t", index=False)
    print(f"[OK] long → {out_long} ({len(full)} rows)", flush=True)

    best = (full.dropna(subset=["top_tag_r2"])
                .sort_values("top_tag_r2", ascending=False)
                .drop_duplicates(subset=["chrom", "sv_index"], keep="first"))
    out_best = OUT / "ld_features_best.tsv"
    best.to_csv(out_best, sep="\t", index=False)
    print(f"[OK] best  → {out_best} ({len(best)} rows)", flush=True)

    feat = pd.read_csv(OUT / "per_sv_features.tsv", sep="\t", low_memory=False)
    print(f"[INFO] base feature table: {len(feat)} SVs × {feat.shape[1]} cols",
          flush=True)
    feat["sv_index"] = feat["sv_index"].astype(int)
    best_keep = best[["chrom", "sv_index", "top_tag_r2", "n_high_ld_r2_05",
                      "max_r2_within_50kb", "max_r2_within_500kb", "max_r2_within_5Mb",
                      "mean_pair_r2_window", "n_indep_snps_window",
                      "n_window_snps_actual"]].copy()
    merged = feat.merge(best_keep, on=["chrom", "sv_index"], how="left")
    out_merged = OUT / "per_sv_features_with_ld.tsv"
    merged.to_csv(out_merged, sep="\t", index=False)
    print(f"[OK] merged → {out_merged} ({len(merged)} SVs × {merged.shape[1]} cols)",
          flush=True)
    cov = merged[["top_tag_r2", "n_indep_snps_window"]].notna().mean()
    print(f"\nLD-feature coverage on best8 SV table:\n{cov.to_string()}", flush=True)


if __name__ == "__main__":
    main()
