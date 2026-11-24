#!/usr/bin/env python3
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
    base = pd.read_csv(OUT / "per_sv_features_v3.tsv", sep="\t", low_memory=False)
    print(f"[INFO] base v3 table: {len(base)} x {base.shape[1]}", flush=True)

    g23 = pd.read_csv(OUT / "per_sv_genome_features.tsv", sep="\t", low_memory=False)
    print(f"[INFO] genome (GC + sat): {len(g23)} rows; "
          f"sat_family value counts:\n{g23['sat_family'].value_counts().to_string()}",
          flush=True)
    base["sv_index"] = base["sv_index"].astype(int)
    g23["sv_index"] = g23["sv_index"].astype(int)
    merged = base.merge(g23, on=["chrom", "sv_index"], how="left",
                         suffixes=("", "_g23"))

    fam_keep = (g23["sat_family"].value_counts()
                                  .pipe(lambda s: s[s >= 50])
                                  .index.tolist())
    print(f"[INFO] sat_family one-hot keep: {fam_keep}", flush=True)
    for f in fam_keep:
        merged[f"sat_{f}"] = (merged["sat_family"] == f).astype(int)

    sv_anc_path = OUT / "sv_ancestry_features.tsv"
    if sv_anc_path.exists():
        sv_anc = pd.read_csv(sv_anc_path, sep="\t", low_memory=False)
        sv_anc["sv_index"] = sv_anc["sv_index"].astype(int)
        print(f"[INFO] sv_ancestry rows: {len(sv_anc)}", flush=True)
        merged = merged.merge(sv_anc, on=["chrom", "sv_index"], how="left",
                              suffixes=("", "_svanc"))
    else:
        print(f"[WARN] {sv_anc_path} not found; skipping SV-level ancestry merge",
              flush=True)

    chunks = sorted(glob.glob(str(OUT / "ancestry_chunks" / "ancestry_features_chr*.tsv")))
    if chunks:
        print(f"[INFO] SNP-window ancestry chunks: {len(chunks)}", flush=True)
        anc_frames = [pd.read_csv(c, sep="\t") for c in chunks]
        anc = pd.concat(anc_frames, ignore_index=True)
        anc["sv_index"] = anc["sv_index"].astype(int)
        merged = merged.merge(anc, on=["chrom", "sv_index"], how="left",
                              suffixes=("", "_snpwin"))

    print(f"[INFO] merged v3+genome: {len(merged)} x {merged.shape[1]}", flush=True)

    cov_cols = ["gc_100kb", "sat_active", "sat_local_complexity",
                "af_AFR_sv", "af_EUR_sv", "af_skew_sv", "max_FST_pair_sv",
                "mean_AF_AFR", "mean_AF_EUR", "max_FST_pair_window"]
    for c in cov_cols:
        if c in merged.columns:
            print(f"  {c}: {merged[c].notna().mean()*100:.1f}% non-NA", flush=True)

    out_path = OUT / "per_sv_features_with_genome.tsv"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"\n[OK] wrote {out_path}  ({len(merged)} x {merged.shape[1]})", flush=True)

    geno_only = merged[["chrom", "sv_index", "pos"] + [c for c in g23.columns
                       if c not in {"chrom", "sv_index"}] + [
                       f"sat_{f}" for f in fam_keep] + [c for c in anc.columns
                       if c not in {"chrom", "sv_index"}]].copy()
    geno_only.to_csv(OUT / "genome_features_only.tsv", sep="\t", index=False)
    print(f"[OK] wrote {OUT/'genome_features_only.tsv'}", flush=True)


if __name__ == "__main__":
    main()
