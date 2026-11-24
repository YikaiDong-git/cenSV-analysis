#!/usr/bin/env python3
"""Phase D verification — diploid vs haploid encoding of SV genotype matrix."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
SUFFIX = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
SV_BASE = ROOT / "sv_extraction" / "decompose" / SUFFIX
BEST8 = ROOT / "results" / "imputation" / "decompose" / "cosine_nj" / f"best8_{SUFFIX}"
OUT = ROOT / "analyses" / "imputability" / "output"


def main():
    chrom = "chr1"
    sv_npz = SV_BASE / chrom / f"{chrom}_sv_genotypes.npz"
    sv_sam = SV_BASE / chrom / f"{chrom}_sv_samples.txt"

    g = np.load(sv_npz)
    arr = g[g.files[0]]
    samples = open(sv_sam).read().splitlines()

    print(f"chr1 SV NPZ: shape {arr.shape}", flush=True)
    print(f"chr1 sample list: {len(samples)} entries", flush=True)
    print(f"value distribution: unique={np.unique(arr)[:20]}", flush=True)
    n_per_value = pd.Series(arr.flatten()).value_counts().sort_index()
    print(f"value counts (first 10):\n{n_per_value.head(10).to_string()}", flush=True)
    print(f"max value: {arr.max()}, min: {arr.min()}", flush=True)
    print(f"matrix dtype: {arr.dtype}", flush=True)

    print(f"\nsv_index=0 first 30 dosages: {arr[0, :30]}", flush=True)
    print(f"sv_index=0 dosage sum: {arr[0].sum()}", flush=True)
    print(f"sv_index=0 nonzero count: {(arr[0] > 0).sum()}", flush=True)
    print(f"sv_index=0 max value: {arr[0].max()}", flush=True)

    summary_path = BEST8 / chrom / f"{chrom}_imputation_summary.tsv"
    summ = pd.read_csv(summary_path, sep="\t")
    sv0 = summ[summ["sv_index"] == 0].iloc[0]
    print(f"\n=== best8 imputation_summary chr1 sv_index=0 ===", flush=True)
    print(sv0[["sv_index", "chrom", "pos", "svtype", "svlen",
                "svfreq", "imputed_mean", "imputed_std"]].to_string(), flush=True)

    af_haploid = arr[0].sum() / arr.shape[1]
    af_diploid = arr[0].sum() / (2 * arr.shape[1])
    print(f"\n=== AF computation for sv_index=0 ===", flush=True)
    print(f"if HAPLOID encoding (n_hap = {arr.shape[1]}): AF = {af_haploid:.4f}",
          flush=True)
    print(f"if DIPLOID encoding (n_diploid = {arr.shape[1]}): AF = {af_diploid:.4f}",
          flush=True)
    print(f"\nsvfreq from imputation_summary: {sv0['svfreq']:.4f}", flush=True)
    print(f"\nsvfreq matches HAPLOID encoding: "
          f"{abs(af_haploid - sv0['svfreq']) < 0.01}",  flush=True)
    print(f"svfreq matches DIPLOID encoding: "
          f"{abs(af_diploid - sv0['svfreq']) < 0.01}",  flush=True)

    if abs(af_haploid - sv0['svfreq']) < 0.01:
        print(f"\n>>> CONCLUSION: SV NPZ is HAPLOID-encoded, "
              f"samples list = {arr.shape[1]} haplotypes (not diploid samples)",
              flush=True)
        print(f">>> My 13b_sv_ancestry.py is computing AF as sum/(2*n_obs), "
              f"which UNDERESTIMATES AF by 2x.", flush=True)
        print(f">>> Fix: change normalization from /2 to /1 for haploid matrix.",
              flush=True)
    elif abs(af_diploid - sv0['svfreq']) < 0.01:
        print(f"\n>>> CONCLUSION: SV NPZ is DIPLOID-encoded, samples list = "
              f"{arr.shape[1]} diploid individuals", flush=True)
        print(f">>> My 13b_sv_ancestry.py treatment is CORRECT.", flush=True)
    else:
        print(f"\n>>> NEITHER encoding matches svfreq exactly", flush=True)
        print(f">>> Need deeper investigation", flush=True)


if __name__ == "__main__":
    main()
