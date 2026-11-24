#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
SUFFIX = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
SV_BASE = ROOT / "sv_extraction" / "decompose" / SUFFIX
PANEL = ROOT / "resources" / "1000G" / "20130606_g1k_3202_samples_ped_population.txt"


def main():
    print("=== chrX ploidy / encoding check ===")
    npz = SV_BASE / "chrX" / "chrX_sv_genotypes.npz"
    sam = SV_BASE / "chrX" / "chrX_sv_samples.txt"

    g = np.load(npz)
    arr = g[g.files[0]]
    samples = open(sam).read().splitlines()

    print(f"shape: {arr.shape}")
    print(f"dtype: {arr.dtype}")
    print(f"unique values: {np.unique(arr)[:20]}")
    print(f"max value: {arr.max()}, min: {arr.min()}")

    panel = pd.read_csv(PANEL, sep=r"\s+", engine="python")
    panel_sex = dict(zip(panel["SampleID"], panel["Sex"]))

    n_male = 0; n_female = 0; n_unknown = 0
    sex_array = []
    for s in samples:
        sex = panel_sex.get(s, "?")
        sex_array.append(sex)
        if sex == 1: n_male += 1
        elif sex == 2: n_female += 1
        else: n_unknown += 1
    print(f"sample sexes: male={n_male} female={n_female} unknown={n_unknown}")

    sample_max = arr.max(axis=0)
    male_idx = [i for i, x in enumerate(sex_array) if x == 1]
    female_idx = [i for i, x in enumerate(sex_array) if x == 2]
    print(f"\nper-sample max value distribution:")
    print(f"  males (n={n_male}): max ranges [{sample_max[male_idx].min()}, "
          f"{sample_max[male_idx].max()}]  unique: {np.unique(sample_max[male_idx])}")
    print(f"  females (n={n_female}): max ranges [{sample_max[female_idx].min()}, "
          f"{sample_max[female_idx].max()}]  unique: {np.unique(sample_max[female_idx])}")

    if sample_max[male_idx].max() <= 1 and sample_max[female_idx].max() <= 2:
        print("\n>>> CONCLUSION: chrX is correctly ploidy-aware "
              "(males capped at 1, females at 2)")
    elif arr.max() == 2:
        print("\n>>> CONCLUSION: chrX uses uniform diploid encoding "
              "(may double-count male chrX dosage — could affect AF calculations)")
    elif arr.max() == 1:
        print("\n>>> CONCLUSION: chrX uses uniform haploid encoding "
              "(treated as haploid for everyone)")


if __name__ == "__main__":
    main()
