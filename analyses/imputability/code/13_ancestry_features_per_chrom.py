#!/usr/bin/env python3
"""Phase B4 — ancestry-stratified features per chrom."""
import os
import sys
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"
SNP_DIR = ROOT / "snp_extraction" / "fksize500000"
PANEL = ROOT / "resources" / "1000G" / "20130606_g1k_3202_samples_ped_population.txt"
WIN_HALF = 500_000

POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def load_pop_map():
    df = pd.read_csv(PANEL, sep=r"\s+", engine="python")
    return dict(zip(df["SampleID"], df["Superpopulation"]))


def chrom_features(chrom):
    snp_npz = SNP_DIR / f"{chrom}_snp_genotypes.npz"
    snp_var = SNP_DIR / f"{chrom}_snp_variants.tsv"
    snp_sam = SNP_DIR / f"{chrom}_snp_samples.txt"
    for p in [snp_npz, snp_var, snp_sam]:
        if not p.exists():
            print(f"[SKIP] {chrom}: missing {p.name}", flush=True)
            return None

    g = np.load(snp_npz)
    arr = g[g.files[0]].astype(np.float32)
    pos = pd.read_csv(snp_var, sep="\t")["pos"].values
    samples = open(snp_sam).read().splitlines()

    pop_map = load_pop_map()
    pop_of_sample = np.array([pop_map.get(s, "OTHER") for s in samples])

    pop_idx = {pop: np.where(pop_of_sample == pop)[0] for pop in POPS}
    n_per_pop = {pop: len(pop_idx[pop]) for pop in POPS}
    print(f"[{chrom}] SNPs={arr.shape[0]} samples={arr.shape[1]} "
          f"per-pop n={n_per_pop}", flush=True)
    if any(n < 10 for n in n_per_pop.values()):
        print(f"[FAIL] {chrom}: at least one super-pop has <10 samples", flush=True)
        return None

    miss_mask = ~np.isfinite(arr) | (arr < 0)
    arr_m = np.where(miss_mask, np.nan, arr)
    p = {}
    n_eff = {}
    for pop in POPS:
        idx = pop_idx[pop]
        sub = arr_m[:, idx]
        n_obs = np.sum(~np.isnan(sub), axis=1)
        sum_d = np.nansum(sub, axis=1)
        af = np.where(n_obs > 0, sum_d / (2 * n_obs), np.nan)
        p[pop] = af
        n_eff[pop] = n_obs

    af_stack = np.vstack([p[pop] for pop in POPS])
    af_skew = np.nanmax(af_stack, axis=0) - np.nanmin(af_stack, axis=0)

    pair_num = {}
    pair_den = {}
    for i, j in combinations(POPS, 2):
        pi = p[i]; pj = p[j]
        ni = np.maximum(n_eff[i], 1); nj = np.maximum(n_eff[j], 1)
        num = (pi - pj) ** 2 - pi * (1 - pi) / np.maximum(ni - 1, 1) \
                              - pj * (1 - pj) / np.maximum(nj - 1, 1)
        den = pi * (1 - pj) + pj * (1 - pi)
        pair_num[(i, j)] = num
        pair_den[(i, j)] = den

    feat = pd.read_csv(OUT / "per_sv_features_v3.tsv", sep="\t",
                       low_memory=False, usecols=["chrom", "sv_index", "pos"])
    feat = feat[feat["chrom"] == chrom].copy()
    rows = []
    for _, sv in feat.iterrows():
        sv_pos = int(sv["pos"])
        mask = (np.abs(pos - sv_pos) <= WIN_HALF)
        n_in = int(mask.sum())
        if n_in == 0:
            rows.append({"chrom": chrom, "sv_index": int(sv["sv_index"]),
                         **{f"mean_AF_{p_}": np.nan for p_ in POPS},
                         "AF_skew_mean_window": np.nan,
                         "max_FST_pair_window": np.nan,
                         "n_pop_snps_window": 0})
            continue
        d = {"chrom": chrom, "sv_index": int(sv["sv_index"])}
        for pop in POPS:
            d[f"mean_AF_{pop}"] = float(np.nanmean(p[pop][mask]))
        d["AF_skew_mean_window"] = float(np.nanmean(af_skew[mask]))
        max_fst = -np.inf
        for ij, num in pair_num.items():
            den = pair_den[ij]
            num_w = num[mask]; den_w = den[mask]
            valid = np.isfinite(num_w) & np.isfinite(den_w)
            if valid.sum() < 5:
                continue
            num_sum = np.nansum(num_w[valid])
            den_sum = np.nansum(den_w[valid])
            if den_sum > 0:
                fst = num_sum / den_sum
                if fst > max_fst:
                    max_fst = fst
        d["max_FST_pair_window"] = float(max_fst) if np.isfinite(max_fst) else np.nan
        d["n_pop_snps_window"] = n_in
        rows.append(d)

    out = pd.DataFrame(rows)
    return out


def main():
    chrom = os.environ.get("CHROM_ARG")
    if chrom is None:
        print("[FAIL] set CHROM_ARG=chrN", flush=True)
        sys.exit(1)
    chunk_dir = OUT / "ancestry_chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    df = chrom_features(chrom)
    if df is None or len(df) == 0:
        print(f"[FAIL] {chrom}: no rows produced", flush=True)
        sys.exit(1)
    out_path = chunk_dir / f"ancestry_features_{chrom}.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {out_path} ({len(df)} rows)", flush=True)
    print(df[["mean_AF_AFR", "mean_AF_EUR", "AF_skew_mean_window",
              "max_FST_pair_window", "n_pop_snps_window"]].describe().to_string(),
          flush=True)


if __name__ == "__main__":
    main()
