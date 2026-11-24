#!/usr/bin/env python3
"""Phase B4 (revised) — SV-level ancestry stratification."""
import sys
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
SUFFIX = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
SV_BASE = ROOT / "sv_extraction" / "decompose" / SUFFIX
OUT = ROOT / "analyses" / "imputability" / "output"
PANEL = ROOT / "resources" / "1000G" / "20130606_g1k_3202_samples_ped_population.txt"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def load_pop_map():
    df = pd.read_csv(PANEL, sep=r"\s+", engine="python")
    return dict(zip(df["SampleID"], df["Superpopulation"]))


def chrom_features(chrom, pop_map):
    sv_npz = SV_BASE / chrom / f"{chrom}_sv_genotypes.npz"
    sv_sam = SV_BASE / chrom / f"{chrom}_sv_samples.txt"
    if not sv_npz.exists() or not sv_sam.exists():
        print(f"[SKIP] {chrom}: missing", flush=True)
        return None

    g = np.load(sv_npz)
    arr = g[g.files[0]].astype(np.float32)
    samples = open(sv_sam).read().splitlines()
    pop_of_sample = np.array([pop_map.get(s, "OTHER") for s in samples])
    pop_idx = {pop: np.where(pop_of_sample == pop)[0] for pop in POPS}
    n_per_pop = {pop: len(pop_idx[pop]) for pop in POPS}
    print(f"[{chrom}] SVs={arr.shape[0]} haps={arr.shape[1]}  per-pop n={n_per_pop}",
          flush=True)
    if min(n_per_pop.values()) < 5:
        print(f"[FAIL] {chrom}: a super-pop has <5 samples ({n_per_pop})", flush=True)
        return None

    arr_m = np.where(np.isfinite(arr) & (arr >= 0), arr, np.nan)

    mx = np.nanmax(arr_m)
    is_diploid = mx > 1.05
    norm = 2.0 if is_diploid else 1.0
    print(f"[{chrom}] arr.max={mx:.2f}  treating as {'diploid' if is_diploid else 'haploid'}",
          flush=True)

    p = {}
    for pop in POPS:
        idx = pop_idx[pop]
        sub = arr_m[:, idx]
        n_obs = np.sum(~np.isnan(sub), axis=1)
        sum_d = np.nansum(sub, axis=1)
        af = np.where(n_obs > 0, sum_d / (norm * n_obs), np.nan)
        p[pop] = af

    af_stack = np.vstack([p[pop] for pop in POPS])
    af_skew = np.nanmax(af_stack, axis=0) - np.nanmin(af_stack, axis=0)

    all_idx = np.concatenate([pop_idx[pop] for pop in POPS])
    sub_all = arr_m[:, all_idx]
    n_all = np.sum(~np.isnan(sub_all), axis=1)
    total_af = np.where(n_all > 0,
                        np.nansum(sub_all, axis=1) / (norm * n_all), np.nan)

    pair_fst = {}
    for i, j in combinations(POPS, 2):
        pi = p[i]; pj = p[j]
        ni = max(int(norm * n_per_pop[i]), 2); nj = max(int(norm * n_per_pop[j]), 2)
        num = (pi - pj) ** 2 - pi * (1 - pi) / (ni - 1) \
                              - pj * (1 - pj) / (nj - 1)
        den = pi * (1 - pj) + pj * (1 - pi)
        fst = np.where(den > 0, num / den, np.nan)
        pair_fst[f"{i}_{j}"] = fst

    fst_stack = np.vstack(list(pair_fst.values()))
    max_fst = np.nanmax(fst_stack, axis=0)

    rows = []
    for i in range(arr.shape[0]):
        d = {"chrom": chrom, "sv_index": int(i), "total_af_sv": float(total_af[i])}
        for pop in POPS:
            d[f"af_{pop}_sv"] = float(p[pop][i])
        d["af_skew_sv"] = float(af_skew[i])
        d["max_FST_pair_sv"] = float(max_fst[i]) if np.isfinite(max_fst[i]) else np.nan
        rows.append(d)
    return pd.DataFrame(rows)


def main():
    pop_map = load_pop_map()
    print(f"[INFO] panel: {len(pop_map)} samples, super-pops: "
          f"{set(pop_map.values())}", flush=True)
    frames = []
    for chrom in CHROMS:
        df = chrom_features(chrom, pop_map)
        if df is not None:
            frames.append(df)
    full = pd.concat(frames, ignore_index=True)
    out_path = OUT / "sv_ancestry_features.tsv"
    full.to_csv(out_path, sep="\t", index=False)
    print(f"\n[OK] wrote {out_path}  ({len(full)} rows)", flush=True)
    print(full[[f"af_{pop}_sv" for pop in POPS] +
                ["af_skew_sv", "max_FST_pair_sv", "total_af_sv"]
                ].describe().to_string(), flush=True)


if __name__ == "__main__":
    main()
