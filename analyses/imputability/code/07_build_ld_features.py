#!/usr/bin/env python3
"""Build SNP-LD features per SV — the missing ~⅓ of CV R² variance."""
import sys
import os
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

warnings.filterwarnings("ignore")

ROOT = Path(CENSV_ROOT + "")
SUFFIX = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
SV_BASE = ROOT / "sv_extraction" / "decompose" / SUFFIX
SNP_BASE = ROOT / "snp_extraction"
OUT = ROOT / "analyses" / "imputability" / "output"

STRATEGY_DIR_MAP = {
    "edge500": "edge500", "edge1000": "edge1000", "edge2000": "edge2000",
    "fksize100k": "fksize100000", "fksize200k": "fksize200000",
    "fksize500k": "fksize500000",
    "midpoint500": "midpoint500", "midpoint1000": "midpoint1000",
}
STRATEGIES = list(STRATEGY_DIR_MAP.keys())
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

MAX_PAIR_SNPS = 5000
MAX_EIGEN_SNPS = 2000
RNG = np.random.default_rng(42)


def safe_corr(A, b):
    A = A.astype(np.float32)
    b = b.astype(np.float32)
    A_mean = A.mean(axis=1, keepdims=True)
    b_mean = b.mean()
    A_c = A - A_mean
    b_c = b - b_mean
    A_norm = np.sqrt((A_c ** 2).sum(axis=1))
    b_norm = np.sqrt((b_c ** 2).sum())
    denom = A_norm * b_norm
    num = A_c @ b_c
    r = np.where(denom > 0, num / denom, 0.0)
    return r


def pairwise_r2_mean(X, max_snps=MAX_PAIR_SNPS):
    n = X.shape[0]
    if n < 2:
        return np.nan
    if n > max_snps:
        idx = RNG.choice(n, max_snps, replace=False)
        X = X[idx]
    Xc = X - X.mean(axis=1, keepdims=True)
    norm = np.sqrt((Xc ** 2).sum(axis=1))
    norm[norm == 0] = 1
    Xn = Xc / norm[:, None]
    R = Xn @ Xn.T
    np.fill_diagonal(R, 0)
    R2 = R ** 2
    n_pairs = R2.shape[0] * (R2.shape[0] - 1)
    return float(R2.sum() / n_pairs) if n_pairs > 0 else np.nan


def n_indep_eigen(X, max_snps=MAX_EIGEN_SNPS):
    n = X.shape[0]
    if n < 2:
        return float(n)
    if n > max_snps:
        idx = RNG.choice(n, max_snps, replace=False)
        X = X[idx]
    Xc = X - X.mean(axis=1, keepdims=True)
    norm = np.sqrt((Xc ** 2).sum(axis=1))
    norm[norm == 0] = 1
    Xn = Xc / norm[:, None]
    R = Xn @ Xn.T
    eigvals = np.linalg.eigvalsh(R)
    eigvals = np.clip(eigvals, 0, None)
    s = eigvals.sum()
    if s <= 0:
        return float(n)
    eigvals = eigvals / s * R.shape[0]
    return float(R.shape[0] / (eigvals ** 2).sum() * R.shape[0])


def process_chrom_strategy(chrom, strategy):
    snp_dir = SNP_BASE / STRATEGY_DIR_MAP[strategy]
    snp_npz = snp_dir / f"{chrom}_snp_genotypes.npz"
    snp_var = snp_dir / f"{chrom}_snp_variants.tsv"
    snp_samp = snp_dir / f"{chrom}_snp_samples.txt"
    sv_npz = SV_BASE / chrom / f"{chrom}_sv_genotypes.npz"
    sv_var = SV_BASE / chrom / f"{chrom}_sv_variants.tsv"
    sv_samp = SV_BASE / chrom / f"{chrom}_sv_samples.txt"
    for p in [snp_npz, snp_var, snp_samp, sv_npz, sv_var, sv_samp]:
        if not p.exists():
            return None

    snp_g = np.load(snp_npz)
    snp_arr = snp_g[snp_g.files[0]]
    sv_g = np.load(sv_npz)
    sv_arr = sv_g[sv_g.files[0]]
    snp_pos = pd.read_csv(snp_var, sep="\t")["pos"].values
    snp_samples = open(snp_samp).read().splitlines()
    sv_samples = open(sv_samp).read().splitlines()

    common = [s for s in sv_samples if s in set(snp_samples)]
    if len(common) < 50:
        print(f"[WARN] {chrom} × {strategy}: only {len(common)} common haps "
              f"(SV={len(sv_samples)}, SNP={len(snp_samples)})", flush=True)
        return None
    sv_idx_map = {s: i for i, s in enumerate(sv_samples)}
    snp_idx_map = {s: i for i, s in enumerate(snp_samples)}
    sv_keep = [sv_idx_map[s] for s in common]
    snp_keep = [snp_idx_map[s] for s in common]
    sv_arr = sv_arr[:, sv_keep]
    snp_arr = snp_arr[:, snp_keep]

    sv_var_df = pd.read_csv(sv_var, sep="\t")
    if "pos" in sv_var_df.columns:
        sv_pos_arr = sv_var_df["pos"].values
    elif "svpos" in sv_var_df.columns:
        sv_pos_arr = sv_var_df["svpos"].values
    else:
        sv_pos_arr = sv_var_df.get("start", pd.Series([np.nan]*len(sv_var_df))).values

    mean_pair = pairwise_r2_mean(snp_arr.astype(np.float32))
    n_indep = n_indep_eigen(snp_arr.astype(np.float32))

    rows = []
    for i in range(sv_arr.shape[0]):
        sv_i = sv_arr[i].astype(np.float32)
        if np.std(sv_i) == 0:
            rows.append({
                "chrom": chrom, "strategy": strategy, "sv_index": i,
                "top_tag_r2": np.nan, "n_high_ld_r2_05": np.nan,
                "max_r2_within_50kb": np.nan,
                "max_r2_within_500kb": np.nan,
                "max_r2_within_5Mb": np.nan,
                "mean_pair_r2_window": mean_pair,
                "n_indep_snps_window": n_indep,
                "n_window_snps_actual": int(snp_arr.shape[0]),
            })
            continue
        r = safe_corr(snp_arr.astype(np.float32), sv_i)
        r2 = r ** 2
        top = float(r2.max())
        n_high = int((r2 >= 0.5).sum())
        sv_pos = sv_pos_arr[i] if i < len(sv_pos_arr) else np.nan
        max_50k = max_500k = max_5Mb = np.nan
        if not np.isnan(sv_pos):
            d = np.abs(snp_pos - sv_pos)
            for cap, key in [(50_000, "50kb"), (500_000, "500kb"), (5_000_000, "5Mb")]:
                mask = d <= cap
                if mask.any():
                    val = float(r2[mask].max())
                    if key == "50kb": max_50k = val
                    elif key == "500kb": max_500k = val
                    else: max_5Mb = val
        rows.append({
            "chrom": chrom, "strategy": strategy, "sv_index": int(i),
            "top_tag_r2": top, "n_high_ld_r2_05": n_high,
            "max_r2_within_50kb": max_50k,
            "max_r2_within_500kb": max_500k,
            "max_r2_within_5Mb": max_5Mb,
            "mean_pair_r2_window": mean_pair,
            "n_indep_snps_window": n_indep,
            "n_window_snps_actual": int(snp_arr.shape[0]),
        })
    return pd.DataFrame(rows)


def main():
    chrom_arg = os.environ.get("CHROM_ARG")
    chroms = [chrom_arg] if chrom_arg else CHROMS

    all_rows = []
    for chrom in chroms:
        for strategy in STRATEGIES:
            print(f"[INFO] {chrom} × {strategy}...", flush=True)
            try:
                df = process_chrom_strategy(chrom, strategy)
            except Exception as e:
                print(f"[ERR] {chrom} × {strategy}: {type(e).__name__}: {e}",
                      flush=True)
                continue
            if df is None or len(df) == 0:
                print(f"[SKIP] {chrom} × {strategy}: missing inputs", flush=True)
                continue
            all_rows.append(df)
            print(f"  → {len(df)} rows  top_tag_r2 mean={df['top_tag_r2'].mean():.3f}",
                  flush=True)

    if not all_rows:
        print("[FAIL] no rows produced", flush=True)
        sys.exit(1)
    full = pd.concat(all_rows, ignore_index=True)
    suffix = f"_{chrom_arg}" if chrom_arg else ""
    chunk_dir = OUT / "ld_chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    out_long = chunk_dir / f"ld_features_per_sv_strategy{suffix}.tsv"
    full.to_csv(out_long, sep="\t", index=False)
    print(f"\n[OK] wrote {out_long} ({len(full)} rows)", flush=True)


if __name__ == "__main__":
    main()
