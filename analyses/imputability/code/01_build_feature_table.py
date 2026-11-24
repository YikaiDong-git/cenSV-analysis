#!/usr/bin/env python3
"""Build a per-SV feature matrix for imputability analysis."""

import os
import sys
import ast
import numpy as np
import pandas as pd
from pathlib import Path
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
SUFFIX = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
STRATEGIES = ["edge500", "edge1000", "edge2000",
              "fksize100k", "fksize200k", "fksize500k",
              "midpoint500", "midpoint1000"]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
ACROCENTRIC = {"chr13", "chr14", "chr15", "chr21", "chr22"}

OUT_DIR = ROOT / "analyses" / "imputability" / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_gt_counts(s):
    out = {}
    if not isinstance(s, str):
        return out
    for tok in s.split("/"):
        if ":" not in tok:
            continue
        k, v = tok.split(":", 1)
        out[int(k)] = int(v)
    return out


def parse_dist_list(s):
    if not isinstance(s, str) or not s.strip():
        return []
    try:
        return [int(x) for x in ast.literal_eval(s)]
    except (ValueError, SyntaxError):
        return []


def load_active_hor():
    bed = ROOT / "resources" / "chm13v2.0_censat_v2.1.active.bed"
    rows = []
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            if len(f) < 4:
                continue
            rows.append({"chrom": f[0], "start": int(f[1]), "end": int(f[2]), "name": f[3]})
    df = pd.DataFrame(rows)
    agg = df.groupby("chrom").agg(
        active_HOR_span_bp=("end", lambda e: int((e - df.loc[e.index, "start"]).sum())),
        n_HOR_intervals=("name", "count"),
    ).reset_index()
    return agg


def load_pericen_lengths():
    out = {}
    for size in [100000, 200000, 500000, 1000000, 2000000]:
        bed = ROOT / "resources" / f"peri_centromeric_regions_grch38_{size}.bed"
        if not bed.exists():
            continue
        rows = []
        with open(bed) as fh:
            for line in fh:
                f = line.rstrip().split("\t")
                if len(f) < 3:
                    continue
                rows.append((f[0], int(f[1]), int(f[2])))
        df = pd.DataFrame(rows, columns=["chrom", "s", "e"])
        df["len"] = df["e"] - df["s"]
        out[f"pericen_len_{size}"] = df.groupby("chrom")["len"].sum().to_dict()
    return out


def load_per_chrom_best8(chrom):
    base = (ROOT / "results" / "imputation" / "decompose" / "cosine_nj"
            / f"best8_{SUFFIX}" / chrom)
    imp_path = base / f"{chrom}_imputation_summary.tsv"
    sel_path = base / f"{chrom}_best_strategy_selection.tsv"
    if not imp_path.exists() or not sel_path.exists():
        return None
    imp = pd.read_csv(imp_path, sep="\t")
    sel = pd.read_csv(sel_path, sep="\t")
    if "chrom" not in imp.columns or imp["chrom"].isna().all():
        imp["chrom"] = chrom
    df = imp.merge(sel, on="sv_index", how="left", suffixes=("", "_sel"))
    return df


def load_per_chrom_per_strategy(chrom, strategy):
    base = (ROOT / "results" / "elastic_net" / "decompose" / "cosine_nj"
            / f"{strategy}_{SUFFIX}" / chrom / "sv_metrics.tsv")
    if not base.exists():
        return {}
    df = pd.read_csv(base, sep="\t")
    out = {}
    for _, r in df.iterrows():
        idx = int(r["sv_index"])
        gt = parse_gt_counts(r.get("gt_counts", ""))
        dists = parse_dist_list(r.get("selected_sv_snp_dist", ""))
        n_carriers = gt.get(1, 0) + gt.get(2, 0)
        n_hom_alt = gt.get(2, 0)
        n_total = gt.get(0, 0) + gt.get(1, 0) + gt.get(2, 0)
        out[idx] = {
            f"l1_ratio_{strategy}": r.get("l1_ratio", np.nan),
            f"alpha_{strategy}": r.get("alpha", np.nan),
            f"n_nonzero_{strategy}": r.get("n_nonzero_coef", np.nan),
            f"snp_dist_mean_{strategy}": np.mean(dists) if dists else np.nan,
            f"snp_dist_min_{strategy}": np.min(dists) if dists else np.nan,
            f"snp_dist_max_{strategy}": np.max(dists) if dists else np.nan,
            f"n_window_snps_{strategy}": len(dists) if dists else np.nan,
            "n_carriers": n_carriers,
            "n_hom_alt": n_hom_alt,
            "n_total_haps": n_total,
        }
    return out


def main():
    hor = load_active_hor()
    pericen = load_pericen_lengths()

    all_chrom_dfs = []
    for chrom in CHROMS:
        print(f"[INFO] {chrom}: loading best8...", flush=True)
        df = load_per_chrom_best8(chrom)
        if df is None:
            print(f"[WARN] {chrom}: missing best8 inputs, skipping", flush=True)
            continue

        strat_cache = {}
        for strat in STRATEGIES:
            strat_cache[strat] = load_per_chrom_per_strategy(chrom, strat)

        rows_extra = []
        for _, r in df.iterrows():
            idx = int(r["sv_index"])
            best_strat = r.get("best_strategy_by_cvr2", "")
            extras = {}
            if best_strat in strat_cache and idx in strat_cache[best_strat]:
                src = strat_cache[best_strat][idx]
                extras["best_l1_ratio"] = src.get(f"l1_ratio_{best_strat}", np.nan)
                extras["best_alpha"] = src.get(f"alpha_{best_strat}", np.nan)
                extras["best_n_nonzero"] = src.get(f"n_nonzero_{best_strat}", np.nan)
                extras["best_snp_dist_mean"] = src.get(f"snp_dist_mean_{best_strat}", np.nan)
                extras["best_snp_dist_min"] = src.get(f"snp_dist_min_{best_strat}", np.nan)
                extras["best_snp_dist_max"] = src.get(f"snp_dist_max_{best_strat}", np.nan)
                extras["n_window_snps"] = src.get(f"n_window_snps_{best_strat}", np.nan)
                extras["n_carriers"] = src.get("n_carriers", np.nan)
                extras["n_hom_alt"] = src.get("n_hom_alt", np.nan)
                extras["n_total_haps"] = src.get("n_total_haps", np.nan)
            else:
                for strat in STRATEGIES:
                    if idx in strat_cache.get(strat, {}):
                        src = strat_cache[strat][idx]
                        extras["n_carriers"] = src.get("n_carriers", np.nan)
                        extras["n_hom_alt"] = src.get("n_hom_alt", np.nan)
                        extras["n_total_haps"] = src.get("n_total_haps", np.nan)
                        break
            rows_extra.append(extras)

        extras_df = pd.DataFrame(rows_extra)
        df = pd.concat([df.reset_index(drop=True), extras_df.reset_index(drop=True)], axis=1)

        cvr_cols = [f"cv_r2_{s}" for s in STRATEGIES]
        existing = [c for c in cvr_cols if c in df.columns]
        if existing:
            mat = df[existing].values
            df["cv_r2_max"] = np.nanmax(mat, axis=1)
            df["cv_r2_min"] = np.nanmin(mat, axis=1)
            df["cv_r2_mean"] = np.nanmean(mat, axis=1)
            df["cv_r2_var"] = np.nanvar(mat, axis=1)
            df["cv_r2_range"] = df["cv_r2_max"] - df["cv_r2_min"]

        h = hor[hor["chrom"] == chrom]
        if len(h):
            df["active_HOR_span_bp"] = int(h["active_HOR_span_bp"].iloc[0])
            df["n_HOR_intervals"] = int(h["n_HOR_intervals"].iloc[0])
        else:
            df["active_HOR_span_bp"] = np.nan
            df["n_HOR_intervals"] = np.nan
        for k, mp in pericen.items():
            df[k] = mp.get(chrom, np.nan)
        df["acrocentric"] = (chrom in ACROCENTRIC)

        all_chrom_dfs.append(df)

    if not all_chrom_dfs:
        print("[FAIL] no chromosomes loaded", flush=True)
        sys.exit(1)

    full = pd.concat(all_chrom_dfs, ignore_index=True)
    print(f"[INFO] total rows: {len(full)} across {full['chrom'].nunique()} chroms",
          flush=True)

    full["log_svlen"] = np.log10(full["svlen"].clip(lower=1))
    full["log_active_HOR_span"] = np.log10(full["active_HOR_span_bp"].clip(lower=1))
    full["maf"] = full["svfreq"].clip(upper=0.5)
    full["log_maf"] = np.log10(full["maf"].clip(lower=1e-3))
    full["n_train_haps"] = np.where(full["chrom"] == "chr13", 190, 196)
    full["log_n_window_snps"] = np.log10(full["n_window_snps"].clip(lower=1))

    for t in ["INS", "DEL", "INV", "DUP"]:
        full[f"svtype_{t}"] = (full["svtype"] == t).astype(int)
    for c in ["monophyletic", "residual", "non_decomposed", "original_decomposed"]:
        full[f"cat_{c}"] = (full["sv_category"] == c).astype(int)

    full["y_cv_r2"] = full["cv_r2"]
    full["y_pass_01"] = (full["cv_r2"] >= 0.1).astype(int)
    full["y_pass_03"] = (full["cv_r2"] >= 0.3).astype(int)

    out_path = OUT_DIR / "per_sv_features.tsv"
    full.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {out_path} ({len(full)} rows × {full.shape[1]} cols)",
          flush=True)
    print("\nSanity summary:", flush=True)
    print(full[["chrom", "y_cv_r2"]].groupby("chrom").agg(["count", "mean"]).to_string(),
          flush=True)
    print(f"\nOverall pass@0.1: {full['y_pass_01'].mean():.3f}", flush=True)
    print(f"Overall pass@0.3: {full['y_pass_03'].mean():.3f}", flush=True)


if __name__ == "__main__":
    main()
