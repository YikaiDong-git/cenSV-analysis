#!/usr/bin/env python
"""build_rescue_manifest.py — Stratified rescue manifest builder for R_AUG_TOO_LARGE stubs."""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

RNG_SEED = 42

PROJ = Path(CENSV_ROOT + "/analyses/sv_catalog")
SWEEP_ROOT = PROJ / "output/ref_compare_proxy_vs_dense/sweep_6519_locked_v3"
MASTER = SWEEP_ROOT / "manifest/MASTER_MANIFEST_6519.tsv"
MASTER_COLS = [
    "pair_id", "sv_locus", "sv_id_builder", "sv_chrom", "svpos",
    "strategy", "arm", "win_start", "win_end",
    "bundle_dir", "sv_category", "cv_r2",
]


def collect_stubs(pheno: str) -> pd.DataFrame:
    diag_root = SWEEP_ROOT / "diag" / pheno
    rows = []
    for chrom_dir in sorted(diag_root.iterdir()):
        if not chrom_dir.is_dir():
            continue
        for diag in sorted(chrom_dir.glob("*__diag.tsv")):
            with open(diag) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                data = fh.readline().rstrip("\n").split("\t")
            if "status" not in header:
                continue
            d = dict(zip(header, data))
            if d.get("status") != "R_AUG_TOO_LARGE":
                continue
            rows.append(dict(
                sv_locus=d["sv_locus"],
                arm=d["arm"],
                sv_chrom=d["chrom"],
                n_panel_snp=int(d["n_panel_snp"]),
            ))
    df = pd.DataFrame(rows)
    print(f"[stubs] {pheno}: {len(df)} R_AUG_TOO_LARGE rows", file=sys.stderr)
    return df


def stratify(stubs: pd.DataFrame, master: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    joined = stubs.merge(
        master, on=["sv_locus", "arm", "sv_chrom"], how="inner", validate="one_to_one",
    )
    assert len(joined) == len(stubs), \
        f"join lost rows: {len(stubs)} stubs -> {len(joined)} joined; master schema drift"

    anchor_mask = (joined["sv_chrom"] == "chr20") & joined["sv_locus"].str.contains("_sv212_")
    anchor = joined[anchor_mask].assign(rescue_group="anchor_sv212")

    WIDEST_CAP = 100
    rest = joined[~anchor_mask].copy()
    decile_cut = rest["n_panel_snp"].quantile(0.90)
    widest_full = rest[rest["n_panel_snp"] >= decile_cut]
    lower_p = rest[rest["n_panel_snp"] < decile_cut]
    if len(widest_full) > WIDEST_CAP:
        widest_idx = rng.choice(len(widest_full), size=WIDEST_CAP, replace=False)
        widest = widest_full.iloc[sorted(widest_idx)].assign(rescue_group="widest_decile")
    else:
        widest = widest_full.assign(rescue_group="widest_decile")

    n_random = min(50, len(lower_p))
    if n_random > 0:
        random_idx = rng.choice(len(lower_p), size=n_random, replace=False)
        random50 = lower_p.iloc[sorted(random_idx)].assign(rescue_group="random50")
    else:
        random50 = lower_p.iloc[:0].assign(rescue_group="random50")

    out = pd.concat([anchor, widest, random50], ignore_index=True)
    out = out.sort_values(["sv_chrom", "pair_id"]).reset_index(drop=True)
    print(
        f"[strat] anchor={len(anchor)} widest_decile={len(widest)} (cut p>={int(decile_cut)}) "
        f"random50={len(random50)} total={len(out)}",
        file=sys.stderr,
    )
    return out


def write_atomic(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    cols = MASTER_COLS + ["n_panel_snp", "rescue_group"]
    df[cols].to_csv(tmp, sep="\t", index=False)
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    print(f"[write] {path}  rows={len(df)}", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phenos", nargs="+", required=True,
                    help="phenotype identifiers (must match diag/<pheno>/ dir names)")
    ap.add_argument("--out-root", type=Path, default=SWEEP_ROOT / "manifest",
                    help="output dir for MISSING_RESCUE_<pheno>.tsv")
    args = ap.parse_args()

    rng = np.random.default_rng(RNG_SEED)
    master = pd.read_csv(MASTER, sep="\t")
    assert list(master.columns) == MASTER_COLS, \
        f"master schema drift: {master.columns.tolist()}"
    print(f"[master] {len(master)} rows", file=sys.stderr)

    for pheno in args.phenos:
        stubs = collect_stubs(pheno)
        rescue = stratify(stubs, master, rng)
        out_path = args.out_root / f"MISSING_RESCUE_{pheno}.tsv"
        write_atomic(rescue, out_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
