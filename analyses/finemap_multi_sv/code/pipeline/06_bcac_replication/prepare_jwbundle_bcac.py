#!/usr/bin/env python3
"""Phase C: prepare jwbundle for chr1_super1 × BCAC ER-neg replication."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent.parent / "lib"))
from raug import build_canonical_raug, build_z_aug, psd_safeguard


BASE = Path(
    CENSV_ROOT + "/"
    "analyses/finemap_multi_sv"
)
DEFAULT_BCAC = BASE / "inputs/bcac_2017/processed_b38/bcac_erneg_chr1_b38.tsv.gz"
DEFAULT_BUNDLE = (
    BASE / "output/full_restart/bundles/chr1/"
    "chr1_super1_118296444_148907004_bundle.npz"
)
TABIX = Path(
    CENSV_WORKSPACE + "/_shared/envs/workspace/v1/.pixi/envs/biotools/bin/tabix"
)

BGZ_COLS = ["CHROM_INT", "POS", "SNP", "MAF", "N", "BETA", "SE", "PVALUE"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    p.add_argument("--bcac-bgz", type=Path, default=DEFAULT_BCAC)
    p.add_argument("--out-dir", type=Path,
                   default=BASE / "output/full_restart/bcac_replication")
    p.add_argument("--pheno-tag", default="bcac_erneg")
    return p.parse_args()


def load_bundle(path: Path) -> dict:
    z = np.load(path, allow_pickle=False)
    return {k: z[k] for k in z.files}


def tabix_fetch(bgz_path: Path, chrom_int: int, start: int, end: int) -> pd.DataFrame:
    region = f"{chrom_int}:{start}-{end}"
    cmd = [str(TABIX), str(bgz_path), region]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    lines = [ln for ln in proc.stdout.strip().split("\n") if ln and not ln.startswith("#")]
    if not lines:
        raise ValueError(f"tabix returned 0 lines for {bgz_path} {region}")
    rows = [ln.split("\t") for ln in lines]
    df = pd.DataFrame(rows, columns=BGZ_COLS)
    df["POS"] = df["POS"].astype(np.int64)
    for col in ("MAF", "N", "BETA", "SE", "PVALUE"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _normalize_gwas_snp_id(s: str) -> str:
    s = s.replace("chr", "")
    parts = s.split("_")
    if len(parts) == 4:
        return f"{parts[0]}:{parts[1]}:{parts[2]}:{parts[3]}"
    return s


def harmonize_gwas(snp_meta: dict, gwas_df: pd.DataFrame):
    panel_ids = snp_meta["snp_id"].astype(str)
    p = len(panel_ids)
    z = np.full(p, np.nan, dtype=np.float64)
    keep = np.zeros(p, dtype=bool)

    gwas_idx: dict[str, float] = {}
    for snp_id, beta, se in zip(gwas_df["SNP"], gwas_df["BETA"], gwas_df["SE"]):
        if pd.isna(beta) or pd.isna(se) or se <= 0:
            continue
        canon = _normalize_gwas_snp_id(str(snp_id))
        gwas_idx[canon] = float(beta) / float(se)

    for i in range(p):
        pid = str(panel_ids[i])
        if pid in gwas_idx:
            z[i] = gwas_idx[pid]
            keep[i] = True
            continue
        parts = pid.split(":")
        if len(parts) == 4:
            rev = f"{parts[0]}:{parts[1]}:{parts[3]}:{parts[2]}"
            if rev in gwas_idx:
                z[i] = -gwas_idx[rev]
                keep[i] = True
    z[~keep] = 0.0
    return z, keep


def atomic_write_npz(path: Path, **arrays) -> None:
    tmp = path.parent / f"{path.stem}.tmp.{os.getpid()}.npz"
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, path)


def atomic_write_text(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    sys.stderr.write(f"[bcac-prepare] bundle={args.bundle} pheno={args.pheno_tag}\n")
    b = load_bundle(args.bundle)
    R = b["R_snp_snp"]
    B = b["B"]
    snp_meta = {k: b[k] for k in ("snp_id", "snp_chrom", "snp_pos", "snp_ref", "snp_alt", "snp_maf_panel", "snp_swap")}
    p, G = B.shape

    chrom_int = int(str(snp_meta["snp_chrom"][0]).replace("chr", ""))
    region_start = int(snp_meta["snp_pos"].min())
    region_end = int(snp_meta["snp_pos"].max())

    gwas_df = tabix_fetch(args.bcac_bgz, chrom_int, region_start, region_end)
    sys.stderr.write(f"[bcac-prepare] BCAC ER-neg region rows: {len(gwas_df)}\n")

    z_snp, keep = harmonize_gwas(snp_meta, gwas_df)
    n_matched = int(keep.sum())
    sys.stderr.write(f"[bcac-prepare] harmonized: matched {n_matched}/{p}\n")

    if n_matched < 3:
        raise ValueError(f"Only {n_matched} SNPs harmonized — insufficient")

    R_aug, D = build_canonical_raug(R, B)
    z_aug = build_z_aug(z_snp, B, D)
    R_aug_safe, info = psd_safeguard(R_aug)

    region_id = args.bundle.stem.replace("_bundle", "")
    out_path = args.out_dir / f"{region_id}__{args.pheno_tag}__jw.npz"
    atomic_write_npz(
        out_path,
        R_aug=R_aug_safe,
        z_aug=z_aug,
        D_normalizer=D,
        z_snp=z_snp,
        snp_matched=keep.astype(np.int32),
        sv_locus=b["sv_locus"],
        sv_index=b["sv_index"],
        sv_category=b["sv_category"],
        sv_strategy=b["sv_strategy"],
        sv_arm=b["sv_arm"],
        sv_cv_r2=b["sv_cv_r2"],
        snp_id=snp_meta["snp_id"],
        snp_chrom=snp_meta["snp_chrom"],
        snp_pos=snp_meta["snp_pos"],
    )

    info_dict = {
        "p_snps": int(p),
        "G_svs": int(G),
        "n_matched": int(n_matched),
        "n_gwas_region": int(len(gwas_df)),
        "eig_min_pre": float(info["eig_min_pre"]),
        "eig_min_post": float(info["eig_min_post"]),
        "pheno": args.pheno_tag,
        "study": "BCAC_Michailidou2017_ERneg_GCST004988",
        "bundle": str(args.bundle),
        "bgz": str(args.bcac_bgz),
        "r_aug_mode": "canonical_diag1_normalized",
        "v2_finemap_recipe": "L=5, max_iter=500",
    }
    atomic_write_text(
        args.out_dir / f"{region_id}__{args.pheno_tag}__info.json",
        json.dumps(info_dict, indent=2) + "\n"
    )
    sys.stderr.write(f"[bcac-prepare] DONE {out_path}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
