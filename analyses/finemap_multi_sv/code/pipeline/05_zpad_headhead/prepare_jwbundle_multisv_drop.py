#!/usr/bin/env python3
"""Phase B: DROP-variant prepare for z=0 padding head-to-head."""
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


TABIX_MANIFEST = Path(
    CENSV_ROOT + "/"
    "analyses/sv_catalog/output/"
    "ref_compare_proxy_vs_dense/gwas_bgzip_tabix/_manifest.tsv"
)

BGZ_COLS = ["CHROM_INT", "POS", "SNP", "MAF", "N", "BETA", "SE", "PVALUE"]

TABIX = Path(
    CENSV_WORKSPACE + "/_shared/envs/workspace/v1/.pixi/envs/biotools/bin/tabix"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--pheno", required=True)
    p.add_argument("--study", default=None)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--tabix-manifest", type=Path, default=TABIX_MANIFEST)
    return p.parse_args()


def load_bundle(path: Path) -> dict:
    z = np.load(path, allow_pickle=False)
    return {k: z[k] for k in z.files}


def resolve_bgz_for_pheno(pheno: str, study: str | None, manifest_tsv: Path) -> tuple[Path, str]:
    mf = pd.read_csv(manifest_tsv, sep="\t")
    rows = mf[(mf["pheno_code"] == pheno) & (mf["study"] == study)] if study else mf[mf["pheno_code"] == pheno]
    if len(rows) != 1:
        raise ValueError(f"manifest lookup: pheno={pheno} study={study} -> {len(rows)} rows")
    return Path(rows.iloc[0]["output_bgz_path"]), str(rows.iloc[0]["study"])


def tabix_fetch(bgz_path: Path, chrom: str, start: int, end: int) -> pd.DataFrame:
    chrom_num = chrom.replace("chr", "")
    region = f"{chrom_num}:{start}-{end}"
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


def harmonize_gwas(snp_meta: dict, gwas_df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
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

    sys.stderr.write(f"[prepare-drop] bundle={args.bundle} pheno={args.pheno}\n")
    b = load_bundle(args.bundle)
    R_full = b["R_snp_snp"]
    B_full = b["B"]
    snp_meta = {k: b[k] for k in ("snp_id", "snp_chrom", "snp_pos", "snp_ref", "snp_alt", "snp_maf_panel", "snp_swap")}
    p_full, G = B_full.shape
    sys.stderr.write(f"[prepare-drop] bundle dims: p_full={p_full} G={G}\n")

    chrom = str(snp_meta["snp_chrom"][0])
    region_start = int(snp_meta["snp_pos"].min())
    region_end = int(snp_meta["snp_pos"].max())

    bgz, study = resolve_bgz_for_pheno(args.pheno, args.study, args.tabix_manifest)
    gwas_df = tabix_fetch(bgz, chrom, region_start, region_end)
    sys.stderr.write(f"[prepare-drop] gwas region rows: {len(gwas_df)}\n")

    z_snp_full, keep = harmonize_gwas(snp_meta, gwas_df)
    n_matched = int(keep.sum())
    sys.stderr.write(f"[prepare-drop] matched {n_matched}/{p_full}; DROPPING {p_full - n_matched} unmatched\n")

    if n_matched < 3:
        raise ValueError(f"Only {n_matched} SNPs harmonized — insufficient for SuSiE-RSS")

    idx = np.where(keep)[0]
    R = R_full[np.ix_(idx, idx)]
    B = B_full[idx, :]
    z_snp = z_snp_full[idx]
    p = len(idx)

    n_per_sv_after = (B != 0).sum(axis=0)
    n_per_sv_full = (B_full != 0).sum(axis=0)
    n_svs_lost_majority = int(((n_per_sv_after / np.maximum(n_per_sv_full, 1)) < 0.5).sum())
    n_svs_lost_all = int((n_per_sv_after == 0).sum())
    sys.stderr.write(
        f"[prepare-drop] p after drop: {p}; SVs lost majority support: {n_svs_lost_majority}; "
        f"SVs lost ALL support: {n_svs_lost_all}\n"
    )

    sv_keep_mask = (n_per_sv_after >= 3)
    n_svs_dropped_no_support = int((~sv_keep_mask).sum())
    if n_svs_dropped_no_support > 0:
        sys.stderr.write(
            f"[prepare-drop] dropping {n_svs_dropped_no_support} SVs with <3 panel-matched support after drop\n"
        )
        B = B[:, sv_keep_mask]
        sv_locus_kept    = b["sv_locus"][sv_keep_mask]
        sv_index_kept    = b["sv_index"][sv_keep_mask]
        sv_category_kept = b["sv_category"][sv_keep_mask]
        sv_strategy_kept = b["sv_strategy"][sv_keep_mask]
        sv_arm_kept      = b["sv_arm"][sv_keep_mask]
        sv_cv_r2_kept    = b["sv_cv_r2"][sv_keep_mask]
    else:
        sv_locus_kept    = b["sv_locus"]
        sv_index_kept    = b["sv_index"]
        sv_category_kept = b["sv_category"]
        sv_strategy_kept = b["sv_strategy"]
        sv_arm_kept      = b["sv_arm"]
        sv_cv_r2_kept    = b["sv_cv_r2"]
    G_kept = B.shape[1]
    sys.stderr.write(f"[prepare-drop] G after support-drop: {G_kept} (from {G})\n")

    sys.stderr.write(f"[prepare-drop] building canonical R_aug (p+G={p+G})\n")
    R_aug, D = build_canonical_raug(R, B)
    z_aug = build_z_aug(z_snp, B, D)

    R_aug_safe, info = psd_safeguard(R_aug)
    sys.stderr.write(f"[prepare-drop] PSD safeguard: eig_min_pre={info['eig_min_pre']:.3e} "
                     f"post={info['eig_min_post']:.3e} tikhonov={info['tikhonov_applied']}\n")

    snp_meta_dropped = {k: v[idx] if hasattr(v, "__getitem__") and len(v) == p_full else v
                        for k, v in snp_meta.items()}

    out_path = args.out_dir / f"{args.bundle.stem.replace('_bundle', '')}__{args.pheno}__jw_drop.npz"
    atomic_write_npz(
        out_path,
        R_aug=R_aug_safe,
        z_aug=z_aug,
        D_normalizer=D,
        z_snp=z_snp,
        sv_locus=sv_locus_kept,
        sv_index=sv_index_kept,
        sv_category=sv_category_kept,
        sv_strategy=sv_strategy_kept,
        sv_arm=sv_arm_kept,
        sv_cv_r2=sv_cv_r2_kept,
        snp_id=snp_meta_dropped["snp_id"],
        snp_chrom=snp_meta_dropped["snp_chrom"],
        snp_pos=snp_meta_dropped["snp_pos"],
    )

    info_dict = {
        "variant": "drop_unmatched",
        "p_snps_full": int(p_full),
        "p_snps_after_drop": int(p),
        "G_svs_full": int(G),
        "G_svs_kept": int(G_kept),
        "n_matched": int(n_matched),
        "n_svs_lost_majority_support": int(n_svs_lost_majority),
        "n_svs_lost_all_support": int(n_svs_lost_all),
        "n_svs_dropped_no_support_after_drop": int(n_svs_dropped_no_support),
        "n_gwas_region": int(len(gwas_df)),
        "eig_min_pre": float(info["eig_min_pre"]),
        "eig_min_post": float(info["eig_min_post"]),
        "tikhonov_applied": bool(info["tikhonov_applied"]),
        "pheno": args.pheno,
        "study": study,
        "bundle": str(args.bundle),
        "bgz": str(bgz),
        "r_aug_mode": "canonical_diag1_normalized",
    }
    atomic_write_text(args.out_dir / f"{args.bundle.stem.replace('_bundle', '')}__{args.pheno}__info_drop.json",
                      json.dumps(info_dict, indent=2) + "\n")
    sys.stderr.write(f"[prepare-drop] DONE wrote {out_path}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
