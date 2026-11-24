#!/usr/bin/env python3
"""Phase E: cTWAS-faithful prepare_jwbundle."""
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
from raug import build_canonical_raug, psd_safeguard


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
    p.add_argument("--min-matched-support", type=int, default=3,
                   help="minimum #matched β-support SNPs per SV (else z_SV = 0)")
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

    n_gwas = len(gwas_df)

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
    return z, keep, n_gwas


def build_z_aug_ctwas(z_snp: np.ndarray, B: np.ndarray, R: np.ndarray,
                      keep: np.ndarray, min_support: int = 3) -> tuple[np.ndarray, dict]:
    p, G = B.shape
    z_aug = np.empty(p + G, dtype=np.float64)
    z_aug[:p] = z_snp

    n_dropped = 0
    n_per_sv_intersect = np.zeros(G, dtype=np.int64)

    for j in range(G):
        support = np.where(B[:, j] != 0)[0]
        intersect = support[keep[support]]
        n_per_sv_intersect[j] = len(intersect)
        if len(intersect) < min_support:
            z_aug[p + j] = 0.0
            n_dropped += 1
            continue
        beta_int = B[intersect, j]
        z_int = z_snp[intersect]
        R_int = R[np.ix_(intersect, intersect)]
        bRb = float(beta_int @ R_int @ beta_int)
        if bRb <= 0:
            z_aug[p + j] = 0.0
            n_dropped += 1
            continue
        D_int = float(np.sqrt(bRb))
        z_aug[p + j] = float(beta_int @ z_int) / D_int

    info = {
        "n_svs_total": int(G),
        "n_svs_z_aug_dropped": int(n_dropped),
        "n_per_sv_matched_intersect_mean": float(n_per_sv_intersect.mean()),
        "n_per_sv_matched_intersect_min": int(n_per_sv_intersect.min()),
        "n_per_sv_matched_intersect_max": int(n_per_sv_intersect.max()),
    }
    return z_aug, info


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

    sys.stderr.write(f"[prepare-ctwas] bundle={args.bundle} pheno={args.pheno}\n")
    b = load_bundle(args.bundle)
    R = b["R_snp_snp"]
    B = b["B"]
    snp_meta = {k: b[k] for k in ("snp_id", "snp_chrom", "snp_pos", "snp_ref", "snp_alt", "snp_maf_panel", "snp_swap")}
    p, G = B.shape

    chrom = str(snp_meta["snp_chrom"][0])
    region_start = int(snp_meta["snp_pos"].min())
    region_end = int(snp_meta["snp_pos"].max())

    bgz, study = resolve_bgz_for_pheno(args.pheno, args.study, args.tabix_manifest)
    gwas_df = tabix_fetch(bgz, chrom, region_start, region_end)
    n_gwas_region = len(gwas_df)
    sys.stderr.write(f"[prepare-ctwas] gwas region rows: {n_gwas_region}\n")

    N_GWAS_MEAN = float(gwas_df["N"].dropna().mean()) if "N" in gwas_df.columns else float("nan")

    z_snp, keep, _ = harmonize_gwas(snp_meta, gwas_df)
    n_matched = int(keep.sum())
    sys.stderr.write(f"[prepare-ctwas] panel SNPs matched: {n_matched}/{p}\n")

    if n_matched < 3:
        raise ValueError(f"Only {n_matched} SNPs harmonized — insufficient")

    R_aug, D_full = build_canonical_raug(R, B)
    R_aug_safe, info_psd = psd_safeguard(R_aug)

    z_aug, info_zaug = build_z_aug_ctwas(z_snp, B, R, keep, min_support=args.min_matched_support)

    out_path = args.out_dir / f"{args.bundle.stem.replace('_bundle', '')}__{args.pheno}__jw_ctwas.npz"
    atomic_write_npz(
        out_path,
        R_aug=R_aug_safe,
        z_aug=z_aug,
        D_normalizer=D_full,
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
        "variant": "ctwas_faithful_intermediate",
        "p_snps": int(p),
        "G_svs": int(G),
        "n_matched": int(n_matched),
        "n_gwas_region": int(n_gwas_region),
        "n_gwas_mean": N_GWAS_MEAN,
        "eig_min_pre": float(info_psd["eig_min_pre"]),
        "eig_min_post": float(info_psd["eig_min_post"]),
        "tikhonov_applied": bool(info_psd["tikhonov_applied"]),
        "pheno": args.pheno,
        "study": study,
        "bundle": str(args.bundle),
        "bgz": str(bgz),
        "r_aug_mode": "canonical_diag1_normalized_full_p",
        "z_aug_mode": "ctwas_matched_intersection",
        **info_zaug,
    }
    atomic_write_text(
        args.out_dir / f"{args.bundle.stem.replace('_bundle', '')}__{args.pheno}__info_ctwas.json",
        json.dumps(info_dict, indent=2) + "\n"
    )
    sys.stderr.write(f"[prepare-ctwas] DONE wrote {out_path}\n")
    sys.stderr.write(f"[prepare-ctwas] {info_zaug['n_svs_z_aug_dropped']}/{G} SVs got z_SV=0 (insufficient matched support)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
