#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
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


BASE = Path(CENSV_ROOT + "")
SHARED = Path(CENSV_WORKSPACE + "/_shared")

ONE_KG_VCF_TEMPLATE = str(
    SHARED / "datasets/1KGP/v_high_cov_phased/1000G/"
    "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)

KG_EUR_525_PEDIGREE = (
    BASE / "analyses/sv_catalog/output/ld_panel_pilots/"
    "pedigree/kg_eur_founders_pedigree.txt"
)

ADMIT_TSV_V2 = (
    BASE / "analyses/finemap_multi_sv/output/full_restart/"
    "admit/admit_decompose_cvr2_0p3.tsv"
)

WEIGHTS_BASE = (
    BASE / "analyses/sv_catalog/output/"
    "shared/weights_decompose_two_arm"
)

BCFTOOLS = SHARED / "envs/workspace/v1/.pixi/envs/biotools/bin/bcftools"

MAF_MIN = 0.01


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chrom", required=True, help="e.g. chr20")
    p.add_argument("--region-id", required=True)
    p.add_argument("--region-start", type=int, required=True)
    p.add_argument("--region-end", type=int, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--admit-tsv", type=Path, default=ADMIT_TSV_V2)
    p.add_argument("--weights-base", type=Path, default=WEIGHTS_BASE)
    p.add_argument("--ld-pedigree", type=Path, default=KG_EUR_525_PEDIGREE)
    p.add_argument("--ref-vcf", type=Path, default=None)
    return p.parse_args()


def read_pedigree_samples(path: Path) -> list[str]:
    samples: list[str] = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split()
        sid_col = header.index("SampleID")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > sid_col and parts[sid_col]:
                samples.append(parts[sid_col])
            else:
                parts2 = line.rstrip("\n").split()
                if len(parts2) > sid_col and parts2[sid_col]:
                    samples.append(parts2[sid_col])
    if not samples:
        raise ValueError(f"No samples parsed from {path}")
    return samples


def load_admit_rows_v2(admit_tsv: Path, chrom: str,
                       region_start: int, region_end: int) -> pd.DataFrame:
    df = pd.read_csv(admit_tsv, sep="\t")
    chrom_num = chrom.replace("chr", "")

    base_mask = (df["primary_chrom"].astype(str) == chrom_num)
    if "status" in df.columns:
        base_mask &= (df["status"] == "ADMIT")

    region_mask = (
        (df["primary_block_start"] >= region_start)
        & (df["primary_block_end"] <= region_end)
    )

    out = df[base_mask & region_mask].copy().reset_index(drop=True)
    if len(out) == 0:
        raise ValueError(
            f"No ADMIT rows for chrom={chrom} region [{region_start:,}-{region_end:,}]"
        )
    return out


def parse_locus_id(locus_id: str) -> tuple[str, str]:
    parts = locus_id.split("_")
    if len(parts) >= 6 and parts[4] == "best8":
        return "best8", "_".join(parts[5:])
    if len(parts) >= 7 and parts[4] == "best6" and parts[5] == "alt":
        return "best6alt", "_".join(parts[6:])
    raise ValueError(f"unrecognized locus_id pattern: {locus_id}")


def _shard_dir_for_index(sv_index: int) -> str:
    base = (sv_index // 500) * 500
    return f"sv{base:04d}_{base + 499:04d}"


def weight_tsv_path(weights_base: Path, chrom: str, sv_index: int, locus_id: str) -> Path:
    manifest = weights_base / chrom / f"{chrom}_MANIFEST.tsv"
    mf = pd.read_csv(manifest, sep="\t")
    arm_tag, strategy = parse_locus_id(locus_id)
    mf_row = mf[(mf["sv_index"] == sv_index)
                & (mf["weight_path"].str.contains(f"_{arm_tag}_{strategy}_weights.tsv", regex=False))]
    if len(mf_row) != 1:
        raise ValueError(f"manifest lookup failed for sv_index={sv_index} arm={arm_tag} strategy={strategy} ({len(mf_row)} matches)")
    fname = Path(mf_row.iloc[0]["weight_path"]).name
    shard = _shard_dir_for_index(sv_index)
    p = weights_base / chrom / shard / fname
    if not p.exists():
        candidates = list((weights_base / chrom).rglob(fname))
        if len(candidates) == 1:
            return candidates[0]
        raise FileNotFoundError(
            f"weight TSV not found at {p} or via rglob ({len(candidates)} matches)"
        )
    return p


def load_weight_tsv(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, sep="\t")
    required = {"snp_id", "snp_pos", "snp_chrom", "weight"}
    if not required.issubset(df.columns):
        raise ValueError(f"weight TSV {p} missing required cols (have {df.columns.tolist()})")
    return df


def parse_snp_id(snp_id: str) -> tuple[str, int, str, str]:
    parts = snp_id.split(":")
    if len(parts) != 4:
        raise ValueError(f"snp_id format unexpected: {snp_id}")
    return parts[0], int(parts[1]), parts[2], parts[3]


def assemble_union_snp_meta(weight_dfs: list[pd.DataFrame],
                            region_start: int, region_end: int) -> pd.DataFrame:
    seen: dict[str, tuple[str, int, str, str]] = {}
    for df in weight_dfs:
        for snp_id, pos in zip(df["snp_id"], df["snp_pos"]):
            if region_start <= pos <= region_end and snp_id not in seen:
                ch, p, r, a = parse_snp_id(snp_id)
                seen[snp_id] = (ch, p, r, a)
    rows = [(sid, ch, p, r, a) for sid, (ch, p, r, a) in seen.items()]
    rows.sort(key=lambda x: x[2])
    out = pd.DataFrame(rows, columns=["snp_id", "snp_chrom", "snp_pos", "snp_ref", "snp_alt"])
    return out


def query_panel_genotypes(chrom: str, snp_meta: pd.DataFrame, samples: list[str],
                          ref_vcf: Path, region_start: int, region_end: int):
    sample_str = ",".join(samples)
    region_str = f"{chrom}:{region_start}-{region_end}"

    cmd = [
        str(BCFTOOLS), "query",
        "-r", region_str,
        "-s", sample_str,
        "-f", "%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
        str(ref_vcf),
    ]
    sys.stderr.write(f"[bcftools] {' '.join(cmd[:6])} ...\n")
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    lines = proc.stdout.strip().split("\n")
    if not lines or not lines[0]:
        raise ValueError(f"bcftools returned 0 lines for {region_str}")

    n_samples = len(samples)
    want_ids = set(snp_meta["snp_id"].tolist())
    rows_gt: list[np.ndarray] = []
    rows_meta: list[dict] = []

    for ln in lines:
        cols = ln.split("\t")
        if len(cols) < 5 + n_samples:
            continue
        vcf_id, vcf_chrom, vcf_pos, vcf_ref, vcf_alt = cols[:5]
        chrom_no_prefix = vcf_chrom.replace("chr", "")
        my_id = f"{chrom_no_prefix}:{vcf_pos}:{vcf_ref}:{vcf_alt}"
        if my_id not in want_ids:
            alt_id = f"{chrom_no_prefix}:{vcf_pos}:{vcf_alt}:{vcf_ref}"
            if alt_id not in want_ids:
                continue
            target_id = alt_id
            swap = True
        else:
            target_id = my_id
            swap = False

        gts = cols[5:5 + n_samples]
        dos = np.full(n_samples, np.nan, dtype=np.float64)
        for i, gt in enumerate(gts):
            if gt in ("0|0", "0/0"):
                dos[i] = 0
            elif gt in ("0|1", "1|0", "0/1", "1/0"):
                dos[i] = 1
            elif gt in ("1|1", "1/1"):
                dos[i] = 2
        if swap:
            dos = 2.0 - dos
        m = np.isnan(dos)
        if m.all():
            continue
        if m.any():
            dos[m] = np.nanmean(dos)
        af = dos.sum() / (2.0 * n_samples)
        maf = min(af, 1.0 - af)
        if maf < MAF_MIN:
            continue
        rows_gt.append(dos)
        rows_meta.append({
            "snp_id": target_id, "snp_chrom": chrom,
            "snp_pos": int(vcf_pos), "snp_ref": vcf_ref, "snp_alt": vcf_alt,
            "maf_panel": maf, "swap": int(swap),
        })

    if not rows_gt:
        raise ValueError("No SNPs survived MAF filter / bcftools query")

    G_panel = np.array(rows_gt).T
    kept_meta = pd.DataFrame(rows_meta).reset_index(drop=True)
    return G_panel, kept_meta


def compute_R_snp_snp(G_panel: np.ndarray) -> np.ndarray:
    n, p = G_panel.shape
    mu = G_panel.mean(axis=0, keepdims=True)
    sd = G_panel.std(axis=0, ddof=0, keepdims=True)
    sd[sd == 0] = 1.0
    Gs = (G_panel - mu) / sd
    R = (Gs.T @ Gs) / n
    np.fill_diagonal(R, 1.0)
    return R


def build_B_matrix(weight_dfs: list[pd.DataFrame],
                   kept_meta: pd.DataFrame) -> np.ndarray:
    p = len(kept_meta)
    G = len(weight_dfs)
    B = np.zeros((p, G), dtype=np.float64)
    id_to_idx = {sid: i for i, sid in enumerate(kept_meta["snp_id"].tolist())}

    for j, df in enumerate(weight_dfs):
        for sid, w in zip(df["snp_id"], df["weight"]):
            if sid in id_to_idx:
                idx = id_to_idx[sid]
                B[idx, j] = w
            else:
                ch, pos, ref, alt = parse_snp_id(sid)
                rev = f"{ch}:{pos}:{alt}:{ref}"
                if rev in id_to_idx:
                    idx = id_to_idx[rev]
                    B[idx, j] = -w
    return B


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

    chrom = args.chrom
    region_id = args.region_id
    rs, re_ = args.region_start, args.region_end

    sys.stderr.write(f"[multisv-v2] chrom={chrom} region={region_id} [{rs:,}-{re_:,}]\n")

    admit = load_admit_rows_v2(args.admit_tsv, chrom, rs, re_)
    sys.stderr.write(f"[multisv-v2] ADMIT rows (region-scoped): {len(admit)}\n")

    weight_dfs: list[pd.DataFrame] = []
    sv_meta_rows: list[dict] = []
    for _, row in admit.iterrows():
        sv_idx = int(row["sv_index"])
        locus_id = row["sv_locus"]
        wp = weight_tsv_path(args.weights_base, chrom, sv_idx, locus_id)
        wdf = load_weight_tsv(wp)
        weight_dfs.append(wdf)
        sv_meta_rows.append({
            "sv_locus": locus_id, "sv_index": sv_idx, "sv_id": row["sv_id"],
            "sv_category": row["sv_category"], "strategy": row["strategy"],
            "arm": row["arm"], "cv_r2": row["cv_r2"],
            "n_weight_snps": int(row["n_weight_snps"]),
            "weight_path": str(wp),
        })
    sv_meta = pd.DataFrame(sv_meta_rows)

    snp_meta_req = assemble_union_snp_meta(weight_dfs, rs, re_)
    sys.stderr.write(f"[multisv-v2] requested union SNPs in region: {len(snp_meta_req)}\n")

    samples = read_pedigree_samples(args.ld_pedigree)
    sys.stderr.write(f"[multisv-v2] LD panel samples: {len(samples)}\n")

    ref_vcf = args.ref_vcf or Path(ONE_KG_VCF_TEMPLATE.format(chrom=chrom))
    G_panel, kept_meta = query_panel_genotypes(
        chrom, snp_meta_req, samples, ref_vcf, rs, re_
    )
    sys.stderr.write(f"[multisv-v2] kept SNPs after panel/MAF: {len(kept_meta)} / requested {len(snp_meta_req)}\n")

    R = compute_R_snp_snp(G_panel)
    sys.stderr.write(f"[multisv-v2] R_snp_snp shape: {R.shape}\n")

    B = build_B_matrix(weight_dfs, kept_meta)
    G = B.shape[1]
    n_nonzero_per_sv = (B != 0).sum(axis=0)
    sys.stderr.write(f"[multisv-v2] B shape: {B.shape}; nonzero per SV: min={n_nonzero_per_sv.min()} max={n_nonzero_per_sv.max()} mean={n_nonzero_per_sv.mean():.1f}\n")

    keep_j = n_nonzero_per_sv >= 3
    if not keep_j.all():
        bad_j = np.where(~keep_j)[0]
        sys.stderr.write(f"[multisv-v2] WARN: {len(bad_j)} SVs <3 panel-matched SNPs; dropping\n")
        B = B[:, keep_j]
        sv_meta = sv_meta[keep_j].reset_index(drop=True)
        G = B.shape[1]
        sys.stderr.write(f"[multisv-v2] G after drop: {G}\n")
    n_panel_matched_per_sv_kept = n_nonzero_per_sv[keep_j]

    bundle_path = args.out_dir / f"{region_id}_bundle.npz"
    atomic_write_npz(
        bundle_path,
        R_snp_snp=R,
        B=B,
        snp_id=kept_meta["snp_id"].to_numpy().astype("U64"),
        snp_chrom=kept_meta["snp_chrom"].to_numpy().astype("U16"),
        snp_pos=kept_meta["snp_pos"].to_numpy(dtype=np.int64),
        snp_ref=kept_meta["snp_ref"].to_numpy().astype("U64"),
        snp_alt=kept_meta["snp_alt"].to_numpy().astype("U64"),
        snp_maf_panel=kept_meta["maf_panel"].to_numpy(dtype=np.float64),
        snp_swap=kept_meta["swap"].to_numpy(dtype=np.int32),
        sv_locus=sv_meta["sv_locus"].to_numpy().astype("U128"),
        sv_index=sv_meta["sv_index"].to_numpy(dtype=np.int64),
        sv_category=sv_meta["sv_category"].to_numpy().astype("U32"),
        sv_strategy=sv_meta["strategy"].to_numpy().astype("U32"),
        sv_arm=sv_meta["arm"].to_numpy().astype("U16"),
        sv_cv_r2=sv_meta["cv_r2"].to_numpy(dtype=np.float64),
    )
    sys.stderr.write(f"[multisv-v2] wrote {bundle_path}\n")

    meta_tsv = args.out_dir / f"{region_id}_meta.tsv"
    sv_meta.assign(n_panel_matched_snps=n_panel_matched_per_sv_kept).to_csv(meta_tsv, sep="\t", index=False)

    manifest = {
        "region_id": region_id,
        "chrom": chrom,
        "region_start": rs,
        "region_end": re_,
        "p_snps": int(len(kept_meta)),
        "G_svs": int(G),
        "n_admit_input": int(len(admit)),
        "ld_panel": str(args.ld_pedigree),
        "ld_panel_n": int(len(samples)),
        "ref_vcf": str(ref_vcf),
        "maf_min": float(MAF_MIN),
        "admit_source": str(args.admit_tsv),
        "admit_filter": f"primary_block_start>={rs} AND primary_block_end<={re_}",
        "version": "v2_region_extent_scoped",
        "bundle_sha256": hashlib.sha256(open(bundle_path, "rb").read()).hexdigest(),
    }
    atomic_write_text(args.out_dir / f"{region_id}_manifest.json",
                      json.dumps(manifest, indent=2) + "\n")
    sys.stderr.write(f"[multisv-v2] DONE p={len(kept_meta)} G={G}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
