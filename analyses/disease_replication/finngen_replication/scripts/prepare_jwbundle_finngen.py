#!/usr/bin/env python3
import sys
import os
import json
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

MS = Path(CENSV_ROOT + "/analyses/finemap_multi_sv")
sys.path.insert(0, str(MS / "code" / "02_prepare_jwbundle_multisv"))
from prepare_jwbundle_multisv import (  # noqa: E402
    harmonize_gwas, load_bundle, atomic_write_npz, atomic_write_text)
sys.path.insert(0, str(MS / "code" / "lib"))
from raug import build_canonical_raug, build_z_aug, psd_safeguard  # noqa: E402


def read_region_tsv(tsv, chrom, rs, re_):
    rows = []
    with open(tsv) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            p = f[0].split("_")
            if len(p) < 4 or p[0] != chrom:
                continue
            pos = int(p[1])
            if pos < rs or pos > re_:
                continue
            rows.append((f[0], float(f[3]), float(f[4])))
    return pd.DataFrame(rows, columns=["SNP", "BETA", "SE"])


def main():
    bundle_p, htsv, pheno, out_dir = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4])
    out_dir.mkdir(parents=True, exist_ok=True)
    b = load_bundle(Path(bundle_p))
    R = b["R_snp_snp"]
    B = b["B"]
    snp_meta = {k: b[k] for k in ("snp_id", "snp_chrom", "snp_pos", "snp_ref", "snp_alt", "snp_maf_panel", "snp_swap")}
    p, G = B.shape
    chrom = str(snp_meta["snp_chrom"][0])
    rs = int(snp_meta["snp_pos"].min())
    re_ = int(snp_meta["snp_pos"].max())

    gwas_df = read_region_tsv(htsv, chrom, rs, re_)
    sys.stderr.write(f"[prep-fg] {pheno}: region {chrom}:{rs}-{re_} gwas_rows={len(gwas_df)}\n")

    z_snp, keep = harmonize_gwas(snp_meta, gwas_df)
    n_matched = int(keep.sum())
    sys.stderr.write(f"[prep-fg] matched {n_matched}/{p}\n")
    if n_matched < 3:
        raise ValueError(f"Only {n_matched} SNPs harmonized")

    R_aug, D = build_canonical_raug(R, B)
    z_aug = build_z_aug(z_snp, B, D)
    R_aug_safe, info = psd_safeguard(R_aug)
    sys.stderr.write(f"[prep-fg] psd eig_min_pre={info['eig_min_pre']:.3e} post={info['eig_min_post']:.3e} "
                     f"tikhonov={info['tikhonov_applied']} n_clipped={info['n_clipped']}\n")

    stem = Path(bundle_p).stem.replace("_bundle", "")
    out_path = out_dir / f"{stem}__{pheno}__jw.npz"
    atomic_write_npz(
        out_path, R_aug=R_aug_safe, z_aug=z_aug, D_normalizer=D, z_snp=z_snp,
        snp_matched=keep.astype(np.int32),
        sv_locus=b["sv_locus"], sv_index=b["sv_index"], sv_category=b["sv_category"],
        sv_strategy=b["sv_strategy"], sv_arm=b["sv_arm"], sv_cv_r2=b["sv_cv_r2"],
        snp_id=snp_meta["snp_id"], snp_chrom=snp_meta["snp_chrom"], snp_pos=snp_meta["snp_pos"])
    info_dict = dict(
        p_snps=int(p), G_svs=int(G), n_matched=n_matched, n_gwas_region=int(len(gwas_df)),
        eig_min_pre=float(info["eig_min_pre"]), eig_min_post=float(info["eig_min_post"]),
        tikhonov_applied=bool(info["tikhonov_applied"]), n_clipped=int(info["n_clipped"]),
        pheno=pheno, study="finngen_R12", bundle=str(bundle_p),
        canonical_R_aug=True, r_aug_mode="canonical_diag1_normalized",
        L=10, max_iter=200, check_prior=False, susieR_version="0.14.2")
    atomic_write_text(out_dir / f"{stem}__{pheno}__info.json", json.dumps(info_dict, indent=2) + "\n")
    sys.stderr.write(f"[prep-fg] DONE wrote {out_path}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
