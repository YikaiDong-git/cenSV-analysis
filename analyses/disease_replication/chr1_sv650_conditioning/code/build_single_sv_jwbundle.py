#!/usr/bin/env python3
import sys, os, argparse, subprocess
import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

LIB = CENSV_ROOT + "/analyses/finemap_multi_sv/code/lib"
sys.path.insert(0, LIB)
from raug import build_canonical_raug, build_z_aug, psd_safeguard  # noqa
TABIX = CENSV_WORKSPACE + "/_shared/envs/workspace/v1/.pixi/envs/biotools/bin/tabix"


def load_gwas(bgz, chrom, start, end):
    out = subprocess.run([TABIX, bgz, f"{chrom}:{start}-{end}"], capture_output=True, text=True, check=True)
    d = {}
    for ln in out.stdout.splitlines():
        f = ln.split("\t")
        try:
            parts = f[2].split("_")
            ref, alt = parts[-2].upper(), parts[-1].upper()
            beta, se = float(f[5]), float(f[6])
            if se > 0:
                d[(int(f[1]), ref, alt)] = beta / se
        except Exception:
            pass
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True); ap.add_argument("--gwas-bgz", required=True)
    ap.add_argument("--chrom", required=True); ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True); ap.add_argument("--pheno", required=True)
    ap.add_argument("--region-id", required=True); ap.add_argument("--out-dir", required=True)
    ap.add_argument("--only-sv-locus", required=True, help="restrict SuSiE SV block to this one sv_locus")
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True)

    b = np.load(a.bundle, allow_pickle=True)
    R = b["R_snp_snp"].astype(np.float64); B = b["B"].astype(np.float64)
    pos = b["snp_pos"].astype(int)
    ref = np.array([str(x).upper() for x in b["snp_ref"]])
    alt = np.array([str(x).upper() for x in b["snp_alt"]])
    p, G = B.shape

    gw = load_gwas(a.gwas_bgz, a.chrom, a.start, a.end)
    z = np.full(p, np.nan)
    for i in range(p):
        k, ks = (int(pos[i]), ref[i], alt[i]), (int(pos[i]), alt[i], ref[i])
        if k in gw: z[i] = gw[k]
        elif ks in gw: z[i] = -gw[ks]
    obs = np.where(np.isfinite(z))[0]
    if len(obs) < 5:
        print(f"[single] {a.region_id} {a.pheno}: only {len(obs)} observed SNPs — SKIP"); return

    Roo = R[np.ix_(obs, obs)]; Boo = B[obs, :]; zoo = z[obs]
    nz = (Boo != 0).sum(axis=0)
    quad = np.array([float(Boo[:, j] @ (Roo @ Boo[:, j])) for j in range(G)])
    sv_keep = np.where((nz >= 3) & (quad > 0))[0]

    svloc = np.array([str(x) for x in b["sv_locus"]])
    want_idx = np.where(svloc == a.only_sv_locus)[0]
    if len(want_idx) != 1:
        raise SystemExit(f"[single] sv_locus {a.only_sv_locus!r} matched {len(want_idx)} columns "
                         f"(expected exactly 1) — non-unique key, ABORT")
    want = set(want_idx.tolist())
    sv_keep = np.array([j for j in sv_keep if j in want], dtype=int)
    if len(sv_keep) == 0:
        print(f"[single] {a.only_sv_locus} NOT kept for {a.pheno} "
              f"(fails >=3 observed support or quad<=0) — SKIP"); return
    assert len(sv_keep) == 1, f"expected exactly 1 SV after restriction, got {len(sv_keep)}"

    Bk = Boo[:, sv_keep]
    R_aug, D = build_canonical_raug(Roo, Bk)
    R_aug, _ = psd_safeguard(R_aug)
    z_aug = build_z_aug(zoo, Bk, D)

    out = dict(R_aug=R_aug, z_aug=z_aug,
               sv_locus=b["sv_locus"][sv_keep], sv_index=b["sv_index"][sv_keep],
               sv_category=b["sv_category"][sv_keep], sv_strategy=b["sv_strategy"][sv_keep],
               sv_arm=b["sv_arm"][sv_keep], sv_cv_r2=b["sv_cv_r2"][sv_keep],
               snp_id=b["snp_id"][obs], snp_chrom=b["snp_chrom"][obs], snp_pos=b["snp_pos"][obs])
    fn = f"{a.out_dir}/{a.region_id}__{a.pheno}__jw.npz"
    np.savez(fn, **out)
    print(f"[single] {a.region_id} {a.pheno}: p_obs={len(obs)}/{p} G_kept=1 (sv_locus={a.only_sv_locus}, "
          f"orig_G={G}, nz={int(nz[sv_keep[0]])}, quad={quad[sv_keep[0]]:.4f}) -> {fn}")


if __name__ == "__main__":
    main()
