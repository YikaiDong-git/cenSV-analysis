#!/usr/bin/env python3
"""Batched jwbundle builder for the single-SV x all-GWAS SuSiE-RSS sweep (ANTI-THUNDERING)."""
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
    out = subprocess.run([TABIX, bgz, f"{chrom}:{start}-{end}"], capture_output=True, text=True)
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


def match_z(pos, ref, alt, gw):
    p = len(pos); z = np.full(p, np.nan)
    for i in range(p):
        k, ks = (int(pos[i]), ref[i], alt[i]), (int(pos[i]), alt[i], ref[i])
        if k in gw: z[i] = gw[k]
        elif ks in gw: z[i] = -gw[ks]
    return z


def impg_fill(z, R, keep, n_ld=525):
    obs = np.where(keep)[0]; mis = np.where(~keep)[0]
    zf = z.copy()
    if mis.size and obs.size:
        lam = 2.0 / np.sqrt(n_ld)
        C = R[np.ix_(obs, obs)]
        A = (1 - lam) * C + lam * np.eye(len(obs))
        try:
            sol = np.linalg.solve(A, z[obs])
            for u in mis:
                c = R[u, obs]
                zf[u] = (1 - lam) * float(c @ sol)
        except np.linalg.LinAlgError:
            zf[mis] = 0.0
    zf[np.isnan(zf)] = 0.0
    return zf


def write_jw(path, R_aug, z_aug, b, sv_keep):
    np.savez(path, R_aug=R_aug, z_aug=z_aug,
             sv_locus=b["sv_locus"][sv_keep], sv_index=b["sv_index"][sv_keep],
             sv_category=b["sv_category"][sv_keep], sv_strategy=b["sv_strategy"][sv_keep],
             sv_arm=b["sv_arm"][sv_keep], sv_cv_r2=b["sv_cv_r2"][sv_keep])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True); ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True); ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--arm", choices=["matched", "impute"], required=True)
    ap.add_argument("--traits", required=True)
    ap.add_argument("--smoke-cols", default="")
    ap.add_argument("--modes", default="allwindow,hitbundle")
    ap.add_argument("--work-dir", required=True); ap.add_argument("--manifest-out", required=True)
    ap.add_argument("--n-ld", type=int, default=525)
    a = ap.parse_args()
    os.makedirs(a.work_dir, exist_ok=True)
    region = os.path.basename(a.bundle).replace("_bundle.npz", "")
    smoke_cols = [int(x) for x in a.smoke_cols.split(",") if x.strip() != ""]
    modes = set(m for m in a.modes.split(",") if m)

    b = np.load(a.bundle, allow_pickle=True)
    R = b["R_snp_snp"].astype(np.float64); B = b["B"].astype(np.float64)
    pos = b["snp_pos"].astype(int)
    ref = np.array([str(x).upper() for x in b["snp_ref"]]); alt = np.array([str(x).upper() for x in b["snp_alt"]])
    p, G = B.shape

    traits = [ln.split("\t") for ln in open(a.traits).read().splitlines() if ln.strip()]
    man = open(a.manifest_out, "a")
    for name, bgz in traits:
        gw = load_gwas(bgz, a.chrom, a.start, a.end)
        z = match_z(pos, ref, alt, gw)
        keep = np.isfinite(z)
        if keep.sum() < 5:
            continue
        if a.arm == "matched":
            obs = np.where(keep)[0]
            Roo = R[np.ix_(obs, obs)]; Boo = B[obs, :]; zoo = z[obs]
        else:
            zf = impg_fill(z, R, keep, a.n_ld)
            obs = np.arange(p); Roo = R; Boo = B; zoo = zf
        nz = (Boo != 0).sum(axis=0)
        quad = np.einsum("ij,ij->j", Boo, Roo @ Boo)
        sv_keep = np.where((nz >= 3) & (quad > 0))[0]
        if sv_keep.size == 0:
            continue
        if "allwindow" in modes:
            Bk = Boo[:, sv_keep]
            R_aug, D = build_canonical_raug(Roo, Bk); R_aug, _ = psd_safeguard(R_aug)
            z_aug = build_z_aug(zoo, Bk, D)
            aw = f"{a.work_dir}/{region}__{name}__{a.arm}__allwindow.npz"
            write_jw(aw, R_aug, z_aug, b, sv_keep)
            man.write(f"{aw}\tallwindow\t{a.arm}\t{name}\t{region}\t*\t*\n")
        if "hitbundle" in modes:
            for j in smoke_cols:
                if j not in sv_keep:
                    continue
                Bk1 = Boo[:, [j]]
                R_aug1, D1 = build_canonical_raug(Roo, Bk1); R_aug1, _ = psd_safeguard(R_aug1)
                z_aug1 = build_z_aug(zoo, Bk1, D1)
                hb = f"{a.work_dir}/{region}__{name}__{a.arm}__hitbundle__sv{j}.npz"
                write_jw(hb, R_aug1, z_aug1, b, np.array([j]))
                man.write(f"{hb}\thitbundle\t{a.arm}\t{name}\t{region}\t{str(b['sv_locus'][j])}\t{j}\n")
    man.close()
    print(f"[build_cells] {region} arm={a.arm} modes={sorted(modes)}: {len(traits)} traits -> {a.work_dir}")


if __name__ == "__main__":
    main()
