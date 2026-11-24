#!/usr/bin/env python3
import sys, os, re
import numpy as np
CODE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, CODE)
import build_cells as bc                                    # noqa (reuse validated load_gwas/match_z)

bundle, traits_tsv, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out_dir, exist_ok=True)
region = os.path.basename(bundle).replace("_bundle.npz", "")
m = re.match(r"chr(\d+)_super\d+_(\d+)_(\d+)", region)
chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))

b = np.load(bundle, allow_pickle=True)
pos = b["snp_pos"].astype(int)
ref = np.array([str(x).upper() for x in b["snp_ref"]])
alt = np.array([str(x).upper() for x in b["snp_alt"]])
R = b["R_snp_snp"].astype(np.float64)
p = len(pos)
traits = [ln.split("\t") for ln in open(traits_tsv).read().splitlines() if ln.strip()]
T = len(traits)

Z = np.full((p, T), np.nan, dtype=np.float32)
Zf = np.zeros((p, T), dtype=np.float32)
names = []
for t, (name, bgz) in enumerate(traits):
    gw = bc.load_gwas(bgz, chrom, start, end)
    z = bc.match_z(pos, ref, alt, gw)
    Z[:, t] = z
    keep = np.isfinite(z)
    if keep.sum() >= 5:
        Zf[:, t] = bc.impg_fill(z, R, keep)
    names.append(name)

nz = np.isfinite(Z).sum(0).astype(np.int32)
out = f"{out_dir}/{region}__zmat.npz"
np.savez(out, z=Z, zf=Zf, nobs=nz, traits=np.array(names), pos=pos,
         region=region, chrom=chrom, start=start, end=end)
print(f"[precompute_z] {region} p={p} T={T} median_obs/trait={int(np.median(nz))} "
      f"min_obs={int(nz.min())} -> {out}")
