#!/usr/bin/env python3
import sys, os, re
import numpy as np
CODE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CODE)
import build_cells as bc                                    # noqa  (reuse validated kernel fns)
from raug import build_canonical_raug, build_z_aug, psd_safeguard  # noqa

cand, work, manifest = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(work, exist_ok=True)
rows = [l.split("\t") for l in open(cand).read().splitlines()]
idx = {c: i for i, c in enumerate(rows[0])}
cache = {}
man = open(manifest, "w")
n = 0
for r in rows[1:]:
    rb, chrom, start, end = r[idx["region_bundle"]], r[idx["chrom"]], int(r[idx["start"]]), int(r[idx["end"]])
    col, arm, trait, bgz = int(r[idx["col"]]), r[idx["arm"]], r[idx["trait"]], r[idx["bgz"]]
    phys, region, reason = r[idx["phys"]], r[idx["region"]], r[idx["reason"]]
    if not bgz:
        continue
    if rb not in cache:
        b = np.load(rb, allow_pickle=True)
        cache[rb] = (b["R_snp_snp"].astype(np.float64), b["B"].astype(np.float64), b["snp_pos"].astype(int),
                     np.array([str(x).upper() for x in b["snp_ref"]]),
                     np.array([str(x).upper() for x in b["snp_alt"]]))
    R, B, pos, ref, alt = cache[rb]
    gw = bc.load_gwas(bgz, chrom, start, end)
    z = bc.match_z(pos, ref, alt, gw)
    keep = np.isfinite(z)
    if keep.sum() < 5:
        continue
    if arm == "matched":
        obs = np.where(keep)[0]; Roo = R[np.ix_(obs, obs)]; Boo = B[obs, :]; zoo = z[obs]
    else:
        zoo = bc.impg_fill(z, R, keep); Roo = R; Boo = B
    Bk1 = Boo[:, [col]]
    R_aug, D = build_canonical_raug(Roo, Bk1); R_aug, _ = psd_safeguard(R_aug)
    z_aug = build_z_aug(zoo, Bk1, D)
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", f"{phys}__{trait}__{arm}")
    npz = f"{work}/{safe}.npz"
    np.savez(npz, R_aug=R_aug, z_aug=z_aug)
    man.write(f"{npz}\t{phys}\t{trait}\t{arm}\t{region}\t{reason}\t{R_aug.shape[0]-1}\n")
    n += 1
man.close()
print(f"[srss_build] {n} candidate cells -> {work} (manifest {manifest})")
