#!/usr/bin/env python3
import glob, numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
BUNDLES = CENSV_ROOT + "/analyses/finemap_multi_sv/output/bundles"

def wall(p):
    return max(0.15, 0.30 + 0.00050 * p)

rows = []
tot_g = 0
phys = set()
cost = 0.0
for npz in sorted(glob.glob(f"{BUNDLES}/chr*/*_bundle.npz")):
    b = np.load(npz, allow_pickle=True)
    p = int(len(b["snp_pos"])); G = int(len(b["sv_locus"]))
    sl = [str(x) for x in b["sv_locus"]]
    for name in sl:
        phys.add("_".join(name.split("_")[:2]))
    tot_g += G
    cost += G * wall(p)
    rows.append((npz.split("/")[-1], p, G))

for name, p, G in rows:
    print(f"{name:48s} p={p:5d} G={G:5d}")
n_traits, n_arms = 464, 2
total_fits = tot_g * n_traits * n_arms
core_sec = cost * n_traits * n_arms
print(f"\nTOTAL: regions={len(rows)}  sum_G(columns)={tot_g}  distinct_physical_SVs={len(phys)}")
print(f"hit_bundle fits = sum_G x {n_traits} x {n_arms} = {total_fits:,}")
print(f"est core-seconds (hit_bundle) = {core_sec:,.0f}  = {core_sec/3600:,.0f} core-hours")
print(f"  on 120 cores: {core_sec/3600/120:.1f} h  |  on 200 cores: {core_sec/3600/200:.1f} h")
