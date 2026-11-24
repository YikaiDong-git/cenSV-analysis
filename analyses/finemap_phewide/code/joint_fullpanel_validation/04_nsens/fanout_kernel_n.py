#!/usr/bin/env python3
import sys, subprocess
import numpy as np
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
sys.path.insert(0, CENSV_SCRATCH + "/ssv_joint")
import fanout_kernel as fk

load_region = fk.load_region
build_cell = fk.build_cell
TABIX = CENSV_WORKSPACE + "/_shared/envs/workspace/v1/.pixi/envs/biotools/bin/tabix"


def region_n(trait, bgz, chrom_gwas, lo, hi):
    cache = fk._S.setdefault("nc", {})
    if trait not in cache:
        out = subprocess.run([TABIX, bgz, f"{chrom_gwas}:{int(lo)}-{int(hi)}"],
                             capture_output=True, text=True)
        ns = []
        for ln in out.stdout.splitlines():
            f = ln.split("\t")
            if len(f) > 5:
                try:
                    ns.append(float(f[4]))
                except ValueError:
                    pass
        cache[trait] = float(np.median(ns)) if ns else float("nan")
    return cache[trait]
