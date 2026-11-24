#!/usr/bin/env python3
"""step94b — prove each step94 extract reproduces the FULL-matrix cv_r2 regate, per arm x stratum."""
import csv
import glob
import gc
import gzip
import json
import os
import sys

import numpy as np

AR = "/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication"
RV = f"{AR}/results_vds"
CVALL = f"{RV}/step28_phewide_ext2/sv_cv_r2_all.tsv"
ARMS = ("step30_phewide_ext2", "step30_phewide_mcc1_ext2")
THRS = (0.10, 0.20, 0.30, 0.40)
RELTOL = 1e-3


def cvmap():
    m = {}
    for r in csv.DictReader(open(CVALL), delimiter="\t"):
        v = (r.get("cv_r2") or "").strip()
        if v not in ("", "NA", "nan"):
            m[r["sv_locus"]] = float(v)
    return m


def check(arm, tag, cv):
    per = "/per_ancestry" if tag != "allanc" else ""
    sfx = "" if tag == "allanc" else f"_{tag}"
    md, mf = f"{RV}/{arm}/disease_matrices{per}", f"{RV}/{arm}/disease_matrices_firth{per}"
    man = json.load(open(f"{RV}/{arm}/cloud_topk/cloud_topk_manifest_{tag}.json"))
    K = man["K"]

    sv = [ln.split("\t")[0].strip() for ln in open(f"{md}/disease_svorder{sfx}.txt") if ln.strip()]
    phe = [r["phecode"] for r in csv.DictReader(open(f"{md}/disease_meta{sfx}.tsv"), delimiter="\t")]
    cva = np.array([cv.get(s, np.nan) for s in sv])

    P = np.load(f"{mf}/disease_Pv_final{sfx}.npy", mmap_mode="r").astype(np.float32)
    tested = np.isfinite(P).any(axis=0)
    P = np.where(np.isfinite(P), P, np.float32(1.0))

    exp = {}
    with gzip.open(f"{RV}/{arm}/cloud_topk/cloud_topk{K}_{tag}.tsv.gz", "rt") as fh:
        next(fh)
        for ln in fh:
            a = ln.rstrip("\n").split("\t")
            exp.setdefault(a[0], []).append((int(a[1]), float(a[3])))

    keep = [j for j in range(len(phe)) if tested[j]]
    assert set(exp) == {phe[j] for j in keep}, f"{arm}/{tag}: exported phecode set != tested phecode set"

    bad, dmax = [], 0.0
    for thr in THRS:
        m = cva >= thr
        truth = P[m, :].min(axis=0)
        for j in keep:
            t = float(truth[j])
            cand = [p for i, p in exp[phe[j]] if m[i]]
            if not cand:
                bad.append((thr, phe[j], "no candidate", t))
                continue
            pe = min(cand)
            if abs(pe - t) > RELTOL * max(t, 1e-300):
                bad.append((thr, phe[j], pe, t))
            dmax = max(dmax, abs(np.log10(max(pe, 1e-300)) - np.log10(max(t, 1e-300))))
        del m, truth
        gc.collect()
    del P
    gc.collect()
    flag = "OK " if not bad else "FAIL"
    print(f"  {flag} {arm.replace('step30_phewide_',''):11s} {tag:7s} "
          f"cols={len(keep):5d}  mismatches={len(bad):3d}  max|d(-log10 p)|={dmax:.2e}")
    for b in bad[:3]:
        print(f"        {b}")
    return len(bad)


if __name__ == "__main__":
    cv = cvmap()
    want = os.environ.get("PHEWIDE_STRATA")
    total = 0
    print(f"validating extracts against the full matrices (reltol {RELTOL:g}, thresholds {THRS})")
    for arm in ARMS:
        tags = [os.path.basename(f).split("_manifest_")[1][:-5]
                for f in sorted(glob.glob(f"{RV}/{arm}/cloud_topk/cloud_topk_manifest_*.json"))]
        if want:
            tags = [t for t in tags if t in want.split(",")]
        for t in tags:
            total += check(arm, t, cv)
    print(f"\n{'ALL CLEAN' if total == 0 else f'{total} GENUINE MISMATCHES'}")
    sys.exit(0 if total == 0 else 1)
