#!/usr/bin/env python3
"""step94 — compact per-phecode top-K cloud extract, so the PheWAS Manhattan can be REDRAWN off-Workbench."""
import csv
import gc
import gzip
import hashlib
import json
import os

import numpy as np

AR = "/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication"
RV = f"{AR}/results_vds"
CVALL = f"{RV}/step28_phewide_ext2/sv_cv_r2_all.tsv"
K = int(os.environ.get("CLOUD_K", "128"))
ARMS = [os.environ["PHEWIDE_ARM"]] if os.environ.get("PHEWIDE_ARM") else [
    "step30_phewide_ext2", "step30_phewide_mcc1_ext2"]
GRID = (0.10, 0.20, 0.30, 0.40, 0.50, 0.60)
ANCS = ("EUR", "AFR", "AMR", "EAS", "MID", "SAS")
MAXBYTES = 4.9 * 1024 * 1024


def load_cv():
    cv = {}
    for r in csv.DictReader(open(CVALL), delimiter="\t"):
        v = (r.get("cv_r2") or "").strip()
        if v not in ("", "NA", "nan"):
            cv[r["sv_locus"]] = float(v)
    return cv


def run(arm, cvmap, anc=None):
    md = f"{RV}/{arm}/disease_matrices" + ("/per_ancestry" if anc else "")
    mf = f"{RV}/{arm}/disease_matrices_firth" + ("/per_ancestry" if anc else "")
    sfx = f"_{anc}" if anc else ""
    tag = anc or "allanc"
    out = f"{RV}/{arm}/cloud_topk"
    os.makedirs(out, exist_ok=True)

    sv = [ln.split("\t")[0].strip() for ln in open(f"{md}/disease_svorder{sfx}.txt") if ln.strip()]
    phe = [r["phecode"] for r in csv.DictReader(open(f"{md}/disease_meta{sfx}.tsv"), delimiter="\t")]

    P = np.load(f"{mf}/disease_Pv_final{sfx}.npy", mmap_mode="r").astype(np.float32)
    B = np.load(f"{mf}/disease_Be_final{sfx}.npy", mmap_mode="r").astype(np.float32)
    tested = np.isfinite(P).any(axis=0)
    P = np.where(np.isfinite(P), P, np.float32(1.0))
    M, J = P.shape
    assert len(sv) == M, f"{arm}/{tag}: svorder {len(sv)} != matrix rows {M}"
    assert len(phe) == J, f"{arm}/{tag}: meta {len(phe)} != matrix cols {J}"
    assert K <= M, f"{arm}/{tag}: K={K} exceeds {M} SVs"
    keep = np.where(tested)[0]
    print(f"\n=== {arm} [{tag}]  {M} SV x {J} phecode = {M*J:,} cells  (K={K}) ===")
    if len(keep) < J:
        print(f"  {J - len(keep)} of {J} phecodes untested in this stratum -> dropped (not exported as padding)")

    cv = np.array([cvmap.get(s, np.nan) for s in sv], dtype=np.float64)
    missing = int(np.isnan(cv).sum())
    if missing:
        print(f"  !! {missing} SVs have no cv_r2; they are excluded from every gate")
    guarantee, exact_to = {}, None
    for thr in GRID:
        m = cv >= thr
        if not m.any():
            continue
        minp = P[m, :].min(axis=0)
        mx = int(((P < minp[None, :]).sum(axis=0) + 1)[keep].max())
        guarantee[f"{thr:.2f}"] = {"passing_svs": int(m.sum()), "max_rank_needed": mx, "covered_by_K": mx <= K}
        if mx <= K:
            exact_to = thr
        print(f"  cv_r2 >= {thr:.2f}: {int(m.sum()):>6d} SVs pass, max rank {mx:>4d} -> "
              f"{'COVERED' if mx <= K else 'NOT COVERED by K=%d' % K}")
        del m, minp
        gc.collect()

    rows = 0
    parts = ["phecode\tsv_idx\trank\tp\tbeta\n"]
    for j in keep:
        col = P[:, j]
        idx = np.argpartition(col, K - 1)[:K]
        idx = idx[np.argsort(col[idx], kind="stable")]
        pj = phe[j]
        for r, i in enumerate(idx, 1):
            parts.append(f"{pj}\t{i}\t{r}\t{col[i]:.4g}\t{B[i, j]:.4g}\n")
        rows += len(idx)
    blob = "".join(parts).encode()
    del parts, P, B
    gc.collect()

    fp = f"{out}/cloud_topk{K}_{tag}.tsv.gz"
    gz = gzip.compress(blob, 9)
    assert len(gz) < MAXBYTES, f"{fp} is {len(gz)/1048576:.2f} MB, over the 4.9 MB limit - lower CLOUD_K"
    with open(fp, "wb") as fh:
        fh.write(gz)

    man = {
        "arm": arm, "stratum": tag, "K": K, "n_sv": M, "n_phecode": J,
        "n_phecode_tested_here": int(len(keep)), "n_phecode_dropped_untested": int(J - len(keep)),
        "cells_in_full_matrix": M * J,
        "rows_exported": rows, "fraction_of_full_matrix": round(rows / (M * J), 6),
        "svorder_md5": hashlib.md5(open(f"{md}/disease_svorder{sfx}.txt", "rb").read()).hexdigest(),
        "gz_bytes": len(gz), "gz_md5": hashlib.md5(gz).hexdigest(),
        "raw_bytes": len(blob), "sv_without_cv_r2": missing,
        "regate_exact_up_to_cv_r2": exact_to, "guarantee": guarantee,
        "cv_r2_source": "results_vds/step28_phewide_ext2/sv_cv_r2_all.tsv",
        "sv_idx_key": "0-based line number in disease_matrices/disease_svorder.txt",
        "note": ("Per phecode, the K SVs with smallest Firth p. Reproduces the whole-cloud cv_r2 regate "
                 "EXACTLY for any threshold whose max_rank_needed <= K. Do not extrapolate past "
                 "regate_exact_up_to_cv_r2."),
    }
    json.dump(man, open(f"{out}/cloud_topk_manifest_{tag}.json", "w"), indent=2)
    print(f"  wrote {fp}  {len(gz)/1048576:.2f} MB gz ({len(blob)/1048576:.2f} MB raw, {rows:,} rows, "
          f"{100*rows/(M*J):.3f}% of the matrix)")
    print(f"  regate exact up to cv_r2 >= {exact_to}")
    gc.collect()


if __name__ == "__main__":
    cvmap = load_cv()
    print(f"cv_r2 map: {len(cvmap)} SVs from {CVALL}")
    want = os.environ.get("PHEWIDE_STRATA", "allanc")
    if want == "allanc":
        strata = [None]
    elif want == "peranc":
        strata = list(ANCS)
    elif want == "all":
        strata = [None] + list(ANCS)
    else:
        strata = [s.strip() or None for s in want.split(",")]
    for a in ARMS:
        for s in strata:
            run(a, cvmap, s)
