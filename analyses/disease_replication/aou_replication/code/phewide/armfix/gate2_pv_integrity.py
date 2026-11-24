#!/usr/bin/env python3
import os
import numpy as np

HOME = os.path.expanduser("~")
RV = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds"
NEWD = f"{RV}/step30_phewide/disease_matrices"
OLDD = f"{RV}/archive_2026-07-11_pre_armfix/step30_phewide"
AFFECTED = {"chr1", "chr9", "chr16"}
TOL = 0.1


def load(d):
    Pv = np.load(f"{d}/disease_Pv.npy")
    svs = [ln.split("\t")[0] for ln in open(f"{d}/disease_svorder.txt") if ln.strip()]
    import csv
    phs = [r["phecode"] for r in csv.DictReader(open(f"{d}/disease_meta.tsv"), delimiter="\t")]
    return Pv, svs, phs


def main():
    nPv, nsv, nph = load(NEWD)
    oPv, osv, oph = load(OLDD)
    print(f"new Pv {nPv.shape} ({len(nsv)} SV x {len(nph)} phe); old Pv {oPv.shape}", flush=True)

    ocol = {p: j for j, p in enumerate(oph)}
    col_new = [j for j, p in enumerate(nph) if p in ocol]
    col_old = [ocol[nph[j]] for j in col_new]
    orow = {s: i for i, s in enumerate(osv)}
    nl10 = lambda P: -np.log10(np.clip(P, 1e-300, 1.0))

    unaff_max = 0.0; unaff_n = 0; aff_changed = 0; aff_n = 0; bad = []
    for i, sv in enumerate(nsv):
        if sv not in orow:
            continue
        a = nl10(nPv[i, col_new]); b = nl10(oPv[orow[sv], col_old])
        fin = np.isfinite(a) & np.isfinite(b)
        if fin.sum() == 0:
            continue
        mx = float(np.max(np.abs(a[fin] - b[fin])))
        chrom = sv.split("_")[0]
        if chrom in AFFECTED:
            aff_n += 1
            if mx > TOL:
                aff_changed += 1
        else:
            unaff_n += 1
            unaff_max = max(unaff_max, mx)
            if mx > TOL:
                bad.append((sv, mx))

    print("\n=== GATE 2 RESULT ===")
    print(f"  UNAFFECTED (18 chroms): {unaff_n} SV rows compared; max |Δ -log10 p| = {unaff_max:.2e} (tol {TOL:.0e})")
    print(f"  AFFECTED  (chr1/9/16):  {aff_changed}/{aff_n} SV rows changed (expect ~all)")
    if bad:
        bad.sort(key=lambda x: -x[1])
        print(f"  !! {len(bad)} unaffected SV rows exceed tol, worst: {bad[:6]}")
    ok = (len(bad) == 0) and (aff_changed > 0)
    print(f"\n  GATE 2: {'PASS ✓' if ok else 'FAIL ✗'}")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
