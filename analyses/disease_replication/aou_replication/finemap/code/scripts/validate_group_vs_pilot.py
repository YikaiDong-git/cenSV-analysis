#!/usr/bin/env python3
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, os.environ["OPENBLAS_NUM_THREADS"])
os.environ.setdefault("FINEMAP_SCRATCH", "/tmp/finemap_scratch")
os.environ.setdefault("PLINK2", "/home/dataproc/bin/plink2")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
import finemap_core as core
import finemap_group as grp

FM = core.FM
TMP = "/tmp/finemap_scratch/_valid_grp"
import csv
row = [r for r in csv.DictReader(open(f"{FM}/targets/target_manifest.tsv"), delimiter="\t")
       if r["fit_id"] == "chr16.MB_280.AMR"][0]
HIT = [s for s in row["sv_ids"].split(",") if s]
OTHERS = [("BI_168.15", ["chr16_sv13_DEL_mon_best8_midpoint1000"]),
          ("DE_660.6", ["chr16_sv1029_DEL_ori_best8_midpoint1000", "chr16_sv465_DEL_non_best8_edge500"])]


def read_vm(d):
    ids, z, p = [], [], []
    with open(f"{d}/variant_meta.tsv") as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t"); ids.append(c[1]); z.append(float(c[4])); p.append(float(c[5]))
    return ids, np.array(z), np.array(p)


def cmp(pilot, grp_dir, label):
    ip, zp, _ = read_vm(pilot); ig, zg, _ = read_vm(grp_dir)
    order_same = ip == ig
    common = [v for v in ip if v in set(ig)]
    zp2 = np.array([zp[ip.index(v)] for v in common]); zg2 = np.array([zg[ig.index(v)] for v in common])
    dz = float(np.max(np.abs(zp2 - zg2))) if common else float("nan")
    Rp = np.load(f"{pilot}/R.npy"); Rg = np.load(f"{grp_dir}/R.npy")
    dR = float(np.max(np.abs(Rp - Rg))) if Rp.shape == Rg.shape else float("nan")
    ok = order_same and len(ip) == len(ig) and dz < 1e-4 and dR < 1e-6
    print(f"  [{label}] order_same={order_same} nP={len(ip)} nG={len(ig)} "
          f"common={len(common)} max|dz|={dz:.2e} max|dR|={dR:.2e} Rshape {Rp.shape}/{Rg.shape} -> {'OK' if ok else 'MISMATCH'}")
    return ok


def main():
    fits = []
    for phe, svs in [("MB_280", HIT)] + OTHERS:
        for mode in ("hit_only", "all_window"):
            fits.append({"phecode": phe, "mcc": "mcc2", "hit_sv_ids": svs, "panel_mode": mode,
                         "snp_select": "assoc_tag", "out_dir": f"{TMP}/{phe}.{mode}",
                         "fit_id": f"chr16.{phe}.AMR"})
    print(f"build_group chr16 AMR: {len(fits)} fits (3 phecodes x 2 panels)")
    res = grp.build_group("chr16", "AMR", fits)
    for r in res:
        print("   ", (r[0].get("fit_id"), r[1], str(r[2])[:50]) if isinstance(r, tuple) else r)
    ok = True
    for mode in ("hit_only", "all_window"):
        pilot = f"{FM}/results/fits/AMR/chr16.MB_280.AMR/{mode}"
        if os.path.exists(f"{pilot}/variant_meta.tsv"):
            ok &= cmp(pilot, f"{TMP}/MB_280.{mode}", f"MB_280.{mode} pilot-vs-group")
        else:
            print(f"  [MB_280.{mode}] no saved pilot to compare")
    print(f"\nRESULT: {'PASS' if ok else 'FAIL — DO NOT run production'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
