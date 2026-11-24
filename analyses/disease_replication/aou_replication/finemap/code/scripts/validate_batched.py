#!/usr/bin/env python3
"""validate_batched — gate the batched-score engine (finemap_group.build_group) before any production run."""
import os
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, os.environ["OPENBLAS_NUM_THREADS"])
os.environ.setdefault("FINEMAP_SCRATCH", "/tmp/finemap_scratch")
os.environ.setdefault("PLINK2", "/home/dataproc/bin/plink2")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import numpy as np  # noqa: E402
import finemap_core as core  # noqa: E402
import finemap_group as grp  # noqa: E402

FM = core.FM
PILOT = f"{FM}/results/fits/AMR/chr16.MB_280.AMR/hit_only"
TMP = "/tmp/finemap_scratch/_validate"

MB280_SVS = ("chr16_sv1029_DEL_ori_best8_midpoint1000,chr16_sv1029_DEL_ori_best6_alt_edge2000,"
             "chr16_sv47_INS_ori_best6_alt_edge2000,chr16_sv14_DEL_ori_best8_fksize100k,"
             "chr16_sv55_DEL_non_best8_midpoint1000,chr16_sv55_DEL_non_best6_alt_edge2000,"
             "chr16_sv37_INS_ori_best6_alt_edge2000,chr16_sv13_DEL_mon_best8_midpoint1000,"
             "chr16_sv46_INS_mon_best6_alt_edge2000,chr16_sv13_DEL_mon_best6_alt_edge2000,"
             "chr16_sv90_DEL_non_best8_midpoint1000,chr16_sv1027_DEL_non_best8_edge2000,"
             "chr16_sv73_INS_non_best8_edge2000,chr16_sv15_DEL_res_best8_fksize100k,"
             "chr16_sv22_DEL_res_best8_edge500,chr16_sv98_DEL_non_best8_fksize200k,"
             "chr16_sv495_INS_non_best8_fksize100k,chr16_sv88_DEL_non_best8_fksize200k,"
             "chr16_sv71_DEL_non_best8_midpoint1000,chr16_sv82_DEL_non_best8_edge500,"
             "chr16_sv66_DEL_non_best8_edge500,chr16_sv412_DEL_non_best8_edge2000,"
             "chr16_sv79_INS_non_best6_alt_edge2000,chr16_sv329_INS_non_best8_edge500,"
             "chr16_sv230_INS_non_best8_midpoint1000,chr16_sv465_DEL_non_best8_edge500,"
             "chr16_sv103_DEL_non_best6_alt_fksize100k,chr16_sv1017_INS_non_best8_fksize200k,"
             "chr16_sv559_DEL_non_best6_alt_edge1000,chr16_sv102_DEL_non_best8_fksize100k,"
             "chr16_sv230_INS_non_best6_alt_edge2000").split(",")


def read_varmeta(d):
    ids, z, p = [], [], []
    with open(f"{d}/variant_meta.tsv") as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            ids.append(c[1]); z.append(float(c[4])); p.append(float(c[5]))
    return ids, np.array(z), np.array(p)


def compare(a, b, label):
    ia, za, pa = read_varmeta(a)
    ib, zb, pb = read_varmeta(b)
    ok = True
    if ia != ib:
        sa, sb = set(ia), set(ib)
        print(f"  [{label}] VARIANT SET DIFF: only_a={len(sa-sb)} only_b={len(sb-sa)} lenA={len(ia)} lenB={len(ib)}")
        common = [v for v in ia if v in sb]
        if not common:
            return False
        za = np.array([za[ia.index(v)] for v in common]); zb = np.array([zb[ib.index(v)] for v in common])
        ok = False
    else:
        print(f"  [{label}] variant order IDENTICAL ({len(ia)} vars)")
    dz = float(np.max(np.abs(za - zb))) if len(za) == len(zb) else float("nan")
    Ra = np.load(f"{a}/R.npy"); Rb = np.load(f"{b}/R.npy")
    dR = float(np.max(np.abs(Ra - Rb))) if Ra.shape == Rb.shape else float("nan")
    print(f"  [{label}] max|dz|={dz:.3e}  max|dR|={dR:.3e}  R.shape a={Ra.shape} b={Rb.shape}")
    return ok and dz < 1e-4 and dR < 1e-6


def main():
    hit = MB280_SVS
    bf = f"{TMP}/buildfit/MB_280.hit_only"
    print("=== build_fit MB_280 hit_only (refactored) ===")
    core.build_fit("chr16", "MB_280", "AMR", "mcc2", hit, bf, panel_mode="hit_only", force=False)

    fits = []
    for phe, svs in [("MB_280", hit),
                     ("BI_168.15", ["chr16_sv13_DEL_mon_best8_midpoint1000"]),
                     ("DE_660.6", ["chr16_sv1029_DEL_ori_best8_midpoint1000", "chr16_sv465_DEL_non_best8_edge500"])]:
        for mode in ("hit_only", "all_window"):
            fits.append({"phecode": phe, "mcc": "mcc2", "hit_sv_ids": svs, "panel_mode": mode,
                         "snp_select": "assoc_tag", "out_dir": f"{TMP}/buildgroup/{phe}.{mode}",
                         "fit_id": f"chr16.{phe}.AMR"})
    print("=== build_group chr16 AMR (3 phecodes x 2 panels) ===")
    res = grp.build_group("chr16", "AMR", fits)
    for r in res:
        print("   ", (r[0].get("fit_id"), r[1], str(r[2])[:60]) if isinstance(r, tuple) else r)

    print("\n=== CHECK 1 refactor: build_fit vs SAVED PILOT (MB_280 hit_only) ===")
    c1 = compare(PILOT, bf, "pilot-vs-buildfit")
    print("=== CHECK 2 engine: build_fit vs build_group (MB_280 hit_only) ===")
    c2 = compare(bf, f"{TMP}/buildgroup/MB_280.hit_only", "buildfit-vs-buildgroup")

    print(f"\nRESULT: refactor_ok={c1}  engine_ok={c2}  => "
          f"{'PASS' if (c1 and c2) else 'FAIL — DO NOT run production'}")
    sys.exit(0 if (c1 and c2) else 1)


if __name__ == "__main__":
    main()
