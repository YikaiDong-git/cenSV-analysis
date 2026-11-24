#!/usr/bin/env python3
"""Step 71 (ext2 arm) — run the whole ext2 disease chain against the MCC=1 phenotypes."""
import os, sys, csv, time, hashlib, argparse

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import step51_disease_scoretest_ext2 as S      # noqa: E402  (loads numpy — pinning above must precede it)

RV = S.RV
MCC1_INP = f"{RV}/inputs_mcc1_ext2"
MCC1_OUT = f"{RV}/step30_phewide_mcc1_ext2"
PRIMARY_OUT = S.OUTA

S.INP = MCC1_INP
S.OUTA = MCC1_OUT
S.MATDIR = f"{MCC1_OUT}/disease_matrices"
S.PSUM = f"{MCC1_OUT}/phecodeX_pheno_summary.tsv"
os.makedirs(S.MATDIR, exist_ok=True)
assert "mcc1_ext2" in S.INP and "mcc1_ext2" in S.OUTA and "mcc1_ext2" in S.MATDIR, "MCC1 patch failed"

import step51b_exact_firth_refit_ext2 as H       # noqa: E402  (binds S.INP/S.OUTA at import — after patch)
import step52_approx_firth_recalibrate_ext2 as H2  # noqa: E402
import step53_build_final_firth_matrix_ext2 as H3  # noqa: E402

for mod, name in ((H, "step51b"), (H2, "step52"), (H3, "step53")):
    for attr in ("INP", "OUTA", "MATDIR", "OUTDIR"):
        v = getattr(mod, attr, None)
        if isinstance(v, str) and "step30_phewide_ext2" in v:
            raise SystemExit(f"ABORT: {name}.{attr} = {v} still points at the PRIMARY ext2 tree — "
                             f"the import happened before the patch and this arm would overwrite it.")


SIG_SKIP_DIRS = {"gates"}


def tree_sig(root):
    h = hashlib.md5(); n = 0
    for dp, _d, fs in os.walk(root):
        _d[:] = [d for d in _d if d not in SIG_SKIP_DIRS or os.path.relpath(dp, root) != "."]
        rel_dir = os.path.relpath(dp, root)
        if rel_dir.split(os.sep)[0] in SIG_SKIP_DIRS:
            continue
        for fn in sorted(fs):
            p = os.path.join(dp, fn)
            if os.path.islink(p) or not os.path.isfile(p):
                continue
            hh = hashlib.md5()
            with open(p, "rb") as f:
                for c in iter(lambda: f.read(1 << 20), b""):
                    hh.update(c)
            h.update(os.path.relpath(p, root).encode()); h.update(hh.hexdigest().encode()); n += 1
    return h.hexdigest(), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--skip-score", action="store_true", help="resume: skip step51 if matrices exist")
    ap.add_argument("--verify-isolation", action="store_true", default=True)
    a = ap.parse_args()
    t0 = time.time()

    for f in ("disease.pheno", "covar.txt"):
        p = f"{MCC1_INP}/{f}"
        if not os.path.exists(p):
            raise SystemExit(f"ABORT: {p} missing — run step69_build_disease_pheno_mcc1_ext2.py first")
    nph = len(open(f"{MCC1_INP}/disease.pheno").readline().rstrip("\n").split("\t")) - 2
    print(f"MCC=1 inputs: {nph} phecodes in disease.pheno", flush=True)
    print(f"  INP    = {S.INP}\n  OUTA   = {S.OUTA}\n  MATDIR = {S.MATDIR}", flush=True)

    pre = tree_sig(PRIMARY_OUT) if a.verify_isolation else None
    if pre:
        print(f"primary ext2 tree fingerprint (pre): {pre[0]} over {pre[1]} files", flush=True)

    if a.skip_score and os.path.exists(f"{S.MATDIR}/disease_Pv.npy"):
        print("\n[1/4] step51 score test — matrices exist, skipping", flush=True)
    else:
        print("\n[1/4] step51 score test (MCC=1)", flush=True)
        sys.argv = ["step51", "--batch", str(a.batch)]
        S.main()
    print(f"\n[2/4] step51b exact Firth ({time.time()-t0:.0f}s)", flush=True)
    sys.argv = ["step51b"]; H.main()
    print(f"\n[3/4] step52 approx Firth ({time.time()-t0:.0f}s)", flush=True)
    sys.argv = ["step52"]; H2.main()
    print(f"\n[4/4] step53 final Firth matrix ({time.time()-t0:.0f}s)", flush=True)
    sys.argv = ["step53"]; H3.main()

    if pre:
        post = tree_sig(PRIMARY_OUT)
        ok = post[0] == pre[0]
        print(f"\nISOLATION: primary ext2 tree {'UNCHANGED' if ok else '*** CHANGED ***'} "
              f"({post[1]} files, {post[0]})", flush=True)
        if not ok:
            raise SystemExit("ABORT: the MCC=1 arm modified the primary ext2 tree.")
    print(f"\nMCC1_EXT2_CHAIN_DONE ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
