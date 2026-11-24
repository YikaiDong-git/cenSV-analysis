#!/usr/bin/env python3
"""Step 83 (ext2 arm) — finish the MCC=1 chain: stages 3-4 only (step52 approx Firth -> step53 merge+BH)."""
import os, sys, time, argparse, hashlib

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S      # noqa: E402  (loads numpy — pinning above must precede it)

RV = S.RV
MCC1_INP = f"{RV}/inputs_mcc1_ext2"
MCC1_OUT = f"{RV}/step30_phewide_mcc1_ext2"
PRIMARY_OUT = S.OUTA

S.INP = MCC1_INP
S.OUTA = MCC1_OUT
S.MATDIR = f"{MCC1_OUT}/disease_matrices"
S.PSUM = f"{MCC1_OUT}/phecodeX_pheno_summary.tsv"
assert S.INP == MCC1_INP and S.OUTA == MCC1_OUT, "engine patch did not take"

import step52_approx_firth_recalibrate_ext2 as H2   # noqa: E402
import step53_build_final_firth_matrix_ext2 as H3   # noqa: E402

_BINDINGS = [("step51.INP", S.INP), ("step51.OUTA", S.OUTA), ("step51.MATDIR", S.MATDIR),
             ("step51.PSUM", S.PSUM), ("step52.OUTDIR", H2.OUTDIR), ("step53.OUT", H3.OUT)]
_bad = [(n, v) for n, v in _BINDINGS if not v.startswith((MCC1_OUT, MCC1_INP))]
if _bad:
    for n, v in _bad:
        print(f"  LEAK {n} -> {v}", flush=True)
    raise SystemExit("ABORT: a downstream module bound a path OUTSIDE the MCC=1 tree (import order).")

SIG_SKIP_DIRS = {"gates"}


def content_sig(root):
    ent = {}
    for dp, _d, fn in os.walk(root):
        parts = os.path.relpath(dp, root).split(os.sep)
        if parts and parts[0] in SIG_SKIP_DIRS:
            continue
        for f in sorted(fn):
            p = os.path.join(dp, f)
            h = hashlib.md5()
            try:
                with open(p, "rb") as fh:
                    for b in iter(lambda: fh.read(1 << 20), b""):
                        h.update(b)
            except OSError:
                continue
            ent[os.path.relpath(p, root)] = h.hexdigest()
    top = hashlib.md5()
    for r in sorted(ent):
        top.update(r.encode()); top.update(ent[r].encode())
    return top.hexdigest(), len(ent), ent


def sig_delta(before, after):
    b, a = before[2], after[2]
    added = sorted(set(a) - set(b))
    removed = sorted(set(b) - set(a))
    modified = sorted(k for k in (set(a) & set(b)) if a[k] != b[k])
    return added, removed, modified


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8,
                    help="fork-Pool workers for step52 (scheduling only; 14 deadlocked on 2026-08-03)")
    a = ap.parse_args()
    t0 = time.time()

    print("=" * 96, flush=True)
    print("STEP 83 — finish the MCC=1 ext2 chain (stages 3-4 only)", flush=True)
    print("=" * 96, flush=True)
    for n, v in _BINDINGS:
        print(f"  {n:<18} {v}", flush=True)

    need = [f"{S.MATDIR}/disease_Pv.npy", f"{S.MATDIR}/disease_Be.npy",
            f"{S.MATDIR}/disease_meta.tsv", f"{S.MATDIR}/disease_svorder.txt",
            f"{MCC1_OUT}/phewide_disease_full_firthrefit.tsv", S.PSUM]
    for p in need:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            raise SystemExit(f"ABORT: required stage-1/2 artefact missing or empty: {p}")
    print(f"\nreusing stage 1+2 artefacts ({len(need)} files, all present)", flush=True)

    H2.NWORK = a.workers
    print(f"step52 workers set to {H2.NWORK} (was 14 at the deadlock)", flush=True)

    pre = content_sig(PRIMARY_OUT)
    print(f"primary ext2 tree fingerprint (pre): {pre[0]} over {pre[1]} files", flush=True)

    print(f"\n[3/4] step52 approx Firth ({time.time()-t0:.0f}s)", flush=True)
    sys.argv = ["step52"]; H2.main()

    print(f"\n[4/4] step53 final Firth matrix ({time.time()-t0:.0f}s)", flush=True)
    sys.argv = ["step53"]; H3.main()

    post = content_sig(PRIMARY_OUT)
    print(f"\nprimary ext2 tree fingerprint (post): {post[0]} over {post[1]} files", flush=True)
    if post[0] != pre[0]:
        add, rem, mod = sig_delta(pre, post)
        for tag, lst in (("ADDED", add), ("REMOVED", rem), ("MODIFIED", mod)):
            for f_ in lst:
                print(f"  {tag:<9} {f_}", flush=True)
        raise SystemExit(f"ABORT: the primary ext2 tree CHANGED during the MCC=1 finish "
                         f"({len(add)} added / {len(rem)} removed / {len(mod)} modified — named above).")
    print("ISOLATION OK — primary ext2 tree byte-unchanged (0 files touched).", flush=True)
    print(f"\nMCC1_EXT2_CHAIN_DONE ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
