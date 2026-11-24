#!/usr/bin/env python3
"""Step 79 (ext2 arm) — PER-ANCESTRY MCC=1 chain driver: step56 (score) -> step57 (Firth)."""
import os, sys, time, argparse, hashlib

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S     # noqa: E402  — engine; patched below BEFORE any other import

RV = S.RV
MCC1_INP = f"{RV}/inputs_mcc1_ext2"
MCC1_OUT = f"{RV}/step30_phewide_mcc1_ext2"
MCC2_OUT = f"{RV}/step30_phewide_ext2"

S.INP = MCC1_INP
S.OUTA = MCC1_OUT
S.MATDIR = f"{MCC1_OUT}/disease_matrices"
S.PSUM = f"{MCC1_OUT}/phecodeX_pheno_summary.tsv"
assert S.INP == MCC1_INP and S.OUTA == MCC1_OUT, "engine patch did not take"

import step56_ancestry_disease_ext2 as P     # noqa: E402  binds MATDIR at import — after the patch
import step57_perancestry_firth_ext2 as F    # noqa: E402  binds FIRTH_MATDIR at import — after the patch

_BINDINGS = [
    ("step51.INP", S.INP), ("step51.OUTA", S.OUTA), ("step51.MATDIR", S.MATDIR), ("step51.PSUM", S.PSUM),
    ("step56.MATDIR", P.MATDIR), ("step57.FIRTH_MATDIR", F.FIRTH_MATDIR),
]
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
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--only", default=None, help="comma-separated ancestries (default: all 6)")
    a = ap.parse_args()
    t0 = time.time()

    print("=" * 96, flush=True)
    print("STEP 79 — per-ancestry MCC=1 chain (ext2 weights)", flush=True)
    print("=" * 96, flush=True)
    for n, v in _BINDINGS:
        print(f"  {n:<22} {v}", flush=True)

    for f in ("disease.pheno", "covar.txt"):
        p = f"{MCC1_INP}/{f}"
        if not os.path.exists(p) or os.path.getsize(p) < (1 << 20):
            raise SystemExit(f"ABORT: missing/empty MCC=1 input {p}")
    if not os.path.exists(S.PSUM):
        raise SystemExit(f"ABORT: missing MCC=1 phecode summary {S.PSUM}")
    nph = len(open(f"{MCC1_INP}/disease.pheno").readline().rstrip("\n").split("\t")) - 2
    print(f"\nMCC=1 disease.pheno carries {nph} phecode columns", flush=True)

    before = content_sig(MCC2_OUT)
    print(f"MCC=2 tree content_sig BEFORE: {before[0]}  ({before[1]} files)", flush=True)

    print(f"\n[stage 1/2] step56 per-ancestry score test ({time.time()-t0:.0f}s)", flush=True)
    sys.argv = ["step56", "--batch", str(a.batch)] + (["--only", a.only] if a.only else [])
    P.main()

    print(f"\n[stage 2/2] step57 per-ancestry Firth ({time.time()-t0:.0f}s)", flush=True)
    sys.argv = ["step57"] + (["--only", a.only] if a.only else [])
    F.main()

    after = content_sig(MCC2_OUT)
    print(f"\nMCC=2 tree content_sig AFTER : {after[0]}  ({after[1]} files)", flush=True)
    if after[0] != before[0]:
        add, rem, mod = sig_delta(before, after)
        for tag, lst in (("ADDED", add), ("REMOVED", rem), ("MODIFIED", mod)):
            for f_ in lst:
                print(f"  {tag:<9} {f_}", flush=True)
        raise SystemExit(f"ABORT: the MCC=2 tree CHANGED during the MCC=1 per-ancestry run "
                         f"({len(add)} added / {len(rem)} removed / {len(mod)} modified — named above). "
                         f"Isolation violated; investigate before trusting either arm.")
    print("ISOLATION OK — MCC=2 tree byte-unchanged (0 files touched).", flush=True)
    print(f"MCC1_PERANCESTRY_DONE ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
