"""_ext2_guard — fail-closed startup assertions for the isolated 11,129 arm."""
import os
import csv
import sys

import finemap_core as core


def check():
    fails = []

    def req(cond, msg):
        if not cond:
            fails.append(msg)

    pub = os.path.realpath(f"{core.FM}/results")
    for name, val in (("KEEP_DIR", core.KEEP_DIR), ("COVAR_DIR", core.COVAR_DIR),
                      ("PHENO_DIR", core.PHENO_DIR)):
        req("inputs_ext2" in val, f"{name} is not the ext2 input tree: {val}")
        req(os.path.isdir(val), f"{name} does not exist: {val}")

    req("step28_phewide_ext2" in core.SVPGEN_DIR,
        f"SVPGEN_DIR is not the ext2 dosage store: {core.SVPGEN_DIR}")
    pfx = core._svpgen_prefixes()
    req(len(pfx) == 48, f"expected 48 ext2 svpgen chunks, found {len(pfx)}")

    mod = sys.modules[__name__]
    here = os.path.dirname(os.path.abspath(mod.__file__))
    req(os.path.dirname(os.path.abspath(core.__file__)) == here,
        f"finemap_core was imported from OUTSIDE this tree: {core.__file__}")

    caller = sys.modules.get("__main__")
    for attr, want in (("RESULTS", "results_ext2"), ("MANIFEST", "targets_ext2")):
        val = getattr(caller, attr, None)
        if val is None:
            continue
        req(want in val, f"{attr} is not ext2-scoped: {val}")
        if attr == "RESULTS":
            req(not os.path.realpath(val).startswith(pub + os.sep)
                and os.path.realpath(val) != pub,
                f"RESULTS resolves inside the PUBLISHED tree: {val}")
        if attr == "MANIFEST":
            req(os.path.exists(val), f"MANIFEST does not exist: {val}")

    man = getattr(caller, "MANIFEST", None)
    if man and os.path.exists(man):
        rows = list(csv.DictReader(open(man), delimiter="\t"))
        for mcc in ("mcc1", "mcc2"):
            need = {r["phecode"] for r in rows if r["primary_mcc"].lower() == mcc}
            if not need:
                continue
            p = f"{core.PHENO_DIR}/disease.{mcc}.pheno"
            if not os.path.exists(p):
                fails.append(f"missing pheno {p}")
                continue
            with open(p) as f:
                cols = set(f.readline().rstrip("\n").split("\t"))
            miss = sorted(need - cols)
            req(not miss, f"{mcc}: {len(miss)} needed phecodes absent from {p}: {miss[:6]}")

    if fails:
        sys.stderr.write("\n*** EXT2 ISOLATION GUARD FAILED ***\n")
        for m in fails:
            sys.stderr.write(f"  - {m}\n")
        sys.stderr.write("Aborting before any fit is written.\n")
        sys.exit(3)
    print(f"[ext2-guard] OK  store={os.path.basename(os.path.dirname(core.SVPGEN_DIR))} "
          f"chunks={len(pfx)} inputs=inputs_ext2 results=results_ext2", flush=True)


check()
