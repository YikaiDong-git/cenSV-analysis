#!/usr/bin/env python3
"""M5-ext2 — build the ISOLATED `finemap/code_ext2/` tree for the 11,129 arm."""
import os
import re
import sys
import shutil
import argparse

AR = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
      "disease_replication/aou_replication")
FM = f"{AR}/finemap"
SRC = f"{FM}/code"
DST = f"{FM}/code_ext2"

LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def patch(path, edits, label):
    with open(path) as f:
        s = f.read()
    for anchor, repl in edits:
        n = s.count(anchor)
        assert n == 1, f"{label}: anchor appears {n}x (need exactly 1): {anchor[:90]!r}"
        s = s.replace(anchor, repl)
    with open(path, "w") as f:
        f.write(s)
    say(f"  patched {label:<32} {len(edits)} edit(s)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="rebuild code_ext2/ from scratch")
    a = ap.parse_args()

    say("=" * 96)
    say("M5-ext2 — isolated code_ext2/ tree")
    say("=" * 96)

    if os.path.exists(DST):
        if not a.force:
            say(f"{DST} already exists — pass --force to rebuild")
            sys.exit(1)
        shutil.rmtree(DST)
    shutil.copytree(SRC, DST, ignore=shutil.ignore_patterns("__pycache__", "prep", "scripts"))
    say(f"copied {SRC} -> {DST}")

    patch(f"{DST}/finemap_core.py", [
        ('SVPGEN_DIR = f"{REPVDS}/step28_phewide/pgen"                    '
         '# svpgen300_{0..21}.{pgen,pvar,psam}',
         'SVPGEN_DIR = f"{REPVDS}/step28_phewide_ext2/pgen"               '
         '# EXT2: svpgen_<chrom>_b<NN>_000.*  (48 chunks)'),
        ('    fs = sorted(glob.glob(f"{SVPGEN_DIR}/svpgen300_*.pgen"),\n'
         '                key=lambda p: int(p.split("svpgen300_")[1].split(".pgen")[0]))',
         '    # EXT2: chunks are svpgen_<chrom>_b<NN>_000.pgen, not svpgen300_<i>.pgen. The old int() sort\n'
         '    # key raises ValueError on these names. Chunk ORDER is not semantically load-bearing --\n'
         '    # _svpgen_index() builds a dict keyed by sv_locus -- so a stable lexicographic sort is\n'
         '    # correct and avoids inventing a numeric ordering the filenames do not have.\n'
         '    fs = sorted(glob.glob(f"{SVPGEN_DIR}/svpgen_*.pgen"), key=os.path.basename)\n'
         '    assert fs, f"no ext2 svpgen chunks under {SVPGEN_DIR}"'),
        ('KEEP_DIR = f"{FM}/inputs/keep"', 'KEEP_DIR = f"{FM}/inputs_ext2/keep"      # EXT2 (incl. AFR)'),
        ('COVAR_DIR = f"{FM}/inputs/covar"', 'COVAR_DIR = f"{FM}/inputs_ext2/covar"    # EXT2 (incl. AFR)'),
        ('PHENO_DIR = f"{FM}/inputs/pheno"', 'PHENO_DIR = f"{FM}/inputs_ext2/pheno"    # EXT2 (290 phecodes)'),
        ('PRED_BFILE = CENTRES + "/vds_predictor_extract/{chrom}/{chrom}.predictor_snps"',
         'PRED_IDS_DIR = f"{FM}/inputs_ext2/pred_ids"        # EXT2: one SNP id per line, per chrom\n'
         '\n'
         '\n'
         'def load_pred_ids(chrom):\n'
         '    # EXT2: EN-predictor SNP ids for `chrom`, derived from the ext2 weight files (m6) rather\n'
         '    # than the stale vds_predictor_extract bim. FAIL-CLOSED: the ext2 arm must never silently\n'
         '    # fall back to an empty set -- an empty predictor set makes value_gate maximally permissive.\n'
         '    p = f"{PRED_IDS_DIR}/{chrom}.pred_ids.txt"\n'
         '    assert os.path.exists(p), f"missing ext2 EN-predictor ids: {p} (run m6_build_pred_ids_ext2.py)"\n'
         '    with open(p) as fh:\n'
         '        ids = {ln.strip() for ln in fh if ln.strip()}\n'
         '    assert ids, f"empty ext2 EN-predictor id set: {p}"\n'
         '    return ids'),
        ('    pred_ids = set()\n'
         '    pred_pref = PRED_BFILE.format(chrom=chrom)\n'
         '    if os.path.exists(pred_pref + ".bim"):\n'
         '        with open(pred_pref + ".bim") as fh:\n'
         '            for line in fh:\n'
         '                pred_ids.add(line.split("\\t")[1])',
         '    pred_ids = load_pred_ids(chrom)          # EXT2: derived set, not the stale bim'),
        ('    scratch = os.environ.get("FINEMAP_SCRATCH", "/tmp/finemap_scratch")',
         '    # EXT2: separate default. The extract is keyed by (stratum, chrom) alone and finemap_group\n'
         '    # deletes it after each group, so a shared root lets the two arms delete each other\'s\n'
         '    # in-use extracts and lets a stale extract be silently reused across arms.\n'
         '    scratch = os.environ.get("FINEMAP_SCRATCH", "/tmp/finemap_scratch_ext2")'),
        ('            "chrom": chrom, "phecode": phecode, "stratum": stratum, "mcc": mcc,',
         '            # EXT2 PROVENANCE. The published schema records no manifest, no dosage store and no\n'
         '            # weight-set, so an ext2 fit_meta.json would be indistinguishable in kind from a\n'
         '            # 6,513 one. Recorded here so the arm of any fit dir is answerable from the file.\n'
         '            "weight_set": "ext2_11129", "svpgen_dir": SVPGEN_DIR,\n'
         '            "keep_dir": KEEP_DIR, "covar_dir": COVAR_DIR, "pheno_dir": PHENO_DIR,\n'
         '            "chrom": chrom, "phecode": phecode, "stratum": stratum, "mcc": mcc,'),
    ], "finemap_core.py")

    patch(f"{DST}/finemap_group.py", [
        ('    ScoreNull, VREL, PRED_BFILE, KEEP_DIR,', '    ScoreNull, VREL, KEEP_DIR,'),
        ('def _load_pred_ids(chrom):\n'
         '    pred_ids = set()\n'
         '    pref = PRED_BFILE.format(chrom=chrom)\n'
         '    if os.path.exists(pref + ".bim"):\n'
         '        with open(pref + ".bim") as fh:\n'
         '            for line in fh:\n'
         '                pred_ids.add(line.split("\\t")[1])\n'
         '    return pred_ids',
         'def _load_pred_ids(chrom):\n'
         '    # EXT2: delegate to finemap_core.load_pred_ids (ext2-derived ids, fail-closed).\n'
         '    return core.load_pred_ids(chrom)'),
        ('    scratch = os.environ.get("FINEMAP_SCRATCH", "/tmp/finemap_scratch")',
         '    # EXT2: separate scratch root. The extract is keyed by (stratum, chrom) ONLY and this module\n'
         '    # DELETES it at group end, so a shared root lets the two arms remove each other\'s in-use\n'
         '    # extracts and lets a stale extract be silently reused across arms.\n'
         '    scratch = os.environ.get("FINEMAP_SCRATCH", "/tmp/finemap_scratch_ext2")'),
    ], "finemap_group.py")

    for fn, extra in (("run_finemap_group_batch.py", []),
                      ("run_finemap_batch.py",
                       [('BATCH_LOG = f"{core.FM}/results/_batch_status.tsv"',
                         'BATCH_LOG = f"{RESULTS}/_batch_status.tsv"   '
                         '# EXT2: was f"{core.FM}/results/..." -- a literal re-spelling of the results\n'
                         '# root that ignores RESULTS, so changing RESULTS alone still wrote to the '
                         'published tree.')]),
                      ("run_finemap_target.py", [])):
        p = f"{DST}/{fn}"
        with open(p) as f:
            s = f.read()
        edits = []
        if '{core.FM}/targets/target_manifest.tsv' in s:
            edits.append(('f"{core.FM}/targets/target_manifest.tsv"',
                          'f"{core.FM}/targets_ext2/target_manifest.tsv"'))
        if '{FM}/targets/target_manifest.tsv' in s:
            edits.append(('f"{FM}/targets/target_manifest.tsv"',
                          'f"{FM}/targets_ext2/target_manifest.tsv"'))
        if 'RESULTS = f"{core.FM}/results"' in s:
            edits.append(('RESULTS = f"{core.FM}/results"', 'RESULTS = f"{core.FM}/results_ext2"'))
        if 'RESULTS = f"{FM}/results"' in s:
            edits.append(('RESULTS = f"{FM}/results"', 'RESULTS = f"{FM}/results_ext2"'))
        edits += extra
        patch(p, edits, fn)

    patch(f"{DST}/collate_finemap.py", [
        ('FITS = os.path.join(FM, "results", "fits")\nOUT = os.path.join(FM, "results")',
         '# EXT2: the published version derives these from __file__, so a copy of this script placed\n'
         '# anywhere under finemap/ resolves them BACK to the published results/ and truncates\n'
         '# finemap_master.tsv / finemap_localization.tsv. Pinned explicitly instead.\n'
         'FITS = os.path.join(FM, "results_ext2", "fits")\n'
         'OUT = os.path.join(FM, "results_ext2")'),
    ], "collate_finemap.py")

    guard = '''"""_ext2_guard — fail-closed startup assertions for the isolated 11,129 arm.

Imported for side effects by every runner in this tree, BEFORE any fit runs. The failure signature of a
partial isolation fix is not a crash — it is silent "cached" skips and quietly-truncated dosages — so every
binding is asserted up front and the process aborts before the first write if any is wrong.
"""
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

    # every phecode a fit needs must be a column in the pheno file that fit will read
    man = getattr(caller, "MANIFEST", None)
    if man and os.path.exists(man):
        rows = list(csv.DictReader(open(man), delimiter="\\t"))
        for mcc in ("mcc1", "mcc2"):
            need = {r["phecode"] for r in rows if r["primary_mcc"].lower() == mcc}
            if not need:
                continue
            p = f"{core.PHENO_DIR}/disease.{mcc}.pheno"
            if not os.path.exists(p):
                fails.append(f"missing pheno {p}")
                continue
            with open(p) as f:
                cols = set(f.readline().rstrip("\\n").split("\\t"))
            miss = sorted(need - cols)
            req(not miss, f"{mcc}: {len(miss)} needed phecodes absent from {p}: {miss[:6]}")

    if fails:
        sys.stderr.write("\\n*** EXT2 ISOLATION GUARD FAILED ***\\n")
        for m in fails:
            sys.stderr.write(f"  - {m}\\n")
        sys.stderr.write("Aborting before any fit is written.\\n")
        sys.exit(3)
    print(f"[ext2-guard] OK  store={os.path.basename(os.path.dirname(core.SVPGEN_DIR))} "
          f"chunks={len(pfx)} inputs=inputs_ext2 results=results_ext2", flush=True)


check()
'''
    with open(f"{DST}/_ext2_guard.py", "w") as f:
        f.write(guard)
    say("  wrote    _ext2_guard.py                  (fail-closed startup assertions)")

    for fn in ("run_finemap_group_batch.py", "run_finemap_batch.py", "run_finemap_target.py"):
        p = f"{DST}/{fn}"
        with open(p) as f:
            s = f.read()
        m = re.search(r"^(FIT_R = .*)$", s, re.M)
        assert m, f"{fn}: no FIT_R line to anchor the guard import"
        s = s.replace(m.group(1), m.group(1) + "\nimport _ext2_guard  # noqa: F401,E402 "
                                              "- fail-closed isolation assertions (side effects on import)")
        with open(p, "w") as f:
            f.write(s)
    say("  wired    guard import into the 3 runners")

    say("")
    say("VERIFY")
    import ast as _ast
    import io as _io
    import tokenize as _tok

    def code_lines(path):
        with open(path) as f:
            src = f.read()
        lines = src.splitlines()
        drop = set()
        try:
            tree = _ast.parse(src)
            for node in _ast.walk(tree):
                body = getattr(node, "body", None)
                if not isinstance(body, list) or not body:
                    continue
                d = body[0]
                if (isinstance(d, _ast.Expr) and isinstance(d.value, _ast.Constant)
                        and isinstance(d.value.value, str)):
                    drop.update(range(d.lineno, (d.end_lineno or d.lineno) + 1))
        except SyntaxError:
            pass
        out = {i + 1: ln for i, ln in enumerate(lines) if (i + 1) not in drop}
        try:
            for t in _tok.generate_tokens(_io.StringIO(src).readline):
                if t.type == _tok.COMMENT and t.start[0] in out:
                    out[t.start[0]] = out[t.start[0]][:t.start[1]]
        except Exception:
            pass
        return out

    PATS = ((r'\{core\.FM\}/results(?!_ext2)', "results root not ext2"),
            (r'\{FM\}/results(?!_ext2)', "results root not ext2"),
            (r'targets/target_manifest', "manifest not ext2"),
            (r'inputs/(keep|covar|pheno)', "input dir not ext2"),
            (r'step28_phewide/pgen', "dosage store not ext2"),
            (r'svpgen300_', "old chunk glob"),
            (r'/tmp/finemap_scratch(?!_ext2)', "shared scratch root"),
            (r'vds_predictor_extract', "stale EN-predictor extract"))
    bad = []
    for root, _, files in os.walk(DST):
        for fn in sorted(files):
            if not fn.endswith(".py"):
                continue
            if fn == "_ext2_guard.py":
                continue
            for ln, txt in code_lines(os.path.join(root, fn)).items():
                for pat, why in PATS:
                    if re.search(pat, txt):
                        bad.append(f"{fn}:{ln}  {why}  ->  {txt.strip()[:80]}")
    if bad:
        say(f"  *** {len(bad)} residual published-tree reference(s):")
        for b in bad:
            say(f"      {b}")
        sys.exit(1)
    say("  0 residual published-tree references in code_ext2/")

    src_files = {f for f in os.listdir(SRC) if f.endswith((".py", ".R"))}
    dst_files = {f for f in os.listdir(DST) if f.endswith((".py", ".R"))}
    say(f"  files: {len(src_files)} source -> {len(dst_files)} in code_ext2 "
        f"(+_ext2_guard.py)")
    with open(f"{DST}/M5_BUILD_LOG.txt", "w") as f:
        f.write("\n".join(LOG) + "\n")
    say("M5_EXT2_DONE")


if __name__ == "__main__":
    main()
