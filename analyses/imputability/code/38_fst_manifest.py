#!/usr/bin/env python3
import argparse
import hashlib
import json
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
HERE = ROOT / "analyses" / "imputability"
OUT_V3 = HERE / "output" / "v3"
MANIFEST = OUT_V3 / "FST_MANIFEST.json"

SUFFIX = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000"
          "_qgap200_rgap200")
SV_BASE = ROOT / "sv_extraction" / "decompose" / SUFFIX
PANEL = ROOT / "resources" / "1000G" / "20130606_g1k_3202_samples_ped_population.txt"

POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

PRODUCER_SV = HERE / "code" / "13b_sv_ancestry.py"
PRODUCER_WINDOW = HERE / "code" / "13_ancestry_features_per_chrom.py"

TABLES = [
    ("output/features/sv_ancestry_features.tsv", ["max_FST_pair_sv"]),
    ("output/features/per_sv_features_with_genome.tsv",
     ["max_FST_pair_sv", "max_FST_pair_window"]),
]
SUPERSEDED = [
    "output/features/_archive/sv_ancestry_features.tsv.individuals_convention_20260907T002007",
    "output/features/_archive/per_sv_features_with_genome.tsv"
    ".individuals_convention_20260907T002007",
]

CONVENTIONS = {"gene_copies": 2, "individuals": 1}
TOL = 1e-12


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_pop_map():
    df = pd.read_csv(PANEL, sep=r"\s+", engine="python")
    return dict(zip(df["SampleID"], df["Superpopulation"]))


def n_per_pop_by_chrom(pop_map):
    out = {}
    for chrom in CHROMS:
        sam = SV_BASE / chrom / f"{chrom}_sv_samples.txt"
        if not sam.exists():
            continue
        samples = sam.read_text().splitlines()
        assigned = [pop_map.get(s, "OTHER") for s in samples]
        out[chrom] = {pop: int(assigned.count(pop)) for pop in POPS}
    return out


def recompute_max_fst(af, n_per_pop, multiplier):
    stack = []
    for i, j in combinations(POPS, 2):
        pi = af[i]
        pj = af[j]
        ni = max(int(multiplier * n_per_pop[i]), 2)
        nj = max(int(multiplier * n_per_pop[j]), 2)
        num = (pi - pj) ** 2 - pi * (1 - pi) / (ni - 1) - pj * (1 - pj) / (nj - 1)
        den = pi * (1 - pj) + pj * (1 - pi)
        with np.errstate(invalid="ignore", divide="ignore"):
            stack.append(np.where(den > 0, num / den, np.nan))
    arr = np.vstack(stack)
    allnan = np.all(np.isnan(arr), axis=0)
    with np.errstate(invalid="ignore"):
        mx = np.nanmax(np.where(np.isnan(arr), -np.inf, arr), axis=0)
    return np.where(allnan, np.nan, mx)


def reproduce(table_path, n_by_chrom):
    cols = ["chrom", "sv_index", "max_FST_pair_sv"] + [f"af_{p}_sv" for p in POPS]
    df = pd.read_csv(table_path, sep="\t", usecols=cols, low_memory=False)
    result = {name: {"n_compared": 0, "n_within_tol": 0, "max_abs_diff": 0.0}
              for name in CONVENTIONS}

    for chrom, sub in df.groupby("chrom", sort=False):
        if chrom not in n_by_chrom:
            sys.exit(f"FAIL no sample list for {chrom}; cannot derive n_per_pop")
        af = {p: sub[f"af_{p}_sv"].to_numpy(dtype=float) for p in POPS}
        shipped = sub["max_FST_pair_sv"].to_numpy(dtype=float)
        keep = np.isfinite(shipped)
        for name, mult in CONVENTIONS.items():
            got = recompute_max_fst(af, n_by_chrom[chrom], mult)
            d = np.abs(got[keep] - shipped[keep])
            d = d[np.isfinite(d)]
            r = result[name]
            r["n_compared"] += int(keep.sum())
            r["n_within_tol"] += int((d <= TOL).sum())
            r["max_abs_diff"] = max(r["max_abs_diff"], float(d.max()) if d.size else 0.0)
    return result


def window_source_is_uncorrected():
    src = PRODUCER_WINDOW.read_text()
    return bool(re.search(r"ni\s*=\s*np\.maximum\(n_eff\[i\],\s*1\)", src))


def verdict(rep):
    matched = [k for k, v in rep.items()
               if v["n_compared"] and v["n_within_tol"] == v["n_compared"]]
    if len(matched) != 1:
        sys.exit("FAIL the reproduction does not discriminate: "
                 f"{len(matched)} of {len(rep)} conventions reproduce the column "
                 f"({matched}). A check that cannot fail proves nothing.")
    return matched[0]


def build(stamp):
    pop_map = load_pop_map()
    n_by_chrom = n_per_pop_by_chrom(pop_map)

    files = []
    for rel, fst_cols in TABLES:
        p = HERE / rel
        if not p.exists():
            sys.exit(f"FAIL missing table {p}")
        files.append({"path": rel, "bytes": p.stat().st_size,
                      "sha256": sha256(p), "fst_columns": fst_cols})

    superseded = []
    for rel in SUPERSEDED:
        p = HERE / rel
        if not p.exists():
            sys.exit(f"FAIL missing archived original {p}")
        superseded.append({"path": rel, "bytes": p.stat().st_size,
                           "sha256": sha256(p)})

    rep = reproduce(HERE / TABLES[0][0], n_by_chrom)
    held = verdict(rep)

    cross = cross_file(n_by_chrom)

    return {
        "stamp": stamp,
        "what": ("Which Hudson F_ST convention the v3 feature tables hold. The tables "
                 "are 24 MB each and are not tracked; this manifest is."),
        "columns": {
            "max_FST_pair_sv": {
                "producer": str(PRODUCER_SV.relative_to(ROOT)),
                "n_is": "gene_copies",
                "formula": ("N_ij = (p_i-p_j)^2 - p_i(1-p_i)/(2n_i-1) "
                            "- p_j(1-p_j)/(2n_j-1); D_ij = p_i(1-p_j) + p_j(1-p_i)"),
                "corrected_in": "3e67323f",
                "evidence": ("recomputed from the table's own af_*_sv columns under "
                             "both conventions"),
                "held_by_file": held,
            },
            "max_FST_pair_window": {
                "producer": str(PRODUCER_WINDOW.relative_to(ROOT)),
                "n_is": "individuals",
                "corrected_in": None,
                "evidence": ("asserted from the producing source at L108-110, not "
                             "recomputed: the SNP genotype matrices are needed"),
                "source_still_uncorrected": window_source_is_uncorrected(),
                "scope": ("a gradient-boosting predictor feature; no reported F_ST "
                          "number derives from it"),
            },
        },
        "files": files,
        "superseded": superseded,
        "n_per_pop": n_by_chrom,
        "reproduction": rep,
        "cross_file": cross,
        "tolerance": TOL,
    }


def cross_file(n_by_chrom):
    a = pd.read_csv(HERE / TABLES[0][0], sep="\t",
                    usecols=["chrom", "sv_index", "max_FST_pair_sv"], low_memory=False)
    b = pd.read_csv(HERE / TABLES[1][0], sep="\t",
                    usecols=["chrom", "sv_index", "max_FST_pair_sv"], low_memory=False)
    m = a.merge(b, on=["chrom", "sv_index"], suffixes=("_src", "_merged"))
    d = (m["max_FST_pair_sv_src"] - m["max_FST_pair_sv_merged"]).abs()
    both_nan = m["max_FST_pair_sv_src"].isna() & m["max_FST_pair_sv_merged"].isna()
    ok = int(((d <= TOL) | both_nan).sum())
    return {"n_joined": int(len(m)), "n_agreeing": ok,
            "max_abs_diff": float(d.max(skipna=True)) if len(m) else 0.0}


def report(man):
    print("stamp: %s" % man["stamp"])
    for f in man["files"]:
        print("  %-46s %10d B  %s" % (f["path"].split("/")[-1], f["bytes"],
                                      f["sha256"][:16]))
    print("\nmax_FST_pair_sv, recomputed from the af_*_sv columns:")
    for name, r in man["reproduction"].items():
        print("  %-12s  %7d of %7d within %g   max |diff| = %.3g"
              % (name, r["n_within_tol"], r["n_compared"], man["tolerance"],
                 r["max_abs_diff"]))
    print("  the file holds: %s"
          % man["columns"]["max_FST_pair_sv"]["held_by_file"])
    c = man["cross_file"]
    print("\ncross-file: %d of %d joined rows agree" % (c["n_agreeing"], c["n_joined"]))
    w = man["columns"]["max_FST_pair_window"]
    print("max_FST_pair_window: n_is=%s, source still uncorrected: %s"
          % (w["n_is"], w["source_still_uncorrected"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--write", action="store_true")
    g.add_argument("--check", action="store_true")
    ap.add_argument("--stamp", help="required with --write, e.g. 2026-09-07")
    args = ap.parse_args()

    if args.write:
        if not args.stamp:
            sys.exit("FAIL --write needs --stamp so a rebuild is reproducible from "
                     "the log")
        man = build(args.stamp)
        report(man)
        MANIFEST.write_text(json.dumps(man, indent=2, sort_keys=True) + "\n")
        print("\nwrote %s" % MANIFEST)
        return

    if not MANIFEST.exists():
        sys.exit(f"FAIL no manifest at {MANIFEST}; run --write first")
    old = json.loads(MANIFEST.read_text())
    new = build(old["stamp"])
    report(new)

    bad = []
    for a, b in zip(old["files"], new["files"]):
        if a["sha256"] != b["sha256"]:
            bad.append("%s: sha256 %s -> %s" % (a["path"], a["sha256"][:16],
                                                b["sha256"][:16]))
    for a, b in zip(old["superseded"], new["superseded"]):
        if a["sha256"] != b["sha256"]:
            bad.append("archived %s changed" % a["path"])

    held_old = old["columns"]["max_FST_pair_sv"]["held_by_file"]
    held_new = new["columns"]["max_FST_pair_sv"]["held_by_file"]
    if held_old != held_new:
        bad.append("max_FST_pair_sv convention %s -> %s" % (held_old, held_new))
    if held_new != "gene_copies":
        bad.append("max_FST_pair_sv holds %s; the corrected convention is gene_copies"
                   % held_new)

    c = new["cross_file"]
    if c["n_agreeing"] != c["n_joined"]:
        bad.append("cross-file: %d of %d rows disagree"
                   % (c["n_joined"] - c["n_agreeing"], c["n_joined"]))

    w_old = old["columns"]["max_FST_pair_window"]["source_still_uncorrected"]
    w_new = new["columns"]["max_FST_pair_window"]["source_still_uncorrected"]
    if w_old != w_new:
        bad.append("13_ancestry_features_per_chrom.py changed: the window column must "
                   "be re-derived and this manifest rewritten")

    print()
    if bad:
        for line in bad:
            print("FAIL %s" % line)
        sys.exit(1)
    print("PASS every check reproduced")


if __name__ == "__main__":
    main()
