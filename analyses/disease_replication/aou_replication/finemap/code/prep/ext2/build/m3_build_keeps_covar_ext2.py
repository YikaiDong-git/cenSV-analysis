#!/usr/bin/env python3
"""M3-ext2 — per-stratum finemap keeps + covars for the 11,129 arm, as a COMPLETE ISOLATED input set."""
import os
import csv
import sys
import hashlib
import argparse

AR = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
      "disease_replication/aou_replication")
FAM = f"{AR}/results_vds/phewide_carryover/sv.fam"
COV3 = os.path.expanduser("~/workspace/working/scripts/covariates/output/covariates_v3_all_ancestries.tsv")
COVBUILT = f"{AR}/results_vds/phewide_carryover/inputs/covar.txt"
REL = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/"
                         "aux/relatedness/samples_relatedness_flagged_samples.tsv")

PUB_K, PUB_C = f"{AR}/finemap/inputs/keep", f"{AR}/finemap/inputs/covar"
OUTK, OUTC = f"{AR}/finemap/inputs_ext2/keep", f"{AR}/finemap/inputs_ext2/covar"

STRATA = ["all-anc", "EUR", "AMR", "AFR"]
ANC_OF = {"EUR": "eur", "AMR": "amr", "AFR": "afr"}
PUBLISHED_STRATA = {"all-anc", "EUR", "AMR"}

LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def guard(path):
    assert not os.path.islink(path), f"refusing to write through symlink: {path}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-gate", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUTK, exist_ok=True)
    os.makedirs(OUTC, exist_ok=True)

    say("=" * 96)
    say("M3-ext2 — isolated per-stratum keeps + covars (adds AFR)")
    say("=" * 96)

    cohort = []
    with open(FAM) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                cohort.append(p[1])
    cohort_set = set(cohort)
    say(f"cohort (sv.fam): {len(cohort):,}")

    anc = {}
    with open(COV3) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            iid = row["IID"]
            if iid in cohort_set:
                anc[iid] = (row.get("ancestry_pred") or "").strip().lower()
    say(f"ancestry_pred resolved for {len(anc):,} cohort samples")

    drop = set()
    with open(REL) as f:
        next(f)
        for line in f:
            drop.add(line.strip())
    say(f"relatedness-flagged: {len(drop):,}")

    def in_stratum(iid, s):
        return True if s == "all-anc" else anc.get(iid) == ANC_OF[s]

    summary = []
    memb = {}
    for s in STRATA:
        full = [i for i in cohort if in_stratum(i, s)]
        unrel = [i for i in full if i not in drop]
        memb[s] = set(full)
        for nm, rows in ((f"{s}.keep.fid_iid.tsv", [f"{i}\t{i}\n" for i in unrel]),
                         (f"{s}.keep.iid.txt", [f"{i}\n" for i in unrel]),
                         (f"{s}.keep_full.fid_iid.tsv", [f"{i}\t{i}\n" for i in full])):
            p = f"{OUTK}/{nm}"
            guard(p)
            with open(p, "w") as f:
                f.writelines(rows)
        summary.append((s, len(full), len(unrel), len(full) - len(unrel)))

    with open(COVBUILT) as f:
        header = f.readline().rstrip("\n").split("\t")
        anc_cols = {i for i, c in enumerate(header) if c.startswith("anc_")}
        keep_idx = [i for i in range(len(header)) if i not in anc_cols]
        iid_col = header.index("IID")
        fh = {}
        for s in STRATA:
            p = f"{OUTC}/{s}.covar.txt"
            guard(p)
            fh[s] = open(p, "w")
            cols = header if s == "all-anc" else [header[i] for i in keep_idx]
            fh[s].write("\t".join(cols) + "\n")
        for line in f:
            p = line.rstrip("\n").split("\t")
            iid = p[iid_col]
            for s in STRATA:
                if iid in memb[s]:
                    fh[s].write((line if line.endswith("\n") else line + "\n") if s == "all-anc"
                                else "\t".join(p[i] for i in keep_idx) + "\n")
        for h in fh.values():
            h.close()
        with open(f"{OUTC}/covar_names.peranc.txt", "w") as f2:
            f2.write(",".join(c for c in header[2:] if not c.startswith("anc_")))
        with open(f"{OUTC}/covar_names.allanc.txt", "w") as f2:
            f2.write(",".join(header[2:]))

    say("")
    say(f"{'stratum':<10}{'N_full':>10}{'N_unrel':>10}{'N_related_dropped':>20}")
    for s, nf, nu, nd in summary:
        tag = "  <-- NEW in ext2" if s not in PUBLISHED_STRATA else ""
        say(f"{s:<10}{nf:>10,}{nu:>10,}{nd:>20,}{tag}")

    if not a.no_gate:
        say("")
        say("GATE — the three published strata must rebuild byte-identically")
        bad = []
        checked = 0
        for s in sorted(PUBLISHED_STRATA):
            for nm, pubdir, outdir in ((f"{s}.keep.fid_iid.tsv", PUB_K, OUTK),
                                       (f"{s}.keep.iid.txt", PUB_K, OUTK),
                                       (f"{s}.keep_full.fid_iid.tsv", PUB_K, OUTK),
                                       (f"{s}.covar.txt", PUB_C, OUTC)):
                pub, new = f"{pubdir}/{nm}", f"{outdir}/{nm}"
                if not os.path.exists(pub):
                    say(f"  {nm:<32} no published copy — SKIPPED")
                    continue
                checked += 1
                if md5(pub) == md5(new):
                    say(f"  {nm:<32} PASS  {md5(new)}")
                else:
                    say(f"  {nm:<32} *** FAIL ***  pub {md5(pub)} vs new {md5(new)}")
                    bad.append(nm)
        say(f"  {checked - len(bad)}/{checked} byte-identical")
        if bad:
            say("\nM3_EXT2_GATE_FAILED — the port does not reproduce the published inputs; "
                "do not build AFR on top of an unverified rebuild.")
            sys.exit(1)

    with open(f"{AR}/finemap/inputs_ext2/m3_ext2_summary.txt", "w") as f:
        f.write("\n".join(LOG) + "\n")
    say("")
    say(f"wrote {OUTK}/  and  {OUTC}/")
    say("M3_EXT2_DONE")


if __name__ == "__main__":
    main()
