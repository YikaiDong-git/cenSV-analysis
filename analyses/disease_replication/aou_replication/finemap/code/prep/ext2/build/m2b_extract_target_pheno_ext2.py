#!/usr/bin/env python3
"""M2b-ext2 — compact target-phecode phenotypes for the ext2 fine-mapping arm."""
import os
import csv
import sys
import subprocess

AR = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
      "disease_replication/aou_replication")
MAN = f"{AR}/finemap/targets_ext2/target_manifest.tsv"
SRC = {"mcc1": f"{AR}/results_vds/inputs_mcc1_ext2/disease.pheno",
       "mcc2": f"{AR}/results_vds/phewide_carryover/inputs/disease.pheno"}
OUT = f"{AR}/finemap/inputs_ext2/pheno"
PUB_OUT = f"{AR}/finemap/inputs/pheno"

LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    assert os.path.realpath(OUT) != os.path.realpath(PUB_OUT), \
        "OUT resolves to the PUBLISHED pheno directory — refusing to run"
    os.makedirs(OUT, exist_ok=True)
    assert not os.path.islink(OUT), f"refusing to write into symlinked dir: {OUT}"

    rows = list(csv.DictReader(open(MAN), delimiter="\t"))
    targets = sorted({r["phecode"] for r in rows})
    need = {m: sorted({r["phecode"] for r in rows if r["primary_mcc"].lower() == m}) for m in SRC}
    say("=" * 96)
    say("M2b-ext2 — target-phecode phenotypes for the ext2 fine-mapping arm")
    say("=" * 96)
    say(f"manifest: {len(rows)} fits · {len(targets)} distinct phecodes")
    for m in SRC:
        say(f"  fits with primary_mcc={m}: {sum(1 for r in rows if r['primary_mcc'].lower() == m):>4}"
            f"   phecodes they need: {len(need[m])}")

    fatal = []
    for mcc, src in SRC.items():
        assert os.path.exists(src), f"missing source for {mcc}: {src}"
        with open(src) as f:
            hdr = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        present = [p for p in targets if p in idx]
        missing_all = [p for p in targets if p not in idx]
        missing_needed = [p for p in need[mcc] if p not in idx]

        cols = [idx["FID"] + 1, idx["IID"] + 1] + [idx[p] + 1 for p in present]
        out = f"{OUT}/disease.{mcc}.pheno"
        assert not os.path.islink(out), f"refusing to write through symlink: {out}"
        with open(out, "w") as fo:
            subprocess.run(["cut", "-f", ",".join(map(str, cols)), src], stdout=fo, check=True)

        say("")
        say(f"{mcc}: source {os.path.basename(os.path.dirname(src))}/{os.path.basename(src)} "
            f"({len(hdr):,} cols)")
        say(f"   wrote {len(present)}/{len(targets)} of ALL manifest phecodes -> {out}")
        say(f"   coverage of the phecodes THIS file must serve: "
            f"{len(need[mcc]) - len(missing_needed)}/{len(need[mcc])}")
        if missing_all:
            say(f"   absent from source but NOT needed here: {len(missing_all) - len(missing_needed)}"
                + (f"  e.g. {missing_all[:5]}" if missing_all else ""))
        if missing_needed:
            say(f"   *** MISSING AND NEEDED: {len(missing_needed)} -> {missing_needed}")
            fatal += [(mcc, p) for p in missing_needed]

    say("")
    if fatal:
        say(f"M2b_EXT2_GATE_FAILED — {len(fatal)} (mcc, phecode) pairs are required by a fit but absent "
            f"from its source pheno.")
        say("A missing phecode raises KeyError inside build_group's unguarded loop and takes down the whole")
        say("(chrom, stratum) group, so this must be resolved before any fit runs.")
        with open(f"{OUT}/m2b_ext2_summary.txt", "w") as f:
            f.write("\n".join(LOG) + "\n")
        sys.exit(1)

    say("COVERAGE GATE PASS — every phecode any fit needs is present in the file that fit will read.")
    for mcc in SRC:
        p = f"{OUT}/disease.{mcc}.pheno"
        with open(p) as f:
            nc = len(f.readline().rstrip("\n").split("\t"))
        say(f"  disease.{mcc}.pheno  {os.path.getsize(p) / 1e6:.1f} MB  {nc:,} cols")
    with open(f"{OUT}/m2b_ext2_summary.txt", "w") as f:
        f.write("\n".join(LOG) + "\n")
    say("M2B_EXT2_DONE")


if __name__ == "__main__":
    main()
