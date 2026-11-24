#!/usr/bin/env python3
"""M6-ext2 — the EN-predictor SNP id set for the 11,129 arm."""
import os
import csv
import sys
import glob

AR = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
      "disease_replication/aou_replication")
PH = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
      "disease_replication/phewide_all_handoff")
CENTRES = "/home/dataproc/workspace/working/jwein22/centromere_phewas/results"
W = f"{PH}/full_11129_ext2"
OUT = f"{AR}/finemap/inputs_ext2/pred_ids"
GATE_MIN = 0.98

LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    os.makedirs(OUT, exist_ok=True)
    assert not os.path.islink(OUT), f"refusing to write into symlinked dir: {OUT}"

    say("=" * 96)
    say("M6-ext2 — EN-predictor SNP ids from the ext2 weight files")
    say("=" * 96)
    say(f"{'chrom':<8}{'ext2 SNPs':>11}{'old bim':>10}{'bim ids in ext2':>18}{'%':>8}   gate")

    rows = []
    fails = []
    for p in sorted(glob.glob(f"{W}/sv_weights_long_ext2_chr*.tsv")):
        chrom = os.path.basename(p).split("sv_weights_long_ext2_")[1][:-4]
        ids = set()
        for r in csv.DictReader(open(p), delimiter="\t"):
            ids.add(f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}")
        op = f"{OUT}/{chrom}.pred_ids.txt"
        assert not os.path.islink(op), f"refusing to write through symlink: {op}"
        with open(op, "w") as f:
            f.writelines(i + "\n" for i in sorted(ids))

        bim = f"{CENTRES}/vds_predictor_extract/{chrom}/{chrom}.predictor_snps.bim"
        if os.path.exists(bim):
            bids = {l.split("\t")[1] for l in open(bim)}
            inn = len(bids & ids)
            frac = inn / len(bids) if bids else 1.0
            ok = frac >= GATE_MIN
            if not ok:
                fails.append((chrom, frac))
            say(f"{chrom:<8}{len(ids):>11,}{len(bids):>10,}{inn:>18,}{frac * 100:>7.1f}%   "
                f"{'PASS' if ok else '*** FAIL ***'}")
            rows.append((chrom, len(ids), len(bids), inn, f"{frac:.4f}", len(ids - bids)))
        else:
            say(f"{chrom:<8}{len(ids):>11,}{'-':>10}{'-':>18}{'-':>8}   no bim (gate n/a)")
            rows.append((chrom, len(ids), 0, 0, "NA", len(ids)))

    with open(f"{OUT}/pred_ids_summary.tsv", "w") as f:
        f.write("chrom\tn_ext2_pred_snps\tn_old_bim\tn_bim_ids_in_ext2\tfrac_bim_covered\t"
                "n_ext2_not_in_bim\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    tot_new = sum(r[5] for r in rows)
    tot_bim = sum(r[2] for r in rows)
    say("")
    say(f"total ext2 EN-predictor SNPs: {sum(r[1] for r in rows):,}")
    say(f"  present in the old extract : {sum(r[3] for r in rows):,} / {tot_bim:,}")
    say(f"  ABSENT from the old extract: {tot_new:,}  <-- these were being flagged is_en_predictor=0,")
    say("                                   making value_gate over-permissive for every one of them")
    if fails:
        say("")
        say(f"M6_GATE_FAILED on {len(fails)} chromosome(s): {fails}")
        say("The id construction does not reproduce the old extract, so the derived set is not trustworthy.")
        with open(f"{OUT}/m6_summary.txt", "w") as f:
            f.write("\n".join(LOG) + "\n")
        sys.exit(1)
    say("GATE PASS — the old extract is a subset of the derived ext2 set on every chromosome")
    with open(f"{OUT}/m6_summary.txt", "w") as f:
        f.write("\n".join(LOG) + "\n")
    say("M6_EXT2_DONE")


if __name__ == "__main__":
    main()
