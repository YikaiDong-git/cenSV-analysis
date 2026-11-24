#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--master-manifest", type=Path, required=True)
    ap.add_argument("--sweep-root", type=Path, required=True,
                    help="dir containing <pheno>/<chrom>/<sv_locus>__<arm>__diag.tsv")
    ap.add_argument("--pheno", required=True,
                    help="phenotype subdir under sweep-root")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    pheno_dir = args.sweep_root / args.pheno
    pheno_dir.mkdir(parents=True, exist_ok=True)

    existing = set()
    for diag in pheno_dir.rglob("*__diag.tsv"):
        stem = diag.name[:-len("__diag.tsv")]
        if stem.endswith("__best8"):
            arm = "best8"
            sv_locus = stem[:-len("__best8")]
        elif stem.endswith("__best6_alt"):
            arm = "best6_alt"
            sv_locus = stem[:-len("__best6_alt")]
        else:
            raise ValueError(f"unknown arm token in diag name: {diag.name}")
        existing.add((sv_locus, arm))

    sys.stderr.write(
        f"[diag] existing files for pheno={args.pheno}: {len(existing)}\n"
    )

    out_rows = []
    out_new_pair_id = 0
    n_master = 0
    with open(args.master_manifest) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        sv_locus_idx = header.index("sv_locus")
        arm_idx = header.index("arm")
        pair_id_idx = header.index("pair_id")
        for line in fh:
            n_master += 1
            parts = line.rstrip("\n").split("\t")
            sv_locus = parts[sv_locus_idx]
            arm = parts[arm_idx]
            if (sv_locus, arm) in existing:
                continue
            out_new_pair_id += 1
            master_pair_id = parts[pair_id_idx]
            new_row = [str(out_new_pair_id)] + parts[1:] + [master_pair_id]
            out_rows.append(new_row)

    sys.stderr.write(
        f"[manifest] master rows={n_master}  missing for {args.pheno}={len(out_rows)}\n"
    )
    if not out_rows:
        sys.stderr.write(f"[done] no missing rows for {args.pheno}\n")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(args.out.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        fh.write("\t".join(header + ["master_pair_id"]) + "\n")
        for row in out_rows:
            fh.write("\t".join(row) + "\n")
    os.replace(tmp, args.out)
    sys.stderr.write(f"[done] wrote {args.out} ({len(out_rows)} rows)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
