#!/usr/bin/env python3
"""Build the finemap manifest for the full-restart sweep."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]


BASE = Path(
    CENSV_ROOT + "/"
    "analyses/finemap_multi_sv"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region-map", type=Path, default=BASE / "inputs/region_map_all_chroms.tsv")
    p.add_argument("--pheno-list", type=Path, default=BASE / "inputs/pheno_top30.tsv")
    p.add_argument("--admit-v2", type=Path,
                   default=BASE / "output/full_restart/admit/admit_decompose_cvr2_0p3.tsv")
    p.add_argument("--out", type=Path,
                   default=BASE / "output/full_restart/manifests/finemap_manifest_v2.tsv")
    p.add_argument("--out-region-map", type=Path,
                   default=BASE / "output/full_restart/manifests/region_map_v2_with_admit.tsv")
    return p.parse_args()


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.parent / f"{path.stem}.tmp.{os.getpid()}{path.suffix}"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def main() -> int:
    args = parse_args()

    region_map = pd.read_csv(args.region_map, sep="\t")
    admit = pd.read_csv(args.admit_v2, sep="\t")
    pheno = pd.read_csv(args.pheno_list, sep="\t")

    sys.stderr.write(f"[manifest-v2] region map: {len(region_map)} regions\n")
    sys.stderr.write(f"[manifest-v2] admit v2: {len(admit)} ADMIT rows\n")
    sys.stderr.write(f"[manifest-v2] phenos: {len(pheno)}\n")

    rmap_rows = []
    for _, row in region_map.iterrows():
        chrom = int(row["chrom"])
        rs = int(row["region_start"])
        re_ = int(row["region_end"])
        chrom_str = str(chrom)
        n_admit_v2 = int(((admit["primary_chrom"].astype(str) == chrom_str)
                          & (admit["primary_block_start"] >= rs)
                          & (admit["primary_block_end"] <= re_)
                          & (admit.get("status", pd.Series("ADMIT", index=admit.index)) == "ADMIT")).sum())
        rmap_rows.append({**row.to_dict(), "n_admit_v2": n_admit_v2})

    rmap_v2 = pd.DataFrame(rmap_rows)
    rmap_v2_keep = rmap_v2[rmap_v2["n_admit_v2"] >= 1].copy().reset_index(drop=True)
    atomic_write_tsv(rmap_v2, args.out_region_map)

    n_dropped = len(rmap_v2) - len(rmap_v2_keep)
    sys.stderr.write(f"[manifest-v2] regions surviving cv_R² >= 0.3: {len(rmap_v2_keep)} (dropped {n_dropped})\n")

    rows = []
    idx = 1
    for _, r in rmap_v2_keep.iterrows():
        for _, p in pheno.iterrows():
            rows.append({
                "idx": idx,
                "chrom": int(r["chrom"]),
                "region_id": r["region_id"],
                "pheno_code": p["pheno_code"],
                "study": p["study"],
                "n_admit_v2": int(r["n_admit_v2"]),
            })
            idx += 1
    manifest = pd.DataFrame(rows)
    atomic_write_tsv(manifest, args.out)
    sys.stderr.write(f"[manifest-v2] wrote {args.out}  ({len(manifest)} tasks)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
