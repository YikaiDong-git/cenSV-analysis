#!/usr/bin/env python3
"""Build the manifest for the Phase B z=0 pad vs drop head-to-head."""
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
    p.add_argument("--region-map-v2", type=Path,
                   default=BASE / "output/full_restart/manifests/region_map_v2_with_admit.tsv")
    p.add_argument("--pheno-list", type=Path, default=BASE / "inputs/pheno_top30.tsv")
    p.add_argument("--out", type=Path,
                   default=BASE / "output/full_restart/manifests/zpad_headhead_manifest.tsv")
    p.add_argument("--chroms", nargs="+", type=int, default=[1, 19])
    return p.parse_args()


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.parent / f"{path.stem}.tmp.{os.getpid()}{path.suffix}"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def main() -> int:
    args = parse_args()
    rmap = pd.read_csv(args.region_map_v2, sep="\t")
    pheno = pd.read_csv(args.pheno_list, sep="\t")
    rmap = rmap[(rmap["chrom"].isin(args.chroms)) & (rmap["n_admit_v2"] >= 1)].reset_index(drop=True)

    rows = []
    idx = 1
    for _, r in rmap.iterrows():
        for _, p in pheno.iterrows():
            for variant in ("pad", "drop"):
                rows.append({
                    "idx": idx,
                    "chrom": int(r["chrom"]),
                    "region_id": r["region_id"],
                    "pheno_code": p["pheno_code"],
                    "study": p["study"],
                    "variant": variant,
                })
                idx += 1

    manifest = pd.DataFrame(rows)
    atomic_write_tsv(manifest, args.out)
    sys.stderr.write(f"[zpad-manifest] wrote {args.out}  ({len(manifest)} tasks)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
