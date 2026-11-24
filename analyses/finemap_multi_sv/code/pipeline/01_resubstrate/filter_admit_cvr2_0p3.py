#!/usr/bin/env python3
"""Re-substrate: filter ADMIT TSV to cv_R² >= 0.3 (was >= 0.1 admit floor)."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]


SRC_ADMIT = Path(
    CENSV_ROOT + "/"
    "analyses/sv_catalog/output/ref_ukbb/admit/admit_decompose_cv01.tsv"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--src", type=Path, default=SRC_ADMIT)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--cv-r2-min", type=float, default=0.30)
    return p.parse_args()


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.parent / f"{path.stem}.tmp.{os.getpid()}{path.suffix}"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def main() -> int:
    args = parse_args()

    src = pd.read_csv(args.src, sep="\t")
    n_total = len(src)
    admit = src[src["status"] == "ADMIT"].copy()
    n_admit = len(admit)
    filt = admit[admit["cv_r2"] >= args.cv_r2_min].copy().reset_index(drop=True)
    n_keep = len(filt)

    sys.stderr.write(f"[resubstrate] source: {n_total} rows ({n_admit} ADMIT)\n")
    sys.stderr.write(f"[resubstrate] cv_r2 >= {args.cv_r2_min}: kept {n_keep} / {n_admit} ({100*n_keep/n_admit:.1f}%)\n")

    per_chrom = filt.groupby("primary_chrom").size().reset_index(name="n_admit_v2")
    sys.stderr.write(f"[resubstrate] per chrom counts:\n{per_chrom.to_string(index=False)}\n")

    atomic_write_tsv(filt, args.out)
    sys.stderr.write(f"[resubstrate] wrote {args.out}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
