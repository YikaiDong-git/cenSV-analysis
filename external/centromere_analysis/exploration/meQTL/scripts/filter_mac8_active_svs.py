#!/usr/bin/env python3
"""Restrict completed trans physical-family metadata to the 68 MAC8-active SVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--physical-families", type=Path, required=True)
    parser.add_argument("--active-families", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    physical = pd.read_csv(args.physical_families, sep="\t")
    active = pd.read_csv(args.active_families, sep="\t")
    keep = set(active["physical_family_id"].astype(str))
    result = physical[physical["physical_family_id"].astype(str).isin(keep)].copy()
    if len(result) != 68 or result["physical_family_id"].nunique() != 68:
        raise ValueError(f"Expected exactly 68 active physical families; found {len(result)}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
