#!/usr/bin/env python3
"""Build a sample manifest for 1KG ONT Modkit methylation BED files."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methylation-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--file-suffix",
        default="_combined.split5hmc.bed.gz",
        help="One suffix or a comma-separated set of required files (for example hp1,hp2 calls).",
    )
    parser.add_argument(
        "--deduplicate-samples",
        choices=["largest_file", "first"],
        default="largest_file",
        help="How to choose when multiple directories resolve to the same sample ID.",
    )
    return parser.parse_args()


def parse_sample_dir(dirname: str) -> dict[str, str]:
    if "-ONT-" in dirname:
        sample_id, rest = dirname.split("-ONT-", 1)
    else:
        sample_id, rest = dirname.split("-", 1)[0], dirname

    tokens = re.split(r"[-_]+", rest)
    reference = next((tok for tok in tokens if tok.lower() in {"hg38", "chm13"}), "")
    pore = next((tok for tok in tokens if re.fullmatch(r"R[0-9]+", tok)), "")
    kit = next((tok for tok in tokens if tok.startswith("LSK")), "")
    caller = next((tok for tok in tokens if tok.startswith("dorado") or tok.startswith("guppy")), "")
    model_tokens = [tok for tok in tokens if tok and tok not in {reference, pore, kit, caller}]
    model = "_".join(model_tokens)
    return {
        "sample_id": sample_id,
        "reference": reference,
        "pore": pore,
        "kit": kit,
        "caller": caller,
        "model": model,
    }


def build_manifest(root: Path, file_suffix: str, deduplicate: str) -> pd.DataFrame:
    file_suffixes = [suffix.strip() for suffix in file_suffix.split(",") if suffix.strip()]
    if not file_suffixes:
        raise ValueError("At least one methylation file suffix is required")
    rows = []
    for sample_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        beds = []
        for suffix in file_suffixes:
            matches = sorted(sample_dir.glob(f"*{suffix}"))
            if not matches:
                beds = []
                break
            beds.append(matches[0])
        if not beds:
            continue
        parsed = parse_sample_dir(sample_dir.name)
        rows.append(
            {
                **parsed,
                "sample_dir": str(sample_dir),
                "directory_name": sample_dir.name,
                "methylation_bed": str(beds[0]),
                "methylation_beds": ";".join(map(str, beds)),
                "n_methylation_beds": len(beds),
                "file_size_bytes": sum(bed.stat().st_size for bed in beds),
                "has_hp1": bool(list(sample_dir.glob("*_hp1.bed.gz"))),
                "has_hp2": bool(list(sample_dir.glob("*_hp2.bed.gz"))),
            }
        )

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    df["duplicate_rank"] = 0
    if df["sample_id"].duplicated().any():
        if deduplicate == "largest_file":
            df = df.sort_values(["sample_id", "file_size_bytes", "directory_name"], ascending=[True, False, True])
        else:
            df = df.sort_values(["sample_id", "directory_name"])
        df["duplicate_rank"] = df.groupby("sample_id").cumcount()
        df = df[df["duplicate_rank"] == 0].copy()
    return df.sort_values("sample_id").reset_index(drop=True)


def main() -> int:
    args = parse_args()
    if not args.methylation_root.exists():
        raise FileNotFoundError(f"Methylation root not found: {args.methylation_root}")
    manifest = build_manifest(args.methylation_root, args.file_suffix, args.deduplicate_samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
