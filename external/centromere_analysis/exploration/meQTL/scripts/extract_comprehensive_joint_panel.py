#!/usr/bin/env python3
"""Range-extract one comprehensive joint-CHM13 chromosome panel."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from comprehensive_conditioning import sha256_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--source-url", required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--index-output", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--bcftools", default="bcftools")
    return parser.parse_args()


def run(command: list[str]) -> str:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(command)}\n{result.stderr[-4000:]}")
    return result.stdout


def main() -> None:
    args = parse_args()
    expected_samples = [value for value in args.samples.read_text().splitlines() if value]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metrics.parent.mkdir(parents=True, exist_ok=True)
    scratch_root = os.environ.get("TMPDIR")
    with tempfile.TemporaryDirectory(prefix=f"joint_panel_{args.chrom}_", dir=scratch_root) as temp_dir:
        local_vcf = Path(temp_dir) / f"{args.chrom}.vcf.gz"
        run([
            args.bcftools, "view", "--regions-file", str(args.regions),
            "--samples-file", str(args.samples), "--types", "snps",
            "--min-alleles", "2", "--max-alleles", "2", "--output-type", "z",
            "--output", str(local_vcf), args.source_url,
        ])
        run([args.bcftools, "index", "--force", "--tbi", str(local_vcf)])
        observed_samples = run([args.bcftools, "query", "--list-samples", str(local_vcf)]).splitlines()
        if observed_samples != expected_samples:
            raise ValueError(f"{args.chrom}: extracted sample order does not match canonical 124-sample order")
        n_records = int(run([args.bcftools, "index", "--nrecords", str(local_vcf)]).strip())
        local_index = Path(str(local_vcf) + ".tbi")
        shutil.move(local_vcf, args.output)
        shutil.move(local_index, args.index_output)
    metrics = {
        "status": "ok", "chrom": args.chrom, "source_url": args.source_url,
        "source_index_url": args.source_url + ".tbi", "range_extraction": True,
        "n_samples": len(expected_samples), "n_biallelic_snps": n_records,
        "regions": str(args.regions), "regions_sha256": sha256_file(args.regions),
        "vcf": str(args.output), "vcf_bytes": args.output.stat().st_size,
        "vcf_sha256": sha256_file(args.output), "index": str(args.index_output),
        "index_bytes": args.index_output.stat().st_size, "index_sha256": sha256_file(args.index_output),
    }
    args.metrics.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
