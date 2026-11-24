#!/usr/bin/env python3
from __future__ import annotations
import argparse
import os
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJECT = Path(CENSV_ROOT + "")
OUT_DIR = PROJECT / "analyses/imputability/output"
V3_TABLE = OUT_DIR / "per_sv_features_with_genome.tsv"
BAM_DIR = PROJECT / "alignments/svim_asm"
WINDOW_HALF = 1000

RNG_SEED = 42


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.parent / f"{path.name}.tmp.{os.getpid()}"
    df.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)


def find_bams() -> list[Path]:
    bams = sorted(BAM_DIR.glob("*_hap?_to_chm13.bam"))
    assert len(bams) > 0, f"No BAMs found under {BAM_DIR}"
    return bams


def fetch_locus_stats(bam: pysam.AlignmentFile, chrom: str, pos: int) -> dict[str, list]:
    mapqs, nms, secondaries = [], [], []
    start = max(0, pos - WINDOW_HALF)
    end = pos + WINDOW_HALF
    for aln in bam.fetch(chrom, start, end):
        if aln.is_unmapped:
            continue
        mapqs.append(aln.mapping_quality)
        try:
            nms.append(aln.get_tag("NM"))
        except KeyError:
            pass
        secondaries.append(bool(aln.is_secondary))
    return {"mapq": mapqs, "nm": nms, "sec": secondaries}


def process_chrom(chrom: str) -> pd.DataFrame:
    df = pd.read_csv(V3_TABLE, sep="\t")
    sub = df[df["chrom"] == chrom].copy()
    if sub.empty:
        print(f"[WARN] no SVs for {chrom} in v3 table", flush=True)
        return pd.DataFrame()

    sub = sub.sort_values(["chrom", "pos"]).reset_index(drop=True)
    print(f"[{chrom}] {len(sub)} SVs, opening {len(find_bams())} BAMs...", flush=True)

    bams = [pysam.AlignmentFile(str(p)) for p in find_bams()]
    print(f"[{chrom}] BAMs opened.", flush=True)

    rows = []
    for _, sv in sub.iterrows():
        all_mapqs, all_nms, all_secs = [], [], []
        n_alns = 0
        for bam in bams:
            stats = fetch_locus_stats(bam, chrom, int(sv["pos"]))
            all_mapqs.extend(stats["mapq"])
            all_nms.extend(stats["nm"])
            all_secs.extend(stats["sec"])
            n_alns += len(stats["mapq"])
        if n_alns == 0:
            rows.append(
                {
                    "chrom": chrom,
                    "sv_index": int(sv["sv_index"]) if "sv_index" in sv else -1,
                    "pos": int(sv["pos"]),
                    "mean_MAPQ_at_locus": np.nan,
                    "median_MAPQ_at_locus": np.nan,
                    "n_contigs_at_locus": 0,
                    "mean_NM_at_locus": np.nan,
                    "frac_secondary_at_locus": np.nan,
                }
            )
        else:
            mapq_arr = np.array(all_mapqs)
            sec_arr = np.array(all_secs, dtype=bool)
            nm_arr = np.array(all_nms) if all_nms else np.array([np.nan])
            rows.append(
                {
                    "chrom": chrom,
                    "sv_index": int(sv["sv_index"]) if "sv_index" in sv else -1,
                    "pos": int(sv["pos"]),
                    "mean_MAPQ_at_locus": float(mapq_arr.mean()),
                    "median_MAPQ_at_locus": float(np.median(mapq_arr)),
                    "n_contigs_at_locus": n_alns,
                    "mean_NM_at_locus": float(np.nanmean(nm_arr)),
                    "frac_secondary_at_locus": float(sec_arr.mean()),
                }
            )

    for bam in bams:
        bam.close()

    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True, help="e.g. chr1")
    ap.add_argument("--out-dir", default=str(OUT_DIR / "mapq_per_chrom"))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"mapq_features_{args.chrom}.tsv"

    if out_path.exists():
        print(f"[{args.chrom}] {out_path} exists, skip", flush=True)
        return 0

    df = process_chrom(args.chrom)
    if df.empty:
        return 0

    df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)
    atomic_write_tsv(df, out_path)
    print(f"[{args.chrom}] wrote {len(df)} rows -> {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
