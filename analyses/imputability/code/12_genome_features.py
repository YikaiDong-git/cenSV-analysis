#!/usr/bin/env python3
"""Phase B2+B3 — locus-orthogonal genome features."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"
FA = ROOT / "resources" / "chm13v2.0.fa.gz"
CENSAT = ROOT / "resources" / "chm13v2.0_censat_v2.1.bed"
CENSAT_ACTIVE = ROOT / "resources" / "chm13v2.0_censat_v2.1.active.bed"

WINDOWS = [("gc_10kb", 10_000), ("gc_100kb", 100_000), ("gc_1Mb", 1_000_000)]

def family_class(name):
    head = name.split("_")[0]
    return head


def gc_in_window(fa, chrom, pos, half):
    start = max(0, pos - half)
    end = pos + half
    seq = fa.fetch(chrom, start, end).upper()
    if len(seq) == 0:
        return np.nan
    a = seq.count("A"); c = seq.count("C")
    g = seq.count("G"); t = seq.count("T")
    denom = a + c + g + t
    if denom == 0:
        return np.nan
    return (g + c) / denom


def b1_coord_check(fa, censat_active):
    feat = pd.read_csv(OUT / "per_sv_features_v3.tsv", sep="\t",
                       low_memory=False, usecols=["chrom", "pos", "sv_index"])
    chr6 = feat[feat["chrom"] == "chr6"].head(5)
    chr6_active = censat_active[censat_active["chrom"] == "chr6"]
    if len(chr6_active) == 0:
        print("[B1 FAIL] no chr6 active HOR in censat-active bed", flush=True)
        return False
    starts = chr6_active["start"].values
    ends = chr6_active["end"].values
    n_in = 0
    print("[B1] coord check: 5 chr6 SVs vs active HOR", flush=True)
    for _, row in chr6.iterrows():
        in_active = ((starts <= row["pos"]) & (row["pos"] < ends)).any()
        n_in += int(in_active)
        print(f"     sv_index={row['sv_index']} pos={row['pos']:>10d} "
              f"in_active_HOR={in_active}", flush=True)
    if n_in >= 4:
        print(f"[B1 PASS] {n_in}/5 chr6 SVs inside active HOR -> coords are CHM13v2",
              flush=True)
        return True
    print(f"[B1 WARN] only {n_in}/5 chr6 SVs in active HOR -> coords may be GRCh38",
          flush=True)
    return False


def main():
    feat_path = OUT / "per_sv_features_v3.tsv"
    feat = pd.read_csv(feat_path, sep="\t", low_memory=False)
    print(f"[INFO] feature table: {len(feat)} SVs", flush=True)

    cs = pd.read_csv(CENSAT, sep="\t", header=None,
                     names=["chrom", "start", "end", "name",
                            "score", "strand", "ts", "te", "rgb"])
    cs["family"] = cs["name"].apply(family_class)
    cs_active = pd.read_csv(CENSAT_ACTIVE, sep="\t", header=None,
                            names=["chrom", "start", "end", "name",
                                   "score", "strand", "ts", "te", "rgb"])
    print(f"[INFO] censat: {len(cs)} rows ({cs['family'].nunique()} families); "
          f"active sub-bed: {len(cs_active)} rows", flush=True)

    fa = pysam.FastaFile(str(FA))
    print(f"[INFO] fasta open: {len(fa.references)} contigs, "
          f"chr1 len={fa.get_reference_length('chr1')}", flush=True)
    if not b1_coord_check(fa, cs_active):
        print("[B1 WARN] coord check non-confirmatory; continuing with CHM13 ref",
              flush=True)

    rows = []
    for i, sv in feat[["sv_index", "chrom", "pos"]].iterrows():
        chrom = sv["chrom"]; pos = int(sv["pos"]); sv_idx = int(sv["sv_index"])
        gcs = {}
        for label, half in WINDOWS:
            gcs[label] = gc_in_window(fa, chrom, pos, half)

        cs_chr = cs[cs["chrom"] == chrom]
        hit = cs_chr[(cs_chr["start"] <= pos) & (pos < cs_chr["end"])]
        sat_family = hit["family"].iloc[0] if len(hit) > 0 else "none"

        cs_act_chr = cs_active[cs_active["chrom"] == chrom]
        in_active = bool(((cs_act_chr["start"] <= pos) &
                          (pos < cs_act_chr["end"])).any())

        local = cs_chr[(cs_chr["end"] > pos - 100_000) &
                       (cs_chr["start"] < pos + 100_000)]
        n_local = local["family"].nunique()

        rows.append({
            "sv_index": sv_idx, "chrom": chrom,
            **gcs,
            "sat_family": sat_family,
            "sat_active": int(in_active),
            "sat_local_complexity": int(n_local),
        })
        if (i + 1) % 1000 == 0:
            print(f"  ...processed {i+1}/{len(feat)}", flush=True)

    out = pd.DataFrame(rows)
    print(f"\n[INFO] per-SV genome features: {len(out)} rows", flush=True)
    print(out["sat_family"].value_counts().head(10).to_string(), flush=True)
    print(f"\nGC summary:\n{out[['gc_10kb', 'gc_100kb', 'gc_1Mb']].describe().to_string()}",
          flush=True)
    print(f"\nsat_active mean: {out['sat_active'].mean():.3f}", flush=True)
    print(f"sat_local_complexity mean: {out['sat_local_complexity'].mean():.2f}",
          flush=True)

    out_path = OUT / "per_sv_genome_features.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\n[OK] wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
