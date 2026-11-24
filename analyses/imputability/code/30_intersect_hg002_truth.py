#!/usr/bin/env python3
"""30_intersect_hg002_truth.py Phase B — intersect HG002 production svim-asm calls vs Q100-truth svim-asm calls."""
from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJECT = Path(CENSV_ROOT + "")
SUFFIX = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"

PROD_VCFS = {
    "hap1": PROJECT / f"variant_calling/svim_asm/HG002_hap1_{SUFFIX}.vcf.gz",
    "hap2": PROJECT / f"variant_calling/svim_asm/HG002_hap2_{SUFFIX}.vcf.gz",
}
TRUTH_VCFS = {
    "PATERNAL": PROJECT / f"analyses/imputability/output/hg002_q100_truth/q100_PATERNAL_{SUFFIX}.vcf.gz",
    "MATERNAL": PROJECT / f"analyses/imputability/output/hg002_q100_truth/q100_MATERNAL_{SUFFIX}.vcf.gz",
}
ACTIVE_BED = PROJECT / "resources/chm13v2.0_censat_v2.1.active.bed"
OUT_TSV = PROJECT / "analyses/imputability/output/hg002_truth_intersect.tsv"

POS_TOL = 100
SVLEN_RATIO_MIN = 0.5


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.parent / f"{path.name}.tmp.{os.getpid()}"
    df.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)


def vcf_to_df(vcf_path: Path, label: str) -> pd.DataFrame:
    rows = []
    n_no_svtype = 0
    vcf = pysam.VariantFile(str(vcf_path))
    for rec in vcf.fetch():
        svtype = rec.info.get("SVTYPE")
        svlen = rec.info.get("SVLEN")
        end = rec.info.get("END", rec.stop)
        if svtype is None:
            n_no_svtype += 1
            continue
        if isinstance(svlen, (tuple, list)):
            svlen = svlen[0]
        if svlen is None:
            svlen_val = float("nan")
        else:
            svlen_val = abs(int(svlen))
        rows.append({
            "chrom": rec.chrom,
            "pos": int(rec.pos),
            "end": int(end) if end is not None else int(rec.pos),
            "svtype": str(svtype),
            "svlen": svlen_val,
            "label": label,
        })
    vcf.close()
    if n_no_svtype > 0:
        print(f"  [vcf_to_df] {vcf_path.name}: skipped {n_no_svtype} records without SVTYPE", flush=True)
    return pd.DataFrame(rows)


def load_active_hor_bed() -> pd.DataFrame:
    bed = pd.read_csv(ACTIVE_BED, sep="\t", header=None,
                       usecols=[0, 1, 2], names=["chrom", "start", "end"])
    return bed


def in_active_hor(chrom: str, pos: int, bed: pd.DataFrame) -> bool:
    sub = bed[bed["chrom"] == chrom]
    return bool(((sub["start"] <= pos) & (sub["end"] > pos)).any())


def match_one(prod_row: pd.Series, truth_df: pd.DataFrame) -> dict:
    chrom_match = truth_df["chrom"] == prod_row["chrom"]
    svtype_match = truth_df["svtype"] == prod_row["svtype"]
    sub = truth_df[chrom_match & svtype_match]
    if sub.empty:
        return {"truth_match": False, "best_truth_pos": -1, "best_truth_svlen": -1}

    pos_diff = (sub["pos"] - prod_row["pos"]).abs()
    sub = sub[pos_diff <= POS_TOL]
    if sub.empty:
        return {"truth_match": False, "best_truth_pos": -1, "best_truth_svlen": -1}

    prod_svlen = prod_row["svlen"]
    if np.isnan(prod_svlen) or prod_row["svtype"] == "BND":
        best = sub.iloc[(sub["pos"] - prod_row["pos"]).abs().argsort()[:1]].iloc[0]
        bp = int(best["pos"])
        bl = float(best["svlen"]) if not pd.isna(best["svlen"]) else -1
        return {"truth_match": True, "best_truth_pos": bp, "best_truth_svlen": bl}

    sub = sub.dropna(subset=["svlen"])
    if sub.empty:
        return {"truth_match": False, "best_truth_pos": -1, "best_truth_svlen": -1}

    ratio = np.minimum(sub["svlen"], prod_svlen) / np.maximum(sub["svlen"], prod_svlen)
    sub = sub[ratio >= SVLEN_RATIO_MIN]
    if sub.empty:
        return {"truth_match": False, "best_truth_pos": -1, "best_truth_svlen": -1}

    pos_diff = (sub["pos"] - prod_row["pos"]).abs()
    svlen_diff = (sub["svlen"] - prod_svlen).abs()
    sub = sub.assign(_score=pos_diff + svlen_diff).sort_values("_score")
    best = sub.iloc[0]
    return {"truth_match": True, "best_truth_pos": int(best["pos"]),
             "best_truth_svlen": int(best["svlen"])}


def main() -> int:
    print(f"[{pd.Timestamp.now()}] loading production VCFs...", flush=True)
    prod = pd.concat([vcf_to_df(p, f"prod_{h}") for h, p in PROD_VCFS.items()],
                      ignore_index=True)
    prod["prod_hap"] = prod["label"].str.replace("prod_", "")
    print(f"  production SVs: {len(prod)}", flush=True)

    print(f"[{pd.Timestamp.now()}] loading Q100 truth VCFs...", flush=True)
    for h, p in TRUTH_VCFS.items():
        if not p.exists():
            sys.exit(f"[ERR] truth VCF missing: {p}\n(Has Phase A3 completed?)")
    truth = pd.concat([vcf_to_df(p, f"truth_{h}") for h, p in TRUTH_VCFS.items()],
                       ignore_index=True)
    print(f"  truth SVs: {len(truth)}", flush=True)

    bed = load_active_hor_bed()
    print(f"  active-HOR BED rows: {len(bed)}", flush=True)

    print(f"[{pd.Timestamp.now()}] matching {len(prod)} production SVs...", flush=True)
    results = []
    for _, row in prod.iterrows():
        m = match_one(row, truth)
        results.append({
            "chrom": row["chrom"],
            "pos": row["pos"],
            "end": row["end"],
            "svtype": row["svtype"],
            "svlen": row["svlen"],
            "prod_hap": row["prod_hap"],
            "truth_match": m["truth_match"],
            "best_truth_pos": m["best_truth_pos"],
            "best_truth_svlen": m["best_truth_svlen"],
            "in_active_HOR": in_active_hor(row["chrom"], row["pos"], bed),
        })

    out = pd.DataFrame(results).sort_values(["chrom", "pos", "svtype", "svlen"]).reset_index(drop=True)
    atomic_write_tsv(out, OUT_TSV)
    print(f"[OK] wrote {len(out)} rows -> {OUT_TSV}", flush=True)

    print("\n=== Summary ===")
    print(f"  production total : {len(out)}")
    print(f"  TP (in truth)    : {int(out['truth_match'].sum())} ({100*out['truth_match'].mean():.1f}%)")
    print(f"  FP (not in truth): {int((~out['truth_match']).sum())} ({100*(~out['truth_match']).mean():.1f}%)")
    print(f"  in active HOR    : {int(out['in_active_HOR'].sum())}")
    print()
    print("  By svtype:")
    print(out.groupby("svtype")["truth_match"].agg(["count", "sum", "mean"]).to_string())
    print()
    print("  By chrom (active HOR only):")
    hor = out[out["in_active_HOR"]]
    print(hor.groupby("chrom")["truth_match"].agg(["count", "sum", "mean"]).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
