#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]


BASE_BCAC = Path(
    CENSV_ROOT + "/"
    "analyses/finemap_multi_sv/inputs/bcac_2017"
)
CHAIN = Path(
    CENSV_WORKSPACE + "/_shared/tools/"
    "process_neale_lab_GWAS/v_master_83578d1_20260514/hg19ToHg38.over.chain"
)
SHARED = Path(CENSV_WORKSPACE + "/_shared")
BGZIP = SHARED / "envs/workspace/v1/.pixi/envs/biotools/bin/bgzip"
TABIX = SHARED / "envs/workspace/v1/.pixi/envs/biotools/bin/tabix"
LIFTOVER = (
    SHARED / "tools/process_neale_lab_GWAS/v_master_83578d1_20260514/"
    ".pixi/envs/default/bin/liftOver"
)

COL_VARNAME    = 0
COL_CHR        = 2
COL_POS_B37    = 3
COL_A0         = 4
COL_A1         = 5
COL_EAF_CTRL   = 46
COL_BETA       = 47
COL_SE         = 48
COL_PVAL       = 49

N_ERNEG_CASES = 21468
N_ERNEG_CTRL  = 100594
N_TOTAL       = N_ERNEG_CASES + N_ERNEG_CTRL


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bcac-gz", type=Path,
                   default=BASE_BCAC / "oncoarray_bcac_public_release_oct17.txt.gz")
    p.add_argument("--chain", type=Path, default=CHAIN)
    p.add_argument("--out-dir", type=Path, default=BASE_BCAC / "processed_b38")
    p.add_argument("--chrom-int", type=int, default=1)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    bed_b37 = args.out_dir / f"bcac_erneg_chr{args.chrom_int}_b37.bed"
    bed_b38 = args.out_dir / f"bcac_erneg_chr{args.chrom_int}_b38.bed"
    bed_unmap = args.out_dir / f"bcac_erneg_chr{args.chrom_int}_unmap.bed"
    tmp_subset = args.out_dir / f"bcac_erneg_chr{args.chrom_int}_subset.tsv"

    sys.stderr.write(f"[bcac] reading {args.bcac_gz} → filtering chr={args.chrom_int}\n")

    n_total = 0
    n_chrom = 0
    n_valid_erneg = 0

    with gzip.open(args.bcac_gz, "rt") as fin, \
         open(tmp_subset, "w") as fout, \
         open(bed_b37, "w") as fbed:
        header = fin.readline().rstrip("\n").split("\t")
        if header[COL_BETA] != "bcac_onco_icogs_gwas_erneg_beta":
            raise ValueError(
                f"COL_BETA index {COL_BETA} does not match expected col name; got '{header[COL_BETA]}'"
            )

        for ln in fin:
            n_total += 1
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 50:
                continue
            ch = parts[COL_CHR]
            if ch != str(args.chrom_int):
                continue
            n_chrom += 1
            try:
                pos = int(parts[COL_POS_B37])
            except ValueError:
                continue
            beta = parts[COL_BETA]
            se = parts[COL_SE]
            pval = parts[COL_PVAL]
            if beta == "NULL" or se == "NULL" or pval == "NULL":
                continue
            try:
                _b = float(beta); _s = float(se); _p = float(pval)
            except ValueError:
                continue
            if _s <= 0:
                continue
            n_valid_erneg += 1
            vname = parts[COL_VARNAME]
            a0 = parts[COL_A0]
            a1 = parts[COL_A1]
            eaf = parts[COL_EAF_CTRL]
            fout.write(f"{vname}\t{ch}\t{pos}\t{a0}\t{a1}\t{eaf}\t{beta}\t{se}\t{pval}\n")
            fbed.write(f"chr{ch}\t{pos-1}\t{pos}\t{vname}\n")

    sys.stderr.write(f"[bcac] total rows scanned: {n_total}\n")
    sys.stderr.write(f"[bcac] chr{args.chrom_int} rows: {n_chrom}\n")
    sys.stderr.write(f"[bcac] valid ER-neg rows: {n_valid_erneg}\n")

    sys.stderr.write(f"[bcac] running UCSC liftOver: {LIFTOVER}\n")
    if not LIFTOVER.exists():
        raise FileNotFoundError(f"liftOver binary not at {LIFTOVER}")
    subprocess.run(
        [str(LIFTOVER), str(bed_b37), str(args.chain), str(bed_b38), str(bed_unmap)],
        check=True,
    )

    sys.stderr.write(f"[bcac] lifted: {bed_b38} ; unmapped: {bed_unmap}\n")

    var_to_b38: dict[str, tuple[str, int]] = {}
    with open(bed_b38) as fb:
        for ln in fb:
            ch38, s38, e38, vname = ln.rstrip("\n").split("\t")
            if ch38 != f"chr{args.chrom_int}":
                continue
            var_to_b38[vname] = (ch38.replace("chr", ""), int(e38))

    sys.stderr.write(f"[bcac] lifted to b38 (same chrom): {len(var_to_b38)}\n")

    out_tsv = args.out_dir / f"bcac_erneg_chr{args.chrom_int}_b38.tsv"
    rows_b38 = []
    with open(tmp_subset) as fin:
        for ln in fin:
            vname, ch, b37pos, a0, a1, eaf, beta, se, pval = ln.rstrip("\n").split("\t")
            if vname not in var_to_b38:
                continue
            chrom_int, pos_b38 = var_to_b38[vname]
            snp_id = f"{chrom_int}_{pos_b38}_{a0}_{a1}"
            try:
                eaf_v = float(eaf) if eaf not in ("NULL", "") else None
            except ValueError:
                eaf_v = None
            maf = min(eaf_v, 1 - eaf_v) if eaf_v is not None else "NA"
            rows_b38.append(
                f"{chrom_int}\t{pos_b38}\t{snp_id}\t{maf}\t{N_TOTAL}\t{beta}\t{se}\t{pval}\n"
            )

    sys.stderr.write(f"[bcac] writing {out_tsv} ({len(rows_b38)} rows)\n")
    rows_b38.sort(key=lambda r: int(r.split("\t")[1]))
    with open(out_tsv, "w") as fo:
        fo.write("CHROM_INT\tPOS\tSNP\tMAF\tN\tBETA\tSE\tPVALUE\n")
        fo.writelines(rows_b38)

    sys.stderr.write(f"[bcac] bgzipping + tabix-indexing\n")
    subprocess.run([str(BGZIP), "-f", str(out_tsv)], check=True)
    out_bgz = out_tsv.with_suffix(".tsv.gz")
    subprocess.run(
        [str(TABIX), "-f", "-s", "1", "-b", "2", "-e", "2", "-S", "1", str(out_bgz)],
        check=True,
    )
    sys.stderr.write(f"[bcac] DONE: {out_bgz}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
