#!/usr/bin/env python3
"""Emit the genotypes a CenSV-framed LD panel needs, selecting nothing."""
import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]
CENSV_LAB = _censv_os.environ["CENSV_LAB"]

ROOT = CENSV_ROOT + ""
SUF = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_"
       "qgap200_rgap200")
PIPE = os.path.join(ROOT, "analyses/ld_maps/code")
WEB = os.path.join(ROOT, "cenSV_web/web/public/data")
KG = CENSV_LAB + "/shared/1KG"
BED = (CENSV_COLLAB_ROOT + "/resources/"
       "chm13v2.0_censat_v2.1.active.bed")
GRCH38_BED = os.path.join(ROOT, "resources/grch38_centromeres.bed")
FLANK = 100_000
STEPS = (1, 2, 5, 10, 30, 100, 300)
SNP_MIN_MAF = 0.01

sys.path.insert(0, PIPE)
import plot_sv_ld_pubfig as P                                              # noqa: E402


def norm(x):
    s = str(x).strip()
    return "" if s in ("nan", "None", "") else s


def key_of(r):
    return (int(r.pos), int(r.end), r.svtype, norm(r.sv_category),
            norm(r.mapped_haplogroup))


def active_span(chrom):
    best = None
    for line in open(BED):
        f = line.split()
        if len(f) >= 3 and f[0] == chrom:
            s, e = int(f[1]), int(f[2])
            if best is None or e - s > best[1] - best[0]:
                best = (s, e)
    return best


def grch38_centromere(chrom):
    for line in open(GRCH38_BED):
        f = line.split()
        if len(f) >= 3 and f[0] == chrom:
            return int(f[1]), int(f[2])
    return None


def contig_length(vcf, chrom):
    for line in P.run_command(["bcftools", "view", "-h", vcf]).splitlines():
        if line.startswith("##contig=") and ("ID=%s," % chrom) in line:
            return int(line.split("length=")[1].split(">")[0].split(",")[0])
    raise ValueError("no contig length for %s in %s" % (chrom, vcf))


def edge_snp(vcf, chrom, start, end, innermost_is_last):
    if start > end:
        return None
    pick = "tail -1" if innermost_is_last else "head -1"
    cmd = ("bcftools query -r %s:%d-%d -f '%%POS\\n' %s | %s"
           % (chrom, start, end, vcf, pick))
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip()
    return int(out) if out else None


def flank_window(vcf, chrom, inner, side, chrlen):
    for k, mult in enumerate(STEPS):
        w = FLANK * mult
        if side == "L":
            lo, hi = max(1, inner - w), inner - 1
        else:
            lo, hi = inner + 1, min(chrlen, inner + w)
        p = edge_snp(vcf, chrom, lo, hi, innermost_is_last=(side == "L"))
        if p is None:
            continue
        if side == "L":
            start, end, gap = max(1, p - FLANK + 1), p, inner - 1 - p
        else:
            start, end, gap = p, min(chrlen, p + FLANK - 1), p - (inner + 1)
        note = None
        if k:
            note = ("the immediate flank held no 1000 Genomes call; this block is anchored "
                    "%.2f Mb out, on the nearest callable sequence" % (max(gap, 0) / 1e6))
        return start, end, max(gap, 0), note
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    chrom = a.chrom

    cat = pd.DataFrame(json.load(open(os.path.join(WEB, "catalog.json"))))
    lay = pd.DataFrame(json.load(open(os.path.join(WEB, "layers.json"))))
    hit = set(zip(lay["chrom"], lay["sv_index"]))
    sub = cat[cat.chrom == chrom]
    if not len(sub):
        sys.exit("%s is not in the catalog" % chrom)

    cdir = os.path.join(ROOT, "sv_extraction/decompose", SUF, chrom)
    g = np.load(os.path.join(cdir, "%s_sv_genotypes.npz" % chrom))["genotypes"]
    var = pd.read_csv(os.path.join(cdir, "%s_sv_variants.tsv" % chrom), sep="\t")
    samples = open(os.path.join(cdir, "%s_sv_samples.txt" % chrom)).read().split()
    variable = g.var(axis=1) > 0

    idx = {}
    for i, r in enumerate(var.itertuples()):
        idx.setdefault(key_of(r), []).append(i)

    rows, ids, marks, lost = [], [], [], []
    for r in sub.sort_values("pos").itertuples():
        hits = [i for i in idx.get(key_of(r), []) if variable[i]]
        row = None
        if len(hits) == 1:
            row = hits[0]
        elif len(hits) > 1:
            d = [abs(float(var.af.iloc[i]) - float(r.svfreq)) for i in hits]
            o = np.argsort(d)
            if d[o[1]] >= max(3 * d[o[0]], d[o[0]] + 0.01):
                row = hits[o[0]]
        if row is None:
            lost.append(r.sv_id)
            continue
        rows.append(row)
        ids.append(r.sv_id)
        marks.append(int((chrom, int(r.sv_index)) in hit))
    if lost:
        print("[%s] %d catalog rows unresolved: %s" % (chrom, len(lost), lost[:5]))

    sv_gt = g[rows, :]
    sv_var = var.iloc[rows].copy().reset_index(drop=True)
    sv_var["sv_id"] = ids
    sv_var["marked"] = marks

    span = active_span(chrom)
    cent = grch38_centromere(chrom)
    vcf = os.path.join(KG, "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_%s."
                            "recalibrated_variants.vcf.gz" % chrom)

    blocks, keep_samples, sides = [], samples, {}
    if cent is not None and os.path.exists(vcf):
        chrlen = contig_length(vcf, chrom)
        for side, inner in (("L", cent[0]), ("R", cent[1])):
            win = flank_window(vcf, chrom, inner, side, chrlen)
            if win is None:
                sides[side] = dict(region=None, n_snp=0, gap_bp=None,
                                   note="1000 Genomes holds no callable SNP within %d Mb of "
                                        "this edge" % (FLANK * STEPS[-1] // 1_000_000))
                continue
            start, end, gap, note = win
            region = "%s:%d-%d" % (chrom, start, end)
            try:
                gt, v, keep = P.load_snp_panel_from_vcf(
                    vcf_path=__import__("pathlib").Path(vcf), region=region,
                    target_samples=samples, min_maf=SNP_MIN_MAF)
            except ValueError as e:
                sides[side] = dict(region=region, n_snp=0, gap_bp=gap, note=str(e))
                continue
            v = v.copy()
            v["side"] = side
            blocks.append((side, gt, v))
            keep_samples = keep
            sides[side] = dict(region=region, n_snp=int(gt.shape[0]), gap_bp=int(gap),
                               note=note)
    elif cent is None:
        sides = {"L": dict(region=None, n_snp=0, gap_bp=None,
                           note="no GRCh38 centromere record"),
                 "R": dict(region=None, n_snp=0, gap_bp=None,
                           note="no GRCh38 centromere record")}

    if blocks:
        snp_gt = np.vstack([b[1] for b in blocks])
        snp_var = pd.concat([b[2] for b in blocks], ignore_index=True)
    else:
        snp_gt, snp_var = None, None

    pos = {s: i for i, s in enumerate(samples)}
    cols = [pos[s] for s in keep_samples]
    sv_gt = sv_gt[:, cols]

    if snp_gt is not None:
        impute = snp_gt.mean(axis=1)
        miss = ~np.isin(snp_gt, (0.0, 1.0, 2.0))
        snp_out = np.where(miss, 3, snp_gt).astype(np.int8)
        snp_var = snp_var.copy()
        snp_var["impute"] = impute
        snp_var["n_missing"] = miss.sum(axis=1).astype(int)
        n_missing_cells = int(miss.sum())
    else:
        snp_out = np.zeros((0, len(cols)), dtype=np.int8)
        n_missing_cells = 0

    out = os.path.join(a.out, chrom)
    os.makedirs(out, exist_ok=True)
    np.savez_compressed(os.path.join(out, "genotypes.npz"),
                        sv=sv_gt.astype(np.int8), snp=snp_out)
    SV_COLS = ["chrom", "pos", "end", "svtype", "svlen", "af", "maf",
               "sv_category", "mapped_haplogroup", "r2_score", "sv_id", "marked"]
    SNP_COLS = ["chrom", "pos", "end", "id", "svtype", "af", "maf", "side", "impute",
                "n_missing"]
    sv_var[[c for c in SV_COLS if c in sv_var.columns]].to_csv(
        os.path.join(out, "sv_variants.tsv"), sep="\t", index=False)
    if snp_var is not None:
        snp_var[[c for c in SNP_COLS if c in snp_var.columns]].to_csv(
            os.path.join(out, "snp_variants.tsv"), sep="\t", index=False)
    json.dump(dict(chrom=chrom, n_samples=len(cols), samples=keep_samples,
                   n_sv=int(sv_gt.shape[0]), n_catalog=int(len(sub)),
                   n_marked=int(sum(marks)), n_unresolved=len(lost),
                   n_snp=int(snp_out.shape[0]),
                   n_snp_missing_cells=n_missing_cells,
                   grch38_centromere=list(cent) if cent else None,
                   flanks=sides,
                   active_start=span[0] if span else None,
                   active_end=span[1] if span else None,
                   snp_min_maf=SNP_MIN_MAF, flank_bp=FLANK,
                   sv_assembly="CHM13v2.0", snp_assembly="GRCh38"),
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print("[%s] SVs %d of %d catalog (%d marked), SNPs %d (L %d, R %d), samples %d, "
          "%d SNP calls imputed"
          % (chrom, sv_gt.shape[0], len(sub), sum(marks), snp_out.shape[0],
             sides.get("L", {}).get("n_snp", 0), sides.get("R", {}).get("n_snp", 0),
             len(cols), n_missing_cells))


if __name__ == "__main__":
    main()
