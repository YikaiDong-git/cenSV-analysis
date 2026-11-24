#!/usr/bin/env python3
"""How many 1000 Genomes SNPs sit in a wider window, and what they would cost to ship."""
import argparse
import json
import os
import subprocess
import sys

import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_LAB = _censv_os.environ["CENSV_LAB"]

ROOT = CENSV_ROOT + ""
PIPE = os.path.join(ROOT, "analyses/ld_maps/code")
FRAMES = os.path.join(ROOT, "analyses/ld_maps/output/_censv_frames")
KG = CENSV_LAB + "/shared/1KG"

TIERS = (100_000, 250_000, 1_000_000, 2_500_000, 5_000_000, 10_000_000)
MAX_TIER = max(TIERS)
MIN_MAF = 0.01

sys.path.insert(0, PIPE)
import plot_sv_ld_pubfig as P                                              # noqa: E402


def contig_length(vcf, chrom):
    for line in P.run_command(["bcftools", "view", "-h", vcf]).splitlines():
        if line.startswith("##contig=") and ("ID=%s," % chrom) in line:
            return int(line.split("length=")[1].split(">")[0].split(",")[0])
    raise ValueError("no contig length for %s in %s" % (chrom, vcf))


def stream_maf(vcf, region, sample_file):
    cmd = ["bcftools", "query", "-r", region, "-S", sample_file,
           "-f", "%POS\t%REF\t%ALT[\t%GT]\n", vcf]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    out = []
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        pos, ref, alt = f[0], f[1], f[2]
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        d = np.array([P.gt_to_dosage(g) for g in f[3:]], dtype=float)
        if np.all(np.isnan(d)):
            continue
        af = float(np.nanmean(d) / 2.0)
        maf = min(af, 1.0 - af)
        if maf < MIN_MAF:
            continue
        obs = d[~np.isnan(d)]
        if obs.min() == obs.max():
            continue
        out.append((int(pos), maf, int(np.isnan(d).sum())))
    proc.stdout.close()
    if proc.wait() != 0:
        raise RuntimeError("bcftools failed on %s" % region)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    chrom = a.chrom

    summary = json.load(open(os.path.join(FRAMES, chrom, "summary.json")))
    samples = summary["samples"]
    vcf = os.path.join(KG, "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_%s."
                            "recalibrated_variants.vcf.gz" % chrom)
    chrlen = contig_length(vcf, chrom)

    sample_file = os.path.join(a.out, "%s.samples.txt" % chrom)
    with open(sample_file, "w") as fh:
        fh.write("\n".join(samples) + "\n")

    rec = {"chrom": chrom, "n_samples": len(samples), "min_maf": MIN_MAF,
           "tiers_bp": list(TIERS), "flanks": {}}

    for side, fl in summary["flanks"].items():
        lo0, hi0 = (int(x) for x in fl["region"].split(":")[1].split("-"))
        inner = hi0 if side == "L" else lo0
        if side == "L":
            lo, hi = max(1, inner - MAX_TIER + 1), inner
        else:
            lo, hi = inner, min(chrlen, inner + MAX_TIER - 1)
        rows = stream_maf(vcf, "%s:%d-%d" % (chrom, lo, hi), sample_file)
        pos = np.array([r[0] for r in rows], dtype=np.int64)
        miss = np.array([r[2] for r in rows], dtype=np.int64)
        dist = (inner - pos) if side == "L" else (pos - inner)

        per_tier = {}
        for t in TIERS:
            sel = dist < t
            per_tier[str(t)] = {
                "n_snp": int(sel.sum()),
                "n_missing_cells": int(miss[sel].sum()),
                "span_bp": int(dist[sel].max() + 1) if sel.any() else 0,
            }
        rec["flanks"][side] = {
            "shipped_region": fl["region"], "inner_anchor": int(inner),
            "search_region": "%s:%d-%d" % (chrom, lo, hi),
            "n_snp_total": len(pos), "per_tier": per_tier,
        }
        print("[%s %s] %s" % (chrom, side,
                              " ".join("%s=%d" % (t, per_tier[t]["n_snp"])
                                       for t in per_tier)), flush=True)

    path = os.path.join(a.out, "%s_tiers.json" % chrom)
    tmp = path + ".tmp"
    json.dump(rec, open(tmp, "w"), indent=1)
    os.replace(tmp, path)
    os.remove(sample_file)
    print("wrote %s" % path)


if __name__ == "__main__":
    main()
