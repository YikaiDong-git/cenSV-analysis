#!/usr/bin/env python3
import os, subprocess, glob, csv, math
import numpy as np

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
SC = "/tmp/phewide"
INP = f"{SC}/inputs"
PGDIR = f"{BASE}/aou_replication/results_vds/step28_phewide/pgen"
PGENS = sorted(glob.glob(f"{PGDIR}/svpgen300_*.pgen"),
               key=lambda p: int(p.split("svpgen300_")[1].split(".pgen")[0]))
OUTA = f"{BASE}/aou_replication/results_vds/step30_phewide"; os.makedirs(OUTA, exist_ok=True)
PLINK2 = f"{HOME}/.local/bin/plink2"
CHUNK = 40
TMP = f"{SC}/glm_tmp"; os.makedirs(TMP, exist_ok=True)


def tiers():
    return {r["sv_locus"]: r["tier"]
            for r in csv.DictReader(open(f"{SC}/phewide_pgen/sv_tier.tsv"), delimiter="\t")}


def pheno_cols(pf):
    with open(pf) as f:
        return f.readline().rstrip("\n").split("\t")[2:]


def write_chunk(pf, cols_all, sel, path):
    keep = [0, 1] + [2 + cols_all.index(c) for c in sel]
    with open(pf) as f, open(path, "w") as fo:
        for line in f:
            p = line.rstrip("\n").split("\t")
            fo.write("\t".join(p[i] for i in keep) + "\n")


def run_arm(pf, kind):
    cols = pheno_cols(pf)
    print(f"[{kind}] {len(cols)} phenotypes -> {((len(cols)-1)//CHUNK)+1} chunks", flush=True)
    glm = (["--glm", "hide-covar", "firth-fallback", "cols=+a1freq,+beta"] if kind == "disease"
           else ["--glm", "hide-covar", "cols=+a1freq,+beta"])
    ext = "glm.logistic.hybrid" if kind == "disease" else "glm.linear"
    rows = []
    nchunk = ((len(cols) - 1) // CHUNK) + 1
    for ci in range(0, len(cols), CHUNK):
        sel = cols[ci:ci + CHUNK]
        cpath = f"{TMP}/chunk.pheno"; write_chunk(pf, cols, sel, cpath)
        for pgen in PGENS:
            pfx = pgen[:-len(".pgen")]
            args = [PLINK2, "--pfile", pfx, *glm, "--pheno", cpath,
                    "--covar", f"{INP}/covar.txt", "--covar-variance-standardize",
                    "--vif", "100000", "--max-corr", "0.9999", "--out", f"{TMP}/g", "--memory", "640"]
            if kind == "disease":
                args.append("--1")
            r = subprocess.run(args, capture_output=True, text=True)
            if r.returncode != 0:
                print(f"  chunk {ci//CHUNK+1} pgen {os.path.basename(pfx)} rc={r.returncode}: {r.stderr[-200:]}", flush=True)
            for fpath in glob.glob(f"{TMP}/g.*.{ext}"):
                ph = os.path.basename(fpath)[len("g."):-len("." + ext)]
                for x in csv.DictReader(open(fpath), delimiter="\t"):
                    if x.get("TEST") != "ADD":
                        continue
                    p = x["P"] if x["P"] not in ("", "NA") else "nan"
                    rows.append([x["ID"], ph, int(x["OBS_CT"]), x.get("A1_FREQ", ""),
                                 x.get("BETA", ""), x.get("SE", ""), float(p)])
            for junk in glob.glob(f"{TMP}/g.*"):
                os.remove(junk)
        print(f"  pheno-chunk {ci//CHUNK+1}/{nchunk}: cumulative {len(rows)} tests", flush=True)
    return rows


def tiered_bh(rows, tiermap):
    qcol = {}
    for tier in ("signal", "null_ext"):
        idx = [i for i, r in enumerate(rows) if tiermap.get(r[0]) == tier and math.isfinite(r[6])]
        if not idx:
            continue
        ps = np.array([rows[i][6] for i in idx]); n = len(ps)
        order = np.argsort(ps); prev = 1.0
        for rank, oi in enumerate(order[::-1]):
            prev = min(prev, ps[oi] * n / (n - rank)); qcol[idx[oi]] = prev
    return qcol


def write_arm(rows, tiermap, path):
    qcol = tiered_bh(rows, tiermap)
    with open(path, "w") as f:
        f.write("sv_locus\ttier\tpheno\tn\ta1_freq\tbeta\tse\tp\tq_bh_tier\n")
        for i, r in enumerate(rows):
            q = qcol.get(i, float("nan"))
            f.write(f"{r[0]}\t{tiermap.get(r[0],'')}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t"
                    f"{'' if not math.isfinite(r[6]) else f'{r[6]:.4g}'}\t"
                    f"{'' if not math.isfinite(q) else f'{q:.4g}'}\n")
    sig = sum(1 for i in range(len(rows)) if math.isfinite(qcol.get(i, float('nan'))) and qcol[i] < 0.05)
    sigs = sum(1 for i, r in enumerate(rows) if tiermap.get(r[0]) == "signal"
               and math.isfinite(qcol.get(i, float('nan'))) and qcol[i] < 0.05)
    print(f"wrote {path}: {len(rows)} tests, q<0.05 {sig} (signal-tier {sigs})", flush=True)


def main():
    tiermap = tiers()
    for pf, kind, out in [(f"{INP}/quant.pheno", "quant", f"{OUTA}/phewide_quant.tsv"),
                          (f"{INP}/disease.pheno", "disease", f"{OUTA}/phewide_disease.tsv")]:
        if not os.path.exists(pf):
            print(f"SKIP {kind}: {pf} missing", flush=True); continue
        rows = run_arm(pf, kind)
        write_arm(rows, tiermap, out)
    print("GLM_DONE", flush=True)


if __name__ == "__main__":
    main()
