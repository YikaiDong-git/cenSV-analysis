#!/usr/bin/env python3
"""Step 60 (ext2 arm) — build WITHIN-AMR principal components, reproducibly."""
import os, sys, csv, time, shutil, subprocess, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S
import step56_ancestry_disease_ext2 as P

HOME = os.path.expanduser("~")
PLINK_SRC = f"{HOME}/workspace/working/resources/plink2"
PLINK = "/tmp/plink2_step60"
ACAF = f"{HOME}/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen"
CHROM = "chr22"
OUT = f"{S.OUTA}/pca"
SCRATCH = "/tmp/step60_pca"


def ensure_plink():
    if not os.path.exists(PLINK) or os.path.getsize(PLINK) != os.path.getsize(PLINK_SRC):
        shutil.copy2(PLINK_SRC, PLINK)
    os.chmod(PLINK, 0o755)
    v = subprocess.run([PLINK, "--version"], capture_output=True, text=True).stdout.strip()
    print(f"plink2: {v}", flush=True)


def psam_has_fid(pfile):
    with open(f"{pfile}.psam") as f:
        for line in f:
            if line.startswith("#"):
                cols = line.lstrip("#").rstrip("\n").split("\t")
                return "FID" in [c.upper() for c in cols]
            break
    return False


def amr_keep_file(path, with_fid):
    with open(f"{S.INP}/covar.txt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        fid_i = hdr.index("FID") if "FID" in hdr else None
        iid_i = hdr.index("IID")
        dcol = {d: hdr.index(d) for d in P.ANC_DUMMIES}
        n = 0
        with open(path, "w") as o:
            if with_fid:
                o.write("#FID\tIID\n")
            else:
                o.write("#IID\n")
            for line in f:
                p = line.rstrip("\n").split("\t")
                if p[dcol["anc_amr"]] not in ("1", "1.0"):
                    continue
                if with_fid:
                    fid = p[fid_i] if fid_i is not None else p[iid_i]
                    o.write(f"{fid}\t{p[iid_i]}\n")
                else:
                    o.write(f"{p[iid_i]}\n")
                n += 1
    return n


def run(cmd, tag):
    t0 = time.time()
    print(f"\n$ {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-4000:] + "\n" + r.stderr[-4000:] + "\n")
        raise SystemExit(f"ABORT: {tag} failed (rc={r.returncode})")
    print(f"  {tag} ok ({time.time()-t0:.0f}s)", flush=True)
    return r.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kpc", type=int, default=20)
    ap.add_argument("--maf", type=float, default=0.01)
    ap.add_argument("--window", type=int, default=200, help="--indep-pairwise window (variants)")
    ap.add_argument("--step", type=int, default=50, help="--indep-pairwise step")
    ap.add_argument("--r2", type=float, default=0.2, help="--indep-pairwise r^2 ceiling")
    ap.add_argument("--min-snps", type=int, default=10000,
                    help="hard floor on pruned SNPs; below this the PCs are too noisy and would "
                         "under-adjust, biasing the rescue toward SURVIVES")
    ap.add_argument("--memory", type=int, default=12000, help="plink2 --memory MB")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True); os.makedirs(SCRATCH, exist_ok=True)
    ensure_plink()

    pfile = f"{ACAF}/acaf_threshold.{CHROM}"
    for ext in (".pgen", ".pvar", ".psam"):
        if not os.path.exists(pfile + ext):
            raise SystemExit(f"ABORT: missing {pfile+ext}")

    with_fid = psam_has_fid(pfile)
    keep = f"{SCRATCH}/amr_keep.txt"
    n = amr_keep_file(keep, with_fid)
    print(f"AMR samples from covar.txt: {n}  (psam {'HAS' if with_fid else 'has NO'} FID column -> "
          f"keep file written with {'FID+IID' if with_fid else 'IID only'})", flush=True)
    if n < 1000:
        raise SystemExit(f"ABORT: only {n} AMR samples parsed — check covar.txt dummy columns")
    psam_ids = set()
    with open(f"{pfile}.psam") as f:
        for line in f:
            if line.startswith("#"):
                continue
            psam_ids.add(line.split("\t")[1 if with_fid else 0].strip())
    keep_ids = {l.split("\t")[-1].strip() for l in open(keep) if not l.startswith("#")}
    inter = len(keep_ids & psam_ids)
    print(f"  keep/psam IID intersection: {inter}/{len(keep_ids)}", flush=True)
    if inter < 0.9 * len(keep_ids):
        raise SystemExit(f"ABORT: only {inter}/{len(keep_ids)} AMR IIDs found in {pfile}.psam — "
                         f"ID convention mismatch, not a filtering problem.")

    VARID = ["--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "100", "missing"]

    local = f"{SCRATCH}/amr_{CHROM}"
    if os.path.exists(local + ".pgen"):
        print(f"local subset pgen already exists ({local}.pgen) — skipping the bucket pass", flush=True)
    else:
        run([PLINK, "--pfile", pfile, "--keep", keep] + VARID + ["--maf", str(a.maf), "--snps-only",
             "just-acgt", "--max-alleles", "2", "--make-pgen",
             "--memory", str(a.memory), "--threads", str(a.threads), "--out", local],
            "extract AMR common SNPs -> local pgen")
    nvar = sum(1 for ln in open(local + ".pvar") if not ln.startswith("#"))
    print(f"local subset: {nvar} common biallelic SNPs x {n} AMR samples", flush=True)

    run([PLINK, "--pfile", local, "--indep-pairwise", str(a.window), str(a.step), str(a.r2),
         "--memory", str(a.memory), "--threads", str(a.threads), "--out", f"{SCRATCH}/amr_prune"],
        "LD prune (local)")
    npr = sum(1 for _ in open(f"{SCRATCH}/amr_prune.prune.in"))
    print(f"pruned-in SNPs: {npr} of {nvar}", flush=True)
    if npr < a.min_snps:
        raise SystemExit(f"ABORT: only {npr} pruned SNPs (floor {a.min_snps}) — too few for stable "
                         f"within-AMR PCs. Noisy PCs under-adjust and would bias the rescue toward "
                         f"SURVIVES. Relax --r2 / shrink --window and re-run (the local pgen is cached, "
                         f"so this costs seconds, not another bucket pass).")

    run([PLINK, "--pfile", local, "--extract", f"{SCRATCH}/amr_prune.prune.in",
         "--pca", str(a.kpc), "approx", "--memory", str(a.memory), "--threads", str(a.threads),
         "--out", f"{SCRATCH}/amr_pca"], "PCA (local)")

    for ext in (".eigenvec", ".eigenval", ".log"):
        src = f"{SCRATCH}/amr_pca{ext}"
        if os.path.exists(src):
            shutil.copy2(src, f"{OUT}/amr_pca{ext}")
    shutil.copy2(f"{SCRATCH}/amr_prune.prune.in", f"{OUT}/amr_pca.prune.in")
    shutil.copy2(keep, f"{OUT}/amr_keep.txt")

    ev = [float(x) for x in open(f"{OUT}/amr_pca.eigenval").read().split()]
    tot = sum(ev)
    print(f"\n=== within-AMR PCA done: {a.kpc} PCs from {npr} pruned {CHROM} SNPs, n={n} ===")
    print("  PC   eigenval   %var(of the K retained)")
    for i, e in enumerate(ev[:10], 1):
        print(f"  PC{i:<3} {e:9.3f}   {100*e/tot:5.2f}%")
    nrow = sum(1 for _ in open(f"{OUT}/amr_pca.eigenvec")) - 1
    print(f"\neigenvec rows: {nrow} (expected ~{n})")
    print(f"-> {OUT}/amr_pca.eigenvec   [ROW-LEVEL PII: bucket only, never commit]")
    print(f"STEP60_AMR_PCA_DONE ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
