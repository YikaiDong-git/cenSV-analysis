#!/usr/bin/env python3
import os, subprocess, glob

HOME = os.path.expanduser("~")
ROOT = "/tmp/phewide"
BUCK = (f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/"
        "aou_replication/results_vds/step28_phewide")
CHUNKDIR = f"{BUCK}/dose_chunks_tmp"
FAM = f"{ROOT}/phewide_pgen/sv.fam"
OUT = f"{BUCK}/pgen"
PLINK2 = f"{HOME}/.local/bin/plink2"
os.makedirs(OUT, exist_ok=True)


def main():
    chunks = sorted(glob.glob(f"{CHUNKDIR}/chunk_*"))
    print(f"{len(chunks)} dose chunks", flush=True)
    total_vars = 0
    for i, ch in enumerate(chunks):
        outp = f"{OUT}/svpgen300_{i}"
        if os.path.exists(outp + ".pgen") and not os.path.exists(outp + "-temporary.pgen"):
            nv = sum(1 for l in open(outp + ".pvar") if not l.startswith("#"))
            total_vars += nv
            print(f"  chunk {i}: already imported ({nv} vars), skip", flush=True); continue
        subprocess.run([PLINK2, "--import-dosage", ch, "noheader", "format=1", "single-chr=1",
                        "--fam", FAM, "--make-pgen", "--out", outp, "--memory", "8000"], check=True)
        if os.path.exists(outp + "-temporary.pgen"):
            raise RuntimeError(f"chunk {i} import INCOMPLETE (temp remains)")
        nv = sum(1 for l in open(outp + ".pvar") if not l.startswith("#"))
        total_vars += nv
        print(f"  chunk {i} -> {outp}.pgen ({nv} vars, {os.path.getsize(outp+'.pgen')/1e6:.0f} MB)", flush=True)
    assert len(chunks) == 22, f"EXPECTED 22 chunks, got {len(chunks)} — check `split -l 300`"
    assert total_vars == 6513, f"EXPECTED 6513 total variants, got {total_vars}"
    print(f"PGEN_DONE: {len(chunks)} chunks, {total_vars} variants -> {OUT}/svpgen300_*.pgen", flush=True)


if __name__ == "__main__":
    main()
