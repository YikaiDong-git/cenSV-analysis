#!/usr/bin/env python3
"""Step 28b (phewide P1c): import sv.dose into plink2 pgen(s)."""
import os, subprocess, glob

HOME = os.path.expanduser("~")
SCRATCH = ("/tmp/claude-1000/-home-dataproc-workspace-working-ydon268/"
           "86d47bc2-43bf-4157-b8ed-79f770cbd7e0/scratchpad")
CHUNKDIR = f"{SCRATCH}/dose_chunks"
FAM = f"{SCRATCH}/phewide_pgen/sv.fam"
OUT = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds/step28_phewide/pgen"
PLINK2 = f"{HOME}/.local/bin/plink2"
os.makedirs(OUT, exist_ok=True)


def main():
    chunks = sorted(glob.glob(f"{CHUNKDIR}/chunk_*"))
    print(f"{len(chunks)} dose chunks", flush=True)
    done = []
    for i, ch in enumerate(chunks):
        outp = f"{OUT}/svpgen_{i}"
        if os.path.exists(outp + ".pgen") and not os.path.exists(outp + "-temporary.pgen"):
            print(f"  chunk {i}: already imported, skip", flush=True); done.append(outp); continue
        subprocess.run([PLINK2, "--import-dosage", ch, "noheader", "format=1", "single-chr=1",
                        "--fam", FAM, "--make-pgen", "--out", outp, "--memory", "640"], check=True)
        if os.path.exists(outp + "-temporary.pgen"):
            raise RuntimeError(f"chunk {i} import INCOMPLETE (temp remains)")
        print(f"  chunk {i} -> {outp}.pgen ({os.path.getsize(outp+'.pgen')/1e6:.0f} MB)", flush=True)
        done.append(outp)
    print(f"PGEN_DONE: {len(done)} chunk pgens -> {OUT}/svpgen_*.pgen", flush=True)


if __name__ == "__main__":
    main()
