#!/usr/bin/env python3
"""M2b - persist a compact target-phecode phenotype (both MCC defs) for finemap."""
import os, csv, subprocess

AR = "/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication"
MAN = f"{AR}/finemap/targets/target_manifest.tsv"
SRC = {"mcc1": "/tmp/phewide_mcc1/inputs/disease.pheno",
       "mcc2": f"{AR}/results_vds/phewide_carryover/inputs/disease.pheno"}
OUT = f"{AR}/finemap/inputs/pheno"; os.makedirs(OUT, exist_ok=True)

targets = sorted({r["phecode"] for r in csv.DictReader(open(MAN), delimiter="\t")})
print(f"target phecodes: {len(targets)}", flush=True)
for mcc, src in SRC.items():
    with open(src) as f:
        hdr = f.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}
    present = [p for p in targets if p in idx]
    missing = [p for p in targets if p not in idx]
    cols = [idx["FID"] + 1, idx["IID"] + 1] + [idx[p] + 1 for p in present]
    out = f"{OUT}/disease.{mcc}.pheno"
    with open(out, "w") as fo:
        subprocess.run(["cut", "-f", ",".join(map(str, cols)), src], stdout=fo, check=True)
    print(f"{mcc}: {len(present)}/{len(targets)} target phecodes present  (missing {len(missing)}: {missing[:10]})  -> {out}", flush=True)
print("PHENO_EXTRACT_DONE", flush=True)
