#!/usr/bin/env python3
"""M3 - per-stratum finemap keeps + within-stratum covar files."""
import os, csv

AR = "/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication"
FAM = f"{AR}/results_vds/phewide_carryover/sv.fam"
COV3 = os.path.expanduser("~/workspace/working/scripts/covariates/output/covariates_v3_all_ancestries.tsv")
COVBUILT = f"{AR}/results_vds/phewide_carryover/inputs/covar.txt"
REL = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/relatedness/samples_relatedness_flagged_samples.tsv")
OUTK = f"{AR}/finemap/inputs/keep"; OUTC = f"{AR}/finemap/inputs/covar"
os.makedirs(OUTK, exist_ok=True); os.makedirs(OUTC, exist_ok=True)
STRATA = ["all-anc", "EUR", "AMR"]
ANC_OF = {"EUR": "eur", "AMR": "amr"}

cohort = []
with open(FAM) as f:
    for line in f:
        p = line.split()
        if len(p) >= 2: cohort.append(p[1])
cohort_set = set(cohort)
print(f"cohort (sv.fam): {len(cohort):,}", flush=True)

anc = {}
with open(COV3) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        iid = row["IID"]
        if iid in cohort_set:
            anc[iid] = (row.get("ancestry_pred") or "").strip().lower()

drop = set()
with open(REL) as f:
    next(f)
    for line in f:
        drop.add(line.strip())
print(f"related flagged: {len(drop):,}", flush=True)

def in_stratum(iid, s):
    return True if s == "all-anc" else anc.get(iid) == ANC_OF[s]

summary = []
memb = {}
for s in STRATA:
    full = [i for i in cohort if in_stratum(i, s)]
    unrel = [i for i in full if i not in drop]
    memb[s] = set(full)
    with open(f"{OUTK}/{s}.keep.fid_iid.tsv", "w") as f:
        f.writelines(f"{i}\t{i}\n" for i in unrel)
    with open(f"{OUTK}/{s}.keep.iid.txt", "w") as f:
        f.writelines(f"{i}\n" for i in unrel)
    with open(f"{OUTK}/{s}.keep_full.fid_iid.tsv", "w") as f:
        f.writelines(f"{i}\t{i}\n" for i in full)
    summary.append((s, len(full), len(unrel), len(full) - len(unrel)))

with open(COVBUILT) as f:
    header = f.readline().rstrip("\n").split("\t")
    anc_cols = {i for i, c in enumerate(header) if c.startswith("anc_")}
    keep_idx = [i for i in range(len(header)) if i not in anc_cols]
    iid_col = header.index("IID")
    fh = {}
    for s in STRATA:
        fh[s] = open(f"{OUTC}/{s}.covar.txt", "w")
        cols = header if s == "all-anc" else [header[i] for i in keep_idx]
        fh[s].write("\t".join(cols) + "\n")
    for line in f:
        p = line.rstrip("\n").split("\t")
        iid = p[iid_col]
        for s in STRATA:
            if iid in memb[s]:
                fh[s].write((line if line.endswith("\n") else line + "\n") if s == "all-anc"
                            else "\t".join(p[i] for i in keep_idx) + "\n")
    for h in fh.values(): h.close()
    with open(f"{OUTC}/covar_names.peranc.txt", "w") as f2:
        f2.write(",".join(c for c in header[2:] if not c.startswith("anc_")))
    with open(f"{OUTC}/covar_names.allanc.txt", "w") as f2:
        f2.write(",".join(header[2:]))

print("\nstratum   N_full   N_unrel  N_related_dropped")
for s, nf, nu, nd in summary:
    print(f"{s:9} {nf:8} {nu:8} {nd:8}")
print("\nwrote keeps -> finemap/inputs/keep/  |  covars -> finemap/inputs/covar/")
