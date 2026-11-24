#!/usr/bin/env python3
import os, csv, subprocess, math
from collections import defaultdict

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
SC = "/tmp/claude-1000/-home-dataproc-workspace-working-ydon268/86d47bc2-43bf-4157-b8ed-79f770cbd7e0/scratchpad"
MED = f"{SC}/phewide_evidence/quant_medians.tsv"
FAM = f"{SC}/phewide_pgen/sv.fam"
INP = f"{SC}/inputs"; os.makedirs(INP, exist_ok=True)
OUTA = f"{BASE}/aou_replication/results_vds/step30_phewide"; os.makedirs(OUTA, exist_ok=True)
CATL = f"{BASE}/aou_replication/results_vds/step21_catalog/loinc_catalog.tsv"
MIN_N = 500


def inv_norm_block(vals):
    n = len(vals)
    order = sorted(range(n), key=lambda i: vals[i][1])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and vals[order[j + 1]][1] == vals[order[i]][1]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    out = []
    for idx in range(n):
        q = (ranks[idx] - 0.5) / n
        out.append((vals[idx][0], _ppf(q)))
    return out


_a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
      1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
_b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
      6.680131188771972e+01, -1.328068155288572e+01]
_c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
      -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
_d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]


def _ppf(p):
    p = min(max(p, 1e-12), 1 - 1e-12)
    pl = 0.02425
    if p < pl:
        q = math.sqrt(-2 * math.log(p))
        return (((((_c[0]*q+_c[1])*q+_c[2])*q+_c[3])*q+_c[4])*q+_c[5]) / ((((_d[0]*q+_d[1])*q+_d[2])*q+_d[3])*q+1)
    if p > 1 - pl:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((_c[0]*q+_c[1])*q+_c[2])*q+_c[3])*q+_c[4])*q+_c[5]) / ((((_d[0]*q+_d[1])*q+_d[2])*q+_d[3])*q+1)
    q = p - 0.5; r = q * q
    return (((((_a[0]*r+_a[1])*r+_a[2])*r+_a[3])*r+_a[4])*r+_a[5])*q / (((((_b[0]*r+_b[1])*r+_b[2])*r+_b[3])*r+_b[4])*r+1)


def main():
    iids = [l.split()[1] for l in open(FAM)]
    fam_set = set(iids)
    name = {r["loinc"]: r["name"] for r in csv.DictReader(open(CATL), delimiter="\t")}

    med_by_loinc = f"{SC}/phewide_evidence/quant_by_loinc.tsv"
    if not os.path.exists(med_by_loinc):
        print("sorting medians by loinc ...", flush=True)
        with open(med_by_loinc, "w") as fo:
            subprocess.run(["sort", "-t", "\t", "-k2,2", "-S", "128M", "--parallel=2", MED],
                           stdout=fo, check=True)
    long_z = f"{SC}/phewide_evidence/quant_int_long.tsv"
    kept_counts = {}
    with open(med_by_loinc) as f, open(long_z, "w") as fo:
        cur = None; block = []
        def emit(loinc, block):
            vals = [(pid, v) for (pid, v) in block if pid in fam_set]
            if len(vals) < MIN_N:
                return
            for pid, z in inv_norm_block(vals):
                fo.write(f"{pid}\t{loinc}\t{z:.5f}\n")
            kept_counts[loinc] = len(vals)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "IID" or len(p) < 3:
                continue
            try:
                v = float(p[2])
            except ValueError:
                continue
            loinc = p[1]
            if loinc != cur:
                if cur is not None:
                    emit(cur, block)
                cur = loinc; block = []
            block.append((p[0], v))
        if cur is not None:
            emit(cur, block)
    kept = sorted(kept_counts)
    print(f"kept {len(kept)} LOINCs (>= {MIN_N} people)", flush=True)

    z_by_pid = f"{SC}/phewide_evidence/quant_int_by_pid.tsv"
    with open(z_by_pid, "w") as fo:
        subprocess.run(["sort", "-t", "\t", "-k1,1n", "-S", "128M", "--parallel=2", long_z],
                       stdout=fo, check=True)
    fam_set = set(iids)
    nwritten = 0
    with open(z_by_pid) as f, open(f"{INP}/quant.pheno", "w") as fo:
        fo.write("FID\tIID\t" + "\t".join(kept) + "\n")
        cur = None; vals = {}
        def emit(pid, vals):
            if pid not in fam_set:
                return 0
            fo.write(f"{pid}\t{pid}\t" + "\t".join(vals.get(lo, "NA") for lo in kept) + "\n")
            return 1
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            if p[0] != cur:
                if cur is not None:
                    nwritten += emit(cur, vals)
                cur = p[0]; vals = {}
            vals[p[1]] = p[2]
        if cur is not None:
            nwritten += emit(cur, vals)
    print(f"quant.pheno rows written (persons with >=1 measured LOINC in cohort): {nwritten}", flush=True)
    with open(f"{OUTA}/quant_pheno_summary.tsv", "w") as f:
        f.write("loinc\tname\tn_people\n")
        for lo in kept:
            f.write(f"{lo}\t{name.get(lo,'')}\t{kept_counts[lo]}\n")
    print(f"wrote quant.pheno ({len(iids)} samples x {len(kept)} LOINCs) + summary -> {OUTA}", flush=True)
    print("QUANT_PHENO_DONE", flush=True)


if __name__ == "__main__":
    main()
