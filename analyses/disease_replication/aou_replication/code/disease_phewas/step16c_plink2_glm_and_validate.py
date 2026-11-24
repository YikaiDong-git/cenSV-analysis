#!/usr/bin/env python3
"""Step 16c: run PLINK2 --glm on the SV-dosage pgen and validate vs the numpy step16."""
import subprocess, glob, os, io, tarfile, csv, math, numpy as np
HOME = os.path.expanduser("~")
RV = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds"
E = f"{RV}/step16_plink2/eur"; INP = f"{E}/inputs"; GLM = f"{E}/glm"
PLINK2 = f"{HOME}/.local/bin/plink2"
NP = f"{RV}/step16_disease/disease_phewas_vds_eur.tsv"
os.makedirs(GLM, exist_ok=True)

subprocess.run([PLINK2, "--import-dosage", f"{INP}/sv.dose", "noheader", "format=1", "single-chr=1",
                "--fam", f"{INP}/sv.fam", "--make-pgen", "--out", f"{INP}/svpgen", "--memory", "640"], check=True)
subprocess.run([PLINK2, "--pfile", f"{INP}/svpgen", "--glm", "hide-covar", "firth-fallback",
                "cols=+a1freq,+beta", "--pheno", f"{INP}/disease.pheno", "--1",
                "--covar", f"{INP}/covar.txt", "--covar-variance-standardize",
                "--vif", "100000", "--max-corr", "0.9999", "--out", f"{GLM}/glm_eur", "--memory", "640"], check=True)

def parse(files_iter):
    rows = {}
    for name, fh in files_iter:
        dis = name[len("glm_eur."):-len(".glm.logistic.hybrid")]
        for x in csv.DictReader(fh, delimiter="\t"):
            if x["TEST"] != "ADD":
                continue
            p = float(x["P"]) if x["P"] not in ("", "NA") else float("nan")
            rows[(x["ID"], dis)] = (int(x["OBS_CT"]), math.exp(float(x["BETA"])), p)
    return rows
live = ((os.path.basename(f), open(f)) for f in glob.glob(f"{GLM}/glm_eur.*.glm.logistic.hybrid"))
rows = parse(live)

ks = list(rows); ps = np.array([rows[k][2] for k in ks])
order = np.argsort(np.where(np.isfinite(ps), ps, 2)); n = int(np.isfinite(ps).sum()); prev = 1.0; q = {}
for rank, i in enumerate(order[::-1]):
    if not np.isfinite(ps[i]):
        continue
    prev = min(prev, ps[i] * n / (n - rank)); q[ks[i]] = prev

npd = {(x["sv_locus"], x["disease"]): float(x["p"]) for x in csv.DictReader(open(NP), delimiter="\t")
       if x["p"] not in ("", "NA")}
a, b = [], []
for k, (nn, orr, p) in rows.items():
    if k in npd and np.isfinite(p) and p > 0 and npd[k] > 0:
        a.append(math.log10(p)); b.append(math.log10(npd[k]))
corr = np.corrcoef(a, b)[0, 1]; mad = max(abs(x - y) for x, y in zip(a, b))
print(f"plink2 vs numpy: {len(a)} common tests, corr(log10 p)={corr:.6f}, max|Δlog10 p|={mad:.4f}")
print(f"BH q<0.05: {sum(1 for k in q if q[k] < 0.05)}")

meta = {(x["sv_locus"], x["disease"]): x.get("category", "") for x in csv.DictReader(open(NP), delimiter="\t")}
outp = f"{RV}/step16_disease/disease_phewas_COMBINED_eur.tsv"
with open(outp, "w") as f:
    f.write("sv_locus\tdisease\tcategory\tn\tOR_rescaled\tp_plink2\tq_bh\tflag\n")
    for k in sorted(rows, key=lambda k: (q.get(k, 9), rows[k][2])):
        nn, orr, p = rows[k]
        f.write(f"{k[0]}\t{k[1]}\t{meta.get(k,'')}\t{nn}\t{orr:.4g}\t{p:.3g}\t{q.get(k,float('nan')):.3g}\t"
                f"{'neg_ctrl' if 'negative_control' in k[1] else ''}\n")
print(f"wrote {outp}")
