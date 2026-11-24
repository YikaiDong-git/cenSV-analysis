#!/usr/bin/env python3
"""Build PLINK2 inputs for the disease PheWAS cross-check (validates step16 numpy)."""
import numpy as np, pandas as pd, gzip, os, sys, csv
HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
RV = f"{BASE}/aou_replication/results_vds"
PHENO = f"{HOME}/workspace/working/PRS/results/phewas/phenotypes.tsv"
COVF = {"eur": f"{HOME}/workspace/working/scripts/covariates/output/covariates_v3_eur.tsv",
        "all": f"{HOME}/workspace/working/scripts/covariates/output/covariates_v3_all_ancestries.tsv"}
DOSFILES = [f"{RV}/step6_222/sv_dosages_222.tsv.gz", f"{RV}/step2/sv_dosages.tsv.gz"]
mode = sys.argv[1]; OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)

ph = pd.read_csv(PHENO, sep="\t", dtype=str)
diseases = [c for c in ph.columns if c != "IID"]
pidx = {i: k for k, i in enumerate(ph["IID"].tolist())}; P = len(pidx)

base = ["age_at_recruitment_spline_1", "age_at_recruitment_spline_2", "age_at_recruitment_spline_3",
        "sex", *[f"gPC{i}" for i in range(1, 17)]]
use = ["IID", *base] + (["ancestry_pred"] if mode == "all" else [])
cov = pd.read_csv(COVF[mode], sep="\t", dtype={"IID": str}, usecols=lambda c: c in use)
cov["sex"] = cov["sex"].map({"Male": 1.0, "Female": 0.0})
for c in base:
    cov[c] = pd.to_numeric(cov[c], errors="coerce")
cov = cov.dropna(subset=base)
cov = cov[cov["IID"].isin(pidx)].reset_index(drop=True)

dosarr = {}; have = set()
for path in DOSFILES:
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t"); cols = header[1:]
        newc = [c for c in cols if c not in have]
        for c in newc:
            dosarr[c] = np.full(P, np.nan, np.float32)
        ci = {c: header.index(c) for c in newc}
        for line in f:
            p = line.rstrip("\n").split("\t"); k = pidx.get(p[0])
            if k is None:
                continue
            for c in newc:
                dosarr[c][k] = float(p[ci[c]])
        have.update(newc)
svs = list(dosarr.keys())

covi = np.array([pidx[i] for i in cov["IID"]], dtype=np.int64)
keepmask = np.ones(len(cov), bool)
for s in svs:
    keepmask &= np.isfinite(dosarr[s][covi])
samp = cov[keepmask].reset_index(drop=True)
si = np.array([pidx[i] for i in samp["IID"]], dtype=np.int64)
iids = samp["IID"].tolist(); N = len(iids)
print(f"mode={mode}: {N} samples (pheno∩cov∩dosage), {len(svs)} SVs, {len(diseases)} diseases", flush=True)

sexmap = samp["sex"].map({1.0: "1", 0.0: "2"}).fillna("0")
with open(f"{OUT}/sv.fam", "w") as f:
    for iid, sx in zip(iids, sexmap):
        f.write(f"{iid}\t{iid}\t0\t0\t{sx}\t-9\n")

with open(f"{OUT}/sv.dose", "w") as f:
    for s in svs:
        d = dosarr[s][si].astype(np.float64)
        lo, hi = np.nanmin(d), np.nanmax(d)
        ds = 2.0 * (d - lo) / (hi - lo) if hi > lo else np.full(N, 1.0)
        f.write(f"{s}\tA\tC\t" + "\t".join(f"{x:.5f}" for x in ds) + "\n")

phs = ph.set_index("IID").loc[iids].reset_index()
with open(f"{OUT}/disease.pheno", "w") as f:
    f.write("FID\tIID\t" + "\t".join(diseases) + "\n")
    for _, r in phs.iterrows():
        vals = [(r[d] if r[d] in ("0", "1") else "NA") for d in diseases]
        f.write(f"{r['IID']}\t{r['IID']}\t" + "\t".join(vals) + "\n")

covcols = list(base)
C = samp[base].copy()
if mode == "all":
    dums = pd.get_dummies(samp["ancestry_pred"], prefix="anc", drop_first=True).astype(int)
    C = pd.concat([C.reset_index(drop=True), dums.reset_index(drop=True)], axis=1)
    covcols = list(C.columns)
with open(f"{OUT}/covar.txt", "w") as f:
    f.write("FID\tIID\t" + "\t".join(covcols) + "\n")
    for iid, (_, row) in zip(iids, C.iterrows()):
        f.write(f"{iid}\t{iid}\t" + "\t".join(f"{row[c]:.6g}" for c in covcols) + "\n")

with open(f"{OUT}/covar_names.txt", "w") as f:
    f.write(",".join(covcols))
print(f"wrote sv.dose ({len(svs)} SV rows), sv.fam, disease.pheno, covar.txt ({len(covcols)} covars) -> {OUT}", flush=True)
