#!/usr/bin/env python3
import numpy as np, pandas as pd, csv, os, gzip, gc
from scipy.stats import norm

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
RV = f"{BASE}/aou_replication/results_vds"
OUT = f"{RV}/step16_disease"; os.makedirs(OUT, exist_ok=True)
PHENO = f"{HOME}/workspace/working/PRS/results/phewas/phenotypes.tsv"
PSUM = f"{HOME}/workspace/working/PRS/results/phewas/phenotype_summary.tsv"
COVDIR = f"{HOME}/workspace/working/scripts/covariates/output"
MODEFILE = {"eur": f"{COVDIR}/covariates_v3_eur.tsv", "all": f"{COVDIR}/covariates_v3_all_ancestries.tsv"}
DOSFILES = [f"{RV}/step6_222/sv_dosages_222.tsv.gz", f"{RV}/step2/sv_dosages.tsv.gz"]
MIN_CASE = 10


def logit_coef(X, y, maxit=100, tol=1e-7):
    beta = np.zeros(X.shape[1]); conv = False
    for _ in range(maxit):
        mu = 1 / (1 + np.exp(-(X @ beta))); W = np.clip(mu * (1 - mu), 1e-9, None)
        try:
            d = np.linalg.solve((X * W[:, None]).T @ X, X.T @ (y - mu))
        except np.linalg.LinAlgError:
            return np.nan, np.nan, np.nan, False
        beta += d
        if np.max(np.abs(d)) < tol:
            conv = True; break
    mu = 1 / (1 + np.exp(-(X @ beta)))
    try:
        se = np.sqrt(np.diag(np.linalg.inv((X * (mu * (1 - mu))[:, None]).T @ X)))
    except np.linalg.LinAlgError:
        return beta[1], np.nan, np.nan, conv
    return beta[1], se[1], 2 * norm.sf(abs(beta[1] / se[1])), conv


def bh(ps):
    ps = np.asarray(ps, float); n = int(np.isfinite(ps).sum()); order = np.argsort(np.where(np.isfinite(ps), ps, 2))
    q = np.full(len(ps), np.nan); prev = 1.0
    for rank, i in enumerate(order[::-1]):
        if not np.isfinite(ps[i]):
            continue
        prev = min(prev, ps[i] * n / (n - rank)); q[i] = prev
    return q


meta = {r["phenotype_id"]: r for r in csv.DictReader(open(PSUM), delimiter="\t")}

ph = pd.read_csv(PHENO, sep="\t", dtype=str)
diseases = [c for c in ph.columns if c != "IID"]
pho_iids = ph["IID"].tolist(); pidx = {i: k for k, i in enumerate(pho_iids)}; P = len(pho_iids)
phmat = {d: pd.to_numeric(ph[d], errors="coerce").to_numpy(np.float32) for d in diseases}
del ph; gc.collect()
print(f"{len(diseases)} diseases, {P} phenotyped IIDs (analysis universe)", flush=True)

dosarr = {}; have = set()
for path in DOSFILES:
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = header[1:]
        newcols = [c for c in cols if c not in have]
        for c in newcols:
            dosarr[c] = np.full(P, np.nan, np.float32)
        ci = {c: header.index(c) for c in newcols}
        for line in f:
            p = line.rstrip("\n").split("\t"); k = pidx.get(p[0])
            if k is None:
                continue
            for c in newcols:
                dosarr[c][k] = float(p[ci[c]])
        have.update(newcols)
svs = list(dosarr.keys())
print(f"{len(svs)} SVs (union, streamed onto {P} pheno IIDs)", flush=True)

HDR = ("sv_locus\tdisease\tlabel\tcategory\tsex_restriction\tmode\tn\tn_cases\tbeta_logOR\tOR\tse\tp\tq_bh\tconverged\n")
for mode, covfile in MODEFILE.items():
    base = ["age_at_recruitment_spline_1", "age_at_recruitment_spline_2", "age_at_recruitment_spline_3",
            "sex", *[f"gPC{i}" for i in range(1, 17)]]
    use = ["IID", *base] + (["ancestry_pred"] if mode == "all" else [])
    cov = pd.read_csv(covfile, sep="\t", dtype={"IID": str}, usecols=lambda c: c in use)
    cov["sex"] = cov["sex"].map({"Male": 1.0, "Female": 0.0})
    for c in base:
        cov[c] = pd.to_numeric(cov[c], errors="coerce")
    cov = cov.dropna(subset=base)
    cov = cov[cov["IID"].isin(pidx)].reset_index(drop=True)
    sel = np.array([pidx[i] for i in cov["IID"]], dtype=np.int64)
    Cbase = cov[base].to_numpy(np.float32)
    if mode == "all":
        Cbase = np.hstack([Cbase, pd.get_dummies(cov["ancestry_pred"], drop_first=True).astype(np.float32).to_numpy()])
    del cov; gc.collect()
    Nmode = len(sel)
    rows = []
    for s in svs:
        dv = dosarr[s][sel]
        for d in diseases:
            yv = phmat[d][sel]
            mask = np.isfinite(dv) & np.isfinite(yv)
            n = int(mask.sum()); yb = yv[mask]; nc = int(yb.sum())
            if nc < MIN_CASE or (n - nc) < MIN_CASE:
                rows.append([s, d, meta[d]["label"], meta[d]["category"], meta[d]["sex_restriction"],
                             mode, n, nc, np.nan, np.nan, np.nan, np.nan, np.nan, "skip_lowcase"]); continue
            dm = dv[mask]
            if not (dm.std() > 1e-9):
                rows.append([s, d, meta[d]["label"], meta[d]["category"], meta[d]["sex_restriction"],
                             mode, n, nc, np.nan, np.nan, np.nan, np.nan, np.nan, "skip_degenerate"]); continue
            X = np.column_stack([np.ones(n, np.float64), dv[mask].astype(np.float64), Cbase[mask].astype(np.float64)])
            keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9
            X = X[:, keep]
            b, se, p, conv = logit_coef(X, yb.astype(np.float64))
            OR = float(np.exp(b)) if np.isfinite(b) else np.nan
            rows.append([s, d, meta[d]["label"], meta[d]["category"], meta[d]["sex_restriction"],
                         mode, n, nc, b, OR, se, p, np.nan, conv])
    qs = bh([r[11] for r in rows])
    for r, q in zip(rows, qs):
        r[12] = q
    outp = f"{OUT}/disease_phewas_vds_{mode}.tsv"
    with open(outp, "w") as f:
        f.write(HDR)
        for r in rows:
            f.write("\t".join("" if (isinstance(x, float) and not np.isfinite(x)) else
                              (f"{x:.4g}" if isinstance(x, float) else str(x)) for x in r) + "\n")
    tested = [r for r in rows if r[13] != "skip_lowcase"]
    sig = [r for r in tested if np.isfinite(r[12]) and r[12] < 0.05]
    print(f"mode={mode}: N={Nmode}, {len(tested)} tests, {len(rows)-len(tested)} skipped low-case, BH q<0.05: {len(sig)} -> {outp}", flush=True)
    del Cbase; gc.collect()
print("DONE")
