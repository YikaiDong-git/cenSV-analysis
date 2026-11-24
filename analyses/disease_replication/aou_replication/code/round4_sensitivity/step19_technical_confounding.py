#!/usr/bin/env python3
import os, gc, gzip, csv, numpy as np
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "round3_conditioning", "vds"))
from step10_snp_conditioning_vds import irnt  # noqa: E402

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
R6 = f"{BASE}/aou_replication/results_vds/step6_222"
COVDIR = f"{HOME}/workspace/working/scripts/covariates/output"
OUT = f"{BASE}/aou_replication/results_vds/step19_technical"; os.makedirs(OUT, exist_ok=True)

HITS = [
    ("chr3_sv20_DEL_non_best6_alt_edge1000", "weight", "29463-7"),
    ("chr4_sv349_DEL_non_best8_fksize200k", "height", "8302-2"),
    ("chr11_sv263_DEL_non_best8_edge500", "BMI", "39156-5"),
    ("chr4_sv398_INS_non_best8_fksize100k", "MPV", "32623-1"),
    ("chr4_sv398_INS_non_best8_fksize100k", "PDW", "776-5"),
]
SVS = sorted({h[0] for h in HITS})
LOINCS = sorted({h[2] for h in HITS})
TECH = ["depth", "coverage", "contamination"]
GPC = [f"gPC{i}" for i in range(1, 17)]
BASECOV = ["age_at_recruitment_spline_1", "age_at_recruitment_spline_2",
           "age_at_recruitment_spline_3", "sex", *GPC]


def ols_full(X, y):
    from scipy.stats import t as tdist
    XtXi = np.linalg.inv(X.T @ X); beta = XtXi @ (X.T @ y); resid = y - X @ beta
    n, p = X.shape; s2 = float(resid @ resid) / (n - p)
    se = np.sqrt(np.diag(XtXi) * s2)
    sst = float(((y - y.mean()) ** 2).sum()); ssr = float(resid @ resid)
    r2 = 1 - ssr / sst if sst > 0 else float("nan")
    tval = beta[1] / se[1]
    return beta[1], se[1], 2 * tdist.sf(abs(tval), n - p), r2


def r2_of(dose, M):
    m = np.isfinite(dose) & np.all(np.isfinite(M), axis=1)
    X = np.column_stack([np.ones(m.sum()), M[m]])
    keep = np.ones(X.shape[1], bool); keep[1:] = X[:, 1:].std(0) > 1e-9
    X = X[:, keep]; y = dose[m].astype(np.float64)
    b = np.linalg.lstsq(X, y, rcond=None)[0]; resid = y - X @ b
    sst = float(((y - y.mean()) ** 2).sum())
    return 1 - float(resid @ resid) / sst if sst > 0 else float("nan"), int(m.sum())


def main():
    import pandas as pd
    use = ["IID", *BASECOV, *TECH]
    cov = pd.read_csv(f"{COVDIR}/covariates_v3_eur.tsv", sep="\t", dtype={"IID": str},
                      usecols=lambda c: c in use)
    cov["sex"] = cov["sex"].map({"Male": 1.0, "Female": 0.0})
    for c in BASECOV + TECH:
        cov[c] = pd.to_numeric(cov[c], errors="coerce")
    cov = cov.dropna(subset=BASECOV + TECH).reset_index(drop=True)
    iids = cov["IID"].tolist(); idx = {i: k for k, i in enumerate(iids)}; N = len(iids)
    Cbase = cov[BASECOV].to_numpy(np.float64)
    Ctech = cov[TECH].to_numpy(np.float64)
    Cgpc = cov[GPC].to_numpy(np.float64)
    print(f"EUR N with complete tech covars = {N}", flush=True)

    dose = {s: np.full(N, np.nan, np.float32) for s in SVS}
    with gzip.open(f"{R6}/sv_dosages_222.tsv.gz", "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t"); col = {s: hdr.index(s) for s in SVS}
        for line in f:
            p = line.rstrip("\n").split("\t"); i = idx.get(p[0])
            if i is not None:
                for s in SVS:
                    dose[s][i] = float(p[col[s]])
    lv = {lc: np.full(N, np.nan, np.float32) for lc in LOINCS}
    for r in csv.DictReader(open(f"{R6}/measurement_medians_222.tsv"), delimiter="\t"):
        i = idx.get(r["IID"])
        if i is not None and r["loinc"] in lv:
            lv[r["loinc"]][i] = float(r["median_value"])
    lirnt = {lc: irnt(v.astype(float)) for lc, v in lv.items()}; del lv; gc.collect()

    rows = []
    print(f"\n{'hit':22} {'R2_seq':>8} {'R2_anc':>8} | {'b_base':>9} {'b_tech':>9} {'atten%':>7} {'p_base':>9} {'p_tech':>9}", flush=True)
    for (sv, tname, loinc) in HITS:
        d = dose[sv].astype(np.float64); y = lirnt[loinc]
        r2_seq, _ = r2_of(d, Ctech)
        r2_anc, _ = r2_of(d, Cgpc)
        m = np.isfinite(y) & np.isfinite(d)
        Xb = np.column_stack([np.ones(m.sum()), d[m], Cbase[m]])
        Xt = np.column_stack([np.ones(m.sum()), d[m], Cbase[m], Ctech[m]])
        bb, seb, pb, _ = ols_full(Xb[:, [True, True, *(Xb[:, 2:].std(0) > 1e-9)]], y[m].astype(np.float64))
        bt, set_, pt, _ = ols_full(Xt[:, [True, True, *(Xt[:, 2:].std(0) > 1e-9)]], y[m].astype(np.float64))
        atten = (abs(bb) - abs(bt)) / abs(bb) * 100 if bb != 0 else float("nan")
        rows.append([sv, tname, loinc, int(m.sum()), f"{r2_seq:.4f}", f"{r2_anc:.4f}",
                     f"{bb:.5g}", f"{pb:.3g}", f"{bt:.5g}", f"{pt:.3g}", f"{atten:.1f}"])
        print(f"{sv[:18]+'/'+tname:22} {r2_seq:8.4f} {r2_anc:8.4f} | {bb:+9.5f} {bt:+9.5f} {atten:7.1f} {pb:9.2e} {pt:9.2e}", flush=True)

    with open(f"{OUT}/technical_confounding_eur.tsv", "w") as fo:
        fo.write("sv\ttrait\tloinc\tn\tr2_seq_depth_cov_contam\tr2_ancestry_gpc16\t"
                 "beta_base\tp_base\tbeta_tech_adj\tp_tech_adj\tpct_atten_tech\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")
    print(f"\nwrote {OUT}/technical_confounding_eur.tsv", flush=True)


if __name__ == "__main__":
    main()
