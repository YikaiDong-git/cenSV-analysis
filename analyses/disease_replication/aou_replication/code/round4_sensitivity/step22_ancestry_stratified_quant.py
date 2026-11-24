#!/usr/bin/env python3
import os, gc, gzip, csv, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "round3_conditioning", "vds"))
from step10_snp_conditioning_vds import irnt, ols_dose  # noqa: E402

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
R6 = f"{BASE}/aou_replication/results_vds/step6_222"
COVDIR = f"{HOME}/workspace/working/scripts/covariates/output"
OUT = f"{BASE}/aou_replication/results_vds/step22_ancestry"; os.makedirs(OUT, exist_ok=True)

HITS = [("chr3_sv20_DEL_non_best6_alt_edge1000", "weight", "29463-7"),
        ("chr4_sv349_DEL_non_best8_fksize200k", "height", "8302-2"),
        ("chr11_sv263_DEL_non_best8_edge500", "BMI", "39156-5"),
        ("chr4_sv398_INS_non_best8_fksize100k", "MPV", "32623-1"),
        ("chr4_sv398_INS_non_best8_fksize100k", "PDW", "776-5")]
SVS = sorted({h[0] for h in HITS}); LOINCS = sorted({h[2] for h in HITS})
ANC = ["eur", "afr", "amr", "eas"]
GPC = [f"gPC{i}" for i in range(1, 17)]
BASECOV = ["age_at_recruitment_spline_1", "age_at_recruitment_spline_2",
           "age_at_recruitment_spline_3", "sex", *GPC]


def load_dosage_and_traits():
    iid2idx = {}
    dcols = {s: [] for s in SVS}
    with gzip.open(f"{R6}/sv_dosages_222.tsv.gz", "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t"); col = {s: hdr.index(s) for s in SVS}
        for line in f:
            p = line.rstrip("\n").split("\t")
            iid2idx[p[0]] = len(iid2idx)
            for s in SVS:
                dcols[s].append(np.float32(p[col[s]]))
    N = len(iid2idx)
    dose = {s: np.asarray(dcols[s], np.float32) for s in SVS}; del dcols; gc.collect()
    trait = {lc: np.full(N, np.nan, np.float32) for lc in LOINCS}
    for r in csv.DictReader(open(f"{R6}/measurement_medians_222.tsv"), delimiter="\t"):
        j = iid2idx.get(r["IID"])
        if j is not None and r["loinc"] in trait:
            trait[r["loinc"]][j] = np.float32(r["median_value"])
    return iid2idx, dose, trait


def fit_anc(anc, iid2idx, dose, trait):
    import pandas as pd
    covfile = f"{COVDIR}/covariates_v3_{anc}.tsv"
    if not os.path.exists(covfile):
        return {}
    use = ["IID", *BASECOV]
    cov = pd.read_csv(covfile, sep="\t", dtype={"IID": str}, usecols=lambda c: c in use)
    cov["sex"] = cov["sex"].map({"Male": 1.0, "Female": 0.0})
    for c in BASECOV:
        cov[c] = pd.to_numeric(cov[c], errors="coerce")
    cov = cov.dropna(subset=BASECOV).reset_index(drop=True)
    gidx = np.array([iid2idx.get(i, -1) for i in cov["IID"]], np.int64)
    have = gidx >= 0
    C = cov[BASECOV].to_numpy(np.float64)[have]; gsel = gidx[have]; del cov; gc.collect()
    out = {}
    for (sv, tname, loinc) in HITS:
        d = dose[sv][gsel].astype(np.float64)
        y = irnt(trait[loinc][gsel].astype(np.float64))
        m = np.isfinite(d) & np.isfinite(y)
        if m.sum() < 200:
            out[(sv, tname)] = (int(m.sum()), np.nan, np.nan, np.nan); continue
        X = np.column_stack([np.ones(m.sum()), d[m], C[m]])
        keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9
        b, se, p = ols_dose(X[:, keep], y[m].astype(np.float64))
        out[(sv, tname)] = (int(m.sum()), b, se, p)
    del C; gc.collect()
    return out


def ivw(betas, ses):
    w = np.array([1 / s**2 for b, s in zip(betas, ses) if np.isfinite(b) and np.isfinite(s) and s > 0])
    bs = np.array([b for b, s in zip(betas, ses) if np.isfinite(b) and np.isfinite(s) and s > 0])
    if len(bs) < 1:
        return np.nan, np.nan, np.nan, np.nan, 0
    bM = float(np.sum(w * bs) / np.sum(w)); seM = float(np.sqrt(1 / np.sum(w)))
    Q = float(np.sum(w * (bs - bM) ** 2)); k = len(bs)
    I2 = max(0.0, (Q - (k - 1)) / Q * 100) if Q > 0 and k > 1 else 0.0
    from scipy.stats import norm
    pM = 2 * norm.sf(abs(bM / seM))
    return bM, seM, pM, I2, k


def main():
    print("loading dosages + traits (shared, array-based)...", flush=True)
    iid2idx, dose, trait = load_dosage_and_traits()
    print(f"  {len(iid2idx)} IIDs in dosage set", flush=True)
    per = {}
    for anc in ANC:
        print(f"fitting {anc}...", flush=True)
        per[anc] = fit_anc(anc, iid2idx, dose, trait)
        for (sv, tname), (n, b, se, p) in sorted(per[anc].items()):
            if np.isfinite(b):
                print(f"  {anc} {tname:7} n={n:6d} b={b:+.4f} se={se:.4f} p={p:.2e}", flush=True)
    with open(f"{OUT}/ancestry_stratified_quant.tsv", "w") as fo:
        fo.write("sv\ttrait\t" + "\t".join(f"{a}_n\t{a}_beta\t{a}_se\t{a}_p" for a in ANC)
                 + "\tmeta_beta\tmeta_se\tmeta_p\tI2_pct\tk_anc\tall_sign_concordant\n")
        for (sv, tname, loinc) in HITS:
            cells, betas, ses, sig_signs = [], [], [], []
            for a in ANC:
                n, b, se, p = per[a].get((sv, tname), (0, np.nan, np.nan, np.nan))
                cells.append(f"{n}\t{b:.5g}\t{se:.5g}\t{p:.3g}")
                betas.append(b); ses.append(se)
                if np.isfinite(b) and np.isfinite(p) and p < 0.05:
                    sig_signs.append(np.sign(b))
            bM, seM, pM, I2, k = ivw(betas, ses)
            if len(sig_signs) <= 1:
                conc = "na"
            elif len(set(sig_signs)) == 1:
                conc = "yes"
            else:
                conc = "DISCORDANT"
            fo.write(f"{sv}\t{tname}\t" + "\t".join(cells)
                     + f"\t{bM:.5g}\t{seM:.5g}\t{pM:.3g}\t{I2:.1f}\t{k}\t{conc}\n")
            print(f"META {tname:7} b={bM:+.4f} se={seM:.4f} p={pM:.2e} I2={I2:.0f}% k={k} concordant={conc}", flush=True)
    print(f"wrote {OUT}/ancestry_stratified_quant.tsv", flush=True)


if __name__ == "__main__":
    main()
