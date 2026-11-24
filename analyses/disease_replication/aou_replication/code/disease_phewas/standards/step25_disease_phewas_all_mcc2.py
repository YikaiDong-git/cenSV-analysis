#!/usr/bin/env python3
import numpy as np, pandas as pd, csv, os, gzip, gc
from scipy.stats import norm

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
RV = f"{BASE}/aou_replication/results_vds"
OUT = f"{RV}/step20_disease_all"
PHENO = f"{OUT}/phenotypes_all_mcc2.tsv"
PSUM = f"{OUT}/phenotype_summary_all_mcc2.tsv"
BASELINE = f"{OUT}/disease_phewas_all_ancestry.tsv"
COVFILE = f"{HOME}/workspace/working/scripts/covariates/output/covariates_v3_all_ancestries.tsv"
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


def main():
    meta = {r["phenotype_id"]: r for r in csv.DictReader(open(PSUM), delimiter="\t")}
    ph = pd.read_csv(PHENO, sep="\t", dtype=str)
    diseases = [c for c in ph.columns if c != "IID"]
    pho_iids = ph["IID"].tolist(); pidx = {i: k for k, i in enumerate(pho_iids)}; P = len(pho_iids)
    phmat = {d: pd.to_numeric(ph[d], errors="coerce").to_numpy(np.float32) for d in diseases}
    del ph; gc.collect()
    print(f"{len(diseases)} diseases, {P} ALL-ANCESTRY phenotyped IIDs (MCC=2)", flush=True)

    dosarr = {}; have = set()
    for path in DOSFILES:
        with gzip.open(path, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            newcols = [c for c in header[1:] if c not in have]
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

    base = ["age_at_recruitment_spline_1", "age_at_recruitment_spline_2", "age_at_recruitment_spline_3",
            "sex", *[f"gPC{i}" for i in range(1, 17)]]
    use = ["IID", *base, "ancestry_pred"]
    cov = pd.read_csv(COVFILE, sep="\t", dtype={"IID": str}, usecols=lambda c: c in use)
    cov["sex"] = cov["sex"].map({"Male": 1.0, "Female": 0.0})
    for c in base:
        cov[c] = pd.to_numeric(cov[c], errors="coerce")
    cov = cov.dropna(subset=base)
    cov = cov[cov["IID"].isin(pidx)].reset_index(drop=True)
    sel = np.array([pidx[i] for i in cov["IID"]], dtype=np.int64)
    Cbase = cov[base].to_numpy(np.float32)
    Cbase = np.hstack([Cbase, pd.get_dummies(cov["ancestry_pred"], drop_first=True).astype(np.float32).to_numpy()])
    anc_counts = cov["ancestry_pred"].value_counts().to_dict()
    del cov; gc.collect()
    print(f"N(all-anc, complete-cov, phenotyped) = {len(sel)}  ancestry: {anc_counts}", flush=True)

    HDR = ("sv_locus\tdisease\tlabel\tcategory\tsex_restriction\tmode\tn\tn_cases\tbeta_logOR\tOR\tse\tp\tq_bh\tconverged\n")
    rows = []
    for s in svs:
        dv = dosarr[s][sel]
        for d in diseases:
            yv = phmat[d][sel]
            mask = np.isfinite(dv) & np.isfinite(yv)
            n = int(mask.sum()); yb = yv[mask]; nc = int(yb.sum())
            if nc < MIN_CASE or (n - nc) < MIN_CASE:
                rows.append([s, d, meta[d]["label"], meta[d]["category"], meta[d]["sex_restriction"],
                             "all", n, nc, np.nan, np.nan, np.nan, np.nan, np.nan, "skip_lowcase"]); continue
            dm = dv[mask]
            if not (dm.std() > 1e-9):
                rows.append([s, d, meta[d]["label"], meta[d]["category"], meta[d]["sex_restriction"],
                             "all", n, nc, np.nan, np.nan, np.nan, np.nan, np.nan, "skip_degenerate"]); continue
            X = np.column_stack([np.ones(n, np.float64), dv[mask].astype(np.float64), Cbase[mask].astype(np.float64)])
            keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9
            X = X[:, keep]
            b, se, p, conv = logit_coef(X, yb.astype(np.float64))
            OR = float(np.exp(b)) if np.isfinite(b) else np.nan
            rows.append([s, d, meta[d]["label"], meta[d]["category"], meta[d]["sex_restriction"],
                         "all", n, nc, b, OR, se, p, np.nan, conv])
    qs = bh([r[11] for r in rows])
    for r, q in zip(rows, qs):
        r[12] = q
    outp = f"{OUT}/disease_phewas_all_mcc2.tsv"
    with open(outp, "w") as f:
        f.write(HDR)
        for r in rows:
            f.write("\t".join("" if (isinstance(x, float) and not np.isfinite(x)) else
                              (f"{x:.4g}" if isinstance(x, float) else str(x)) for x in r) + "\n")
    tested = [r for r in rows if r[13] not in ("skip_lowcase", "skip_degenerate")]
    sig = [r for r in tested if np.isfinite(r[12]) and r[12] < 0.05]
    print(f"\nN={len(sel)}, {len(tested)} tests, MCC=2 BH q<0.05: {len(sig)} -> {outp}", flush=True)
    for r in sorted(sig, key=lambda r: r[12])[:20]:
        print(f"  {r[0][:22]:22} {r[2][:30]:30} OR={r[9]:.3f} p={r[11]:.2e} q={r[12]:.3f}", flush=True)

    if os.path.exists(BASELINE):
        b1 = {(r["sv_locus"], r["disease"]): r for r in csv.DictReader(open(BASELINE), delimiter="\t")}
        now = {(r[0], r[1]): r for r in rows}
        def q(x):
            try:
                return float(x)
            except (TypeError, ValueError):
                return np.nan
        sig1 = {k for k, r in b1.items() if np.isfinite(q(r["q_bh"])) and q(r["q_bh"]) < 0.05}
        sig2 = {(r[0], r[1]) for r in sig}
        diffp = f"{OUT}/mcc2_vs_mcc1_hit_diff.tsv"
        with open(diffp, "w") as f:
            f.write("sv_locus\tdisease\tlabel\tq_mcc1\tq_mcc2\tOR_mcc1\tOR_mcc2\tn_cases_mcc1\tn_cases_mcc2\tstatus\n")
            for k in sorted(sig1 | sig2, key=lambda k: q(now.get(k, {12: "nan"})[12]) if k in now else 9):
                r2 = now.get(k); r1 = b1.get(k)
                st = ("both" if k in sig1 and k in sig2 else "MCC1_only" if k in sig1 else "MCC2_only")
                q1 = q(r1["q_bh"]) if r1 else np.nan
                q2 = r2[12] if r2 else np.nan
                or1 = q(r1["OR"]) if r1 else np.nan
                or2 = r2[9] if r2 else np.nan
                nc1 = r1["n_cases"] if r1 else ""
                nc2 = r2[7] if r2 else ""
                lab = (r2[2] if r2 else r1.get("label", "")) if (r2 or r1) else ""
                f.write(f"{k[0]}\t{k[1]}\t{lab}\t{q1:.4g}\t{q2:.4g}\t{or1:.4g}\t{or2:.4g}\t{nc1}\t{nc2}\t{st}\n")
        print(f"\nMCC=1 hits={len(sig1)}  MCC=2 hits={len(sig2)}  "
              f"both={len(sig1 & sig2)}  MCC1-only(lost)={len(sig1 - sig2)}  MCC2-only(new)={len(sig2 - sig1)}", flush=True)
        for tag in ("sv778", "sv121", "sv114"):
            for k in sorted(sig1 | sig2):
                if tag in k[0]:
                    st = ("both" if k in sig1 & sig2 else "MCC1_only" if k in sig1 else "MCC2_only")
                    print(f"    {k[0][:22]:22} {k[1]:20} {st}", flush=True)
        print(f"wrote {diffp}", flush=True)


if __name__ == "__main__":
    main()
