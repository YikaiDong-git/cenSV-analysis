#!/usr/bin/env python3
import os, sys, csv, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
import pgenlib
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S
from firthmodels import FirthLogisticRegression

INP = S.INP; OUTA = S.OUTA
HITS = f"{OUTA}/phewide_disease_full.tsv"
OUT = f"{OUTA}/phewide_disease_full_firthrefit.tsv"
MAXHITS = 100000
NWORK = 6
FIRTH_MAXIT = 25; FIRTH_PEN = 0.5

_DOS = None; _YV = None; _XCOV = None


def load_sv_dosages(loci, iidx):
    prefixes = S.chunk_prefixes(); want = set(loci); out = {}
    N = len(iidx); buf = np.empty(N, np.float64)
    for pfx in prefixes:
        ids = [l.split("\t")[2] for l in open(f"{pfx}.pvar") if not l.startswith("#")]
        local = {sid: i for i, sid in enumerate(ids) if sid in want}
        if not local:
            continue
        r = pgenlib.PgenReader(f"{pfx}.pgen".encode())
        for sid, v in local.items():
            r.read_dosages(v, buf); out[sid] = buf.astype(np.float32).copy()
        r.close()
    return out


def load_pheno_cols(phecodes, iidx):
    N = len(iidx)
    hdr = open(f"{INP}/disease.pheno").readline().rstrip("\n").split("\t")
    colof = {ph: hdr.index(ph) for ph in phecodes if ph in hdr}
    Y = {ph: np.full(N, -1, np.int8) for ph in colof}
    with open(f"{INP}/disease.pheno") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[1])
            if kk is None:
                continue
            for ph, c in colof.items():
                if p[c] == "0":
                    Y[ph][kk] = 0
                elif p[c] == "1":
                    Y[ph][kk] = 1
    return Y


def firth(Xcov, g, y, present):
    Xp = Xcov[present]; gp = g[present]; yp = y[present].astype(int)
    keep = np.array([np.ptp(Xp[:, j]) > 0 for j in range(Xp.shape[1])])
    Xf = np.column_stack([Xp[:, keep], gp]); vi = int(keep.sum())
    m = FirthLogisticRegression(fit_intercept=True, max_iter=FIRTH_MAXIT, penalty_weight=FIRTH_PEN)
    try:
        m.fit(Xf, yp)
    except Exception:
        return np.nan, np.nan, False
    beta = float(m.coef_[vi])
    try:
        m.lrt(vi); p = float(m.lrt_pvalues_[vi])
    except Exception:
        p = float(m.pvalues_[vi])
    return beta, p, bool(getattr(m, "converged_", True))


def _refit_one(task):
    idx, sid, ph, beta_screen = task
    g = _DOS.get(sid); y = _YV.get(ph)
    if g is None or y is None:
        return idx, (np.nan, np.nan, False, 0)
    present = y >= 0
    bf, pf, conv = firth(_XCOV, g, y, present)
    confirmed = int(np.isfinite(pf) and pf < 0.05 and np.sign(bf) == np.sign(beta_screen))
    return idx, (bf, pf, conv, confirmed)


def main():
    global _DOS, _YV, _XCOV
    t0 = time.time()
    if not os.path.exists(HITS):
        raise SystemExit(f"ABORT: {HITS} not found — run step30g first.")
    hits = list(csv.DictReader(open(HITS), delimiter="\t"))
    hits.sort(key=lambda r: float(r["chi2_p"]) if r["chi2_p"] not in ("", "nan") else 1.0)
    refit = hits[:MAXHITS]
    loci = sorted({r["sv_locus"] for r in refit}); phes = sorted({r["phecode"] for r in refit})
    print(f"{len(hits)} hits -> refitting {len(refit)} ({len(loci)} SVs x {len(phes)} phecodes) on {NWORK} workers",
          flush=True)

    iids, iidx = S.load_order(); X, _ = S.load_X(iidx); _XCOV = X[:, 1:]
    _DOS = load_sv_dosages(loci, iidx); print(f"loaded {len(_DOS)} SV dosages ({time.time()-t0:.0f}s)", flush=True)
    _YV = load_pheno_cols(phes, iidx); print(f"loaded {len(_YV)} pheno cols ({time.time()-t0:.0f}s)", flush=True)

    tasks = [(i, r["sv_locus"], r["phecode"], float(r["beta_logOR"])) for i, r in enumerate(refit)]
    results = [None] * len(refit)
    done = 0
    with Pool(NWORK) as pool:
        for idx, res in pool.imap_unordered(_refit_one, tasks, chunksize=16):
            results[idx] = res; done += 1
            if done % 1000 == 0:
                print(f"  refit {done}/{len(refit)} ({time.time()-t0:.0f}s)", flush=True)

    n_conf = 0
    with open(OUT, "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tbeta_score\tp_score\tq_bh_tier\t"
                "beta_firth\tOR_firth\tp_firth\tfirth_converged\tfirth_confirmed\n")
        for r, res in zip(refit, results):
            bf, pf, conv, confirmed = res
            n_conf += confirmed
            orf = np.exp(bf) if np.isfinite(bf) else np.nan
            f.write(f"{r['sv_locus']}\t{r['tier']}\t{r['phecode']}\t{r['label']}\t{r['cases']}\t{r['n']}\t"
                    f"{r['beta_logOR']}\t{r['chi2_p']}\t{r['q_bh_tier']}\t{'' if not np.isfinite(bf) else f'{bf:.4g}'}\t"
                    f"{'' if not np.isfinite(orf) else f'{orf:.4g}'}\t{'' if not np.isfinite(pf) else f'{pf:.3g}'}\t"
                    f"{int(conv)}\t{confirmed}\n")
    print(f"FIRTH_REFIT_DONE: {len(refit)} hits refit, {n_conf} Firth-confirmed (p<0.05 & concordant) "
          f"= {100*n_conf/max(len(refit),1):.1f}% -> {OUT} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
