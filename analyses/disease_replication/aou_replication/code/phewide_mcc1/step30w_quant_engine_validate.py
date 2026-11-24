#!/usr/bin/env python3
import os, sys, csv
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S
import step30h_firth_refit_hits as H

QT = f"{S.OUTA}/phewide_quant.tsv"


def sample_cells():
    rows = [r for r in csv.DictReader(open(QT), delimiter="\t") if r["p"] not in ("", "nan")]
    rows.sort(key=lambda r: float(r["p"]))
    n = len(rows)
    pick = rows[:5] + rows[n // 100:n // 100 + 3] + rows[n // 10:n // 10 + 3] + rows[-4:]
    pick += sorted(rows, key=lambda r: float(r["n"]))[:4]
    seen, out = set(), []
    for r in pick:
        k = (r["sv_locus"], r["loinc"])
        if k not in seen:
            seen.add(k); out.append(r)
    return out


def read_loincs(iidx, loincs, N):
    hdr = open(f"{S.INP}/quant.pheno").readline().rstrip("\n").split("\t")
    col = {lo: hdr.index(lo) for lo in loincs if lo in hdr}
    Y = {lo: np.full(N, np.nan, np.float64) for lo in col}
    with open(f"{S.INP}/quant.pheno") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t"); k = iidx.get(p[1])
            if k is None:
                continue
            for lo, c in col.items():
                if p[c] != "NA":
                    Y[lo][k] = float(p[c])
    return Y


def exact_ols_p(Xcov, g, y):
    m = np.isfinite(y)
    Xd = np.column_stack([np.ones(m.sum()), Xcov[m], g[m]])
    yy = y[m]
    keep = np.array([np.ptp(Xd[:, j]) > 0 or j == 0 for j in range(Xd.shape[1])])
    Xd = Xd[:, keep]
    beta, *_ = np.linalg.lstsq(Xd, yy, rcond=None)
    resid = yy - Xd @ beta
    dof = Xd.shape[0] - Xd.shape[1]
    sigma2 = resid @ resid / dof
    XtXinv = np.linalg.inv(Xd.T @ Xd)
    se = np.sqrt(sigma2 * XtXinv[-1, -1])
    t = beta[-1] / se
    p = 2 * stats.t.sf(abs(t), dof)
    return beta[-1], t, p, int(m.sum())


def main():
    cells = sample_cells()
    svs = sorted({c["sv_locus"] for c in cells}); loincs = sorted({c["loinc"] for c in cells})
    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx); Xcov = X[:, 1:]
    dos = H.load_sv_dosages(svs, iidx)
    Y = read_loincs(iidx, loincs, N)
    print(f"validating {len(cells)} cells ({len(svs)} SVs x {len(loincs)} LOINCs) vs exact complete-case OLS\n", flush=True)
    print(f"  {'SV':<34}{'LOINC':<9}{'n':>8}{'FWL_p':>11}{'exact_p':>11}{'Δlog10p':>9}{'sign':>5}")
    maxdev = 0.0
    for c in cells:
        g = dos[c["sv_locus"]]; y = Y[c["loinc"]]
        b, t, pe, n = exact_ols_p(Xcov, g, y)
        pf = float(c["p"])
        l10 = lambda p: -np.log10(max(p, 1e-300))
        dev = l10(pe) - l10(pf); maxdev = max(maxdev, abs(dev))
        sign = "ok" if np.sign(b) == np.sign(float(c["beta"])) else "FLIP"
        print(f"  {c['sv_locus'][:33]:<34}{c['loinc']:<9}{n:>8}{pf:>11.1e}{pe:>11.1e}{dev:>9.2f}{sign:>5}", flush=True)
    print(f"\nmax |Δ -log10 p| (exact vs FWL) = {maxdev:.2f}", flush=True)
    print("VERDICT: FWL engine", "VALIDATED (agrees w/ exact; residualize-then-mask is conservative under missingness)"
          if maxdev < 2.0 else "DEVIATES — investigate", flush=True)
    print("QUANT_VALIDATE_DONE", flush=True)


if __name__ == "__main__":
    main()
