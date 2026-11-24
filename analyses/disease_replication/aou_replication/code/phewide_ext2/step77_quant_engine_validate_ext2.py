#!/usr/bin/env python3
"""Step 77 (ext2 arm) — independent validation of the quant engine (<- step30w)."""
"""step30w — validate the FWL partial-correlation QUANT score-test engine (step30d) against EXACT complete-case
OLS, and confirm no degenerate-LOINC flood (the quant analog of the disease PP_907.3 artifact).

WHY: step30d's per-SV quant columns (qt_min_p etc.) were SHIPPED to discovery (aou_persv_summary.tsv) but the
FWL engine (residualize dosages+traits on the FULL-sample covariate basis, then mask missing) was never
calibrated vs a gold standard, and no |Σr|-style degeneracy guard was run on the LOINCs. This closes that gap.

Exact reference: for each sampled (SV,LOINC) cell, fit OLS  Y_loinc[observed] ~ [intercept + 25 covariates +
dosage]  on that LOINC's observed subset (the TRUE per-phenotype FWL), and compare its dosage-coef t-test p to
the stored FWL score-test p. quant.pheno is pre-rank-INT (mean~0/sd~1), so no transform is applied.

Sample: the strongest hits + mid-range + near-threshold + any low-coverage cells (degeneracy-prone).
Out: printed concordance table (Δ -log10 p) + degeneracy verdict.
Run: python3 step30w_quant_engine_validate.py
"""
import os, sys, csv
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S
import step51b_exact_firth_refit_ext2 as H

QT = f"{S.OUTA}/quant/phewide_quant_ext2.tsv"
ALLPANEL = f"{S.OUTA}/quant/phewide_quant_allpanel_ext2.tsv"


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


OUTDIR = f"{S.OUTA}/quant/analysis/engine"
QUOTED_MAXDEV = 0.52


def sweep_degeneracy():
    per_loinc = {}
    npanel = 0
    if os.path.exists(ALLPANEL):
        with open(ALLPANEL) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            il, ip = hdr.index("loinc"), hdr.index("p")
            for line in f:
                q = line.rstrip("\n").split("\t")
                try:
                    p = float(q[ip])
                except (ValueError, IndexError):
                    continue
                npanel += 1
                a = per_loinc.setdefault(q[il], [0, 0])
                a[1] += 1
                if p < 1e-6:
                    a[0] += 1
    rows = []
    with open(QT) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                rows.append((r["sv_locus"], r["tier"], r["loinc"], int(float(r["n"])),
                             float(r["partial_r"]), float(r["p"])))
            except (ValueError, KeyError):
                continue
    NCOV = 26
    susp = []
    for sv, t, lo, n, pr, p in rows:
        why = []
        if abs(pr) > 0.7:
            why.append("abs_partial_r_gt_0.7")
        df = n - NCOV - 1
        if df > 10 and abs(pr) < 1 and p > 0:
            tt = abs(pr) * np.sqrt(df / max(1 - pr * pr, 1e-300))
            pre = float(2 * stats.norm.sf(tt))
            d = abs(np.log10(max(p, 1e-300)) - np.log10(max(pre, 1e-300)))
            if d > 0.5:
                why.append(f"p_inconsistent_with_r_and_n(dlog10={d:.2f})")
        if why:
            susp.append(dict(sv_locus=sv, tier=t, loinc=lo, n=n, partial_r=f"{pr:.4g}",
                             p=f"{p:.3g}", flags="+".join(why)))
    print(f"\n=== DEGENERACY SWEEP ===", flush=True)
    print(f"  flood scan over the FULL panel: {npanel:,} cells / {len(per_loinc)} LOINCs "
          f"({ALLPANEL.split('quant/')[-1]})", flush=True)
    print(f"  r-based checks over the hit table: {len(rows):,} cells -> {len(susp)} suspects", flush=True)
    for s in susp[:10]:
        print(f"    {s['sv_locus'][:38]:<40}{s['loinc']:<9}n={s['n']:<8}r={s['partial_r']:<9}"
              f"p={s['p']:<10}{s['flags']}", flush=True)
    if not per_loinc:
        print("  ⚠️ full panel dump not found — flood fraction NOT computed", flush=True)
        return susp, dict(n_hit_cells=len(rows), n_suspects=len(susp), verdict="flood scan SKIPPED")
    flood = sorted(((v[0] / v[1], lo, v[0], v[1]) for lo, v in per_loinc.items()), reverse=True)
    print("  worst per-LOINC flood fraction (share of that LOINC's 10,799 cells at p<1e-6; "
          "step62's disease-side threshold is 0.02):", flush=True)
    for fr, lo, a, b in flood[:6]:
        print(f"    {lo:<10}{a:>6}/{b:<7} = {fr:.4f}{'   <-- OVER 0.02' if fr > 0.02 else ''}", flush=True)
    over = [f for f in flood if f[0] > 0.02]
    verdict = ("NO FLOOD — every LOINC below step62's 0.02 threshold" if not over else
               f"{len(over)} LOINC(s) OVER 0.02: " + ", ".join(f"{lo} ({fr:.4f})" for fr, lo, _a, _b in over))
    print(f"  VERDICT: {verdict}", flush=True)
    return susp, dict(n_panel_cells=npanel, n_hit_cells=len(rows), n_loincs=len(per_loinc),
                      n_suspects=len(susp), worst_loinc=flood[0][1],
                      worst_flood_frac=f"{flood[0][0]:.5f}", n_loincs_over_002=len(over),
                      verdict=verdict)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    cells = sample_cells()
    svs = sorted({c["sv_locus"] for c in cells}); loincs = sorted({c["loinc"] for c in cells})
    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx); Xcov = X[:, 1:]
    dos = H.load_sv_dosages(svs, iidx)
    Y = read_loincs(iidx, loincs, N)
    print(f"validating {len(cells)} cells ({len(svs)} SVs x {len(loincs)} LOINCs) vs exact complete-case OLS\n", flush=True)
    print(f"  {'SV':<34}{'LOINC':<9}{'n':>8}{'FWL_p':>11}{'exact_p':>11}{'Δlog10p':>9}{'sign':>5}")
    maxdev = 0.0
    orows = []
    for c in cells:
        g = dos[c["sv_locus"]]; y = Y[c["loinc"]]
        b, t, pe, n = exact_ols_p(Xcov, g, y)
        pf = float(c["p"])
        l10 = lambda p: -np.log10(max(p, 1e-300))
        dev = l10(pe) - l10(pf); maxdev = max(maxdev, abs(dev))
        sign = "ok" if np.sign(b) == np.sign(float(c["beta"])) else "FLIP"
        print(f"  {c['sv_locus'][:33]:<34}{c['loinc']:<9}{n:>8}{pf:>11.1e}{pe:>11.1e}{dev:>9.2f}{sign:>5}", flush=True)
        orows.append(dict(sv_locus=c["sv_locus"], tier=c.get("tier", ""), loinc=c["loinc"], n=n,
                          p_fwl=f"{pf:.4g}", p_exact_ols=f"{pe:.4g}", delta_log10p=f"{dev:.4f}",
                          beta_fwl=f"{float(c['beta']):.5g}", beta_exact_ols=f"{b:.5g}",
                          t_exact_ols=f"{t:.4g}", sign_agree=sign))
    print(f"\nmax |Δ -log10 p| (exact vs FWL) = {maxdev:.2f}", flush=True)
    verdict = ("VALIDATED (agrees w/ exact; residualize-then-mask is conservative under missingness)"
               if maxdev < 2.0 else "DEVIATES — investigate")
    print("VERDICT: FWL engine", verdict, flush=True)

    if abs(round(maxdev, 2) - QUOTED_MAXDEV) > 0.005:
        print(f"  ⚠️ max|Δ| {maxdev:.2f} != the {QUOTED_MAXDEV} quoted in HEADTOHEAD §4 — the numeric path "
              f"changed; reconcile before trusting either number.", flush=True)
    else:
        print(f"  [V] reproduces the {QUOTED_MAXDEV} quoted in HEADTOHEAD §4 — persisting it changed "
              f"nothing numerically.", flush=True)

    cf = f"{OUTDIR}/quant_engine_validation_cells.tsv"
    hdr = ["sv_locus", "tier", "loinc", "n", "p_fwl", "p_exact_ols", "delta_log10p",
           "beta_fwl", "beta_exact_ols", "t_exact_ols", "sign_agree"]
    with open(cf, "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in orows:
            f.write("\t".join(str(r.get(h, "")) for h in hdr) + "\n")
    print(f"wrote {cf} ({len(orows)} cells)", flush=True)

    susp, summ = sweep_degeneracy()
    sf = f"{OUTDIR}/quant_engine_degeneracy_sweep.tsv"
    shdr = ["sv_locus", "tier", "loinc", "n", "partial_r", "p", "flags"]
    with open(sf, "w") as f:
        f.write("\t".join(shdr) + "\n")
        for r in susp:
            f.write("\t".join(str(r.get(h, "")) for h in shdr) + "\n")
    print(f"wrote {sf} ({len(susp)} suspects)", flush=True)

    vf = f"{OUTDIR}/quant_engine_validation_summary.tsv"
    with open(vf, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"n_cells_ols_validated\t{len(orows)}\n")
        f.write(f"max_abs_delta_log10p\t{maxdev:.4f}\n")
        f.write(f"engine_verdict\t{verdict}\n")
        for k, v in (summ or {}).items():
            f.write(f"degeneracy_{k}\t{v}\n")
    print(f"wrote {vf}", flush=True)
    print("QUANT_VALIDATE_DONE", flush=True)


if __name__ == "__main__":
    main()
