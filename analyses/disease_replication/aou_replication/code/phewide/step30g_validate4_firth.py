#!/usr/bin/env python3
import os, sys, time
import numpy as np
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S
from firthmodels import FirthLogisticRegression

INP = S.INP; OUTA = S.OUTA; PSUM = S.PSUM
N_SV = 40
FIRTH_MAXIT = 25; FIRTH_PEN = 0.5
CARRIER = 0.5


def pick_phenos(colnames):
    import csv
    rows = [r for r in csv.DictReader(open(PSUM), delimiter="\t") if r["phecode"] in set(colnames)]
    by = {r["phecode"]: (int(r["cases_mcc2"]), r.get("sex", "Both")) for r in rows}
    have = [(ph, by[ph][0]) for ph in colnames if ph in by]
    have.sort(key=lambda x: x[1])
    chosen = []
    for target in (60, 120, 300, 800, 2000, 8000, 30000, 90000):
        cand = min(have, key=lambda x: abs(x[1] - target))
        if cand[0] not in [c[0] for c in chosen]:
            chosen.append(cand)
    for ph, (nc, sx) in by.items():
        if sx != "Both" and ph in colnames and ph not in [c[0] for c in chosen] and nc >= 100:
            chosen.append((ph, nc)); break
    return [c[0] for c in chosen]


def firth_fit(Xcov, g, y, present, verbose=False):
    Xp = Xcov[present]; gp = g[present]; yp = y[present].astype(int)
    keep = np.array([np.ptp(Xp[:, j]) > 0 for j in range(Xp.shape[1])])
    Xf = np.column_stack([Xp[:, keep], gp]); var_index = int(keep.sum())
    m = FirthLogisticRegression(fit_intercept=True, max_iter=FIRTH_MAXIT, penalty_weight=FIRTH_PEN)
    try:
        m.fit(Xf, yp)
    except Exception as e:
        msg = str(e).lower()
        return None, None, ("separation" if ("separat" in msg or "singular" in msg or "converg" in msg) else "fail")
    beta = float(m.coef_[var_index])
    try:
        m.lrt(var_index); p = float(m.lrt_pvalues_[var_index])
    except Exception:
        p = float(m.pvalues_[var_index])
    status = "ok" if getattr(m, "converged_", True) else "separation"
    return beta, p, status


def self_checks(G, X):
    rng = np.random.default_rng(0); N = G.shape[1]
    y = (rng.random(N) < 0.05).astype(np.int8)
    g = G[0]
    Be1, Pv1, *_ = S.score_test_batch(g[None], (g * g)[None], X, y[:, None], min_case=5)
    g2 = 2.0 * (g - g.min()) / (np.ptp(g) + 1e-12)
    Be2, Pv2, *_ = S.score_test_batch(g2[None], (g2 * g2)[None], X, y[:, None], min_case=5)
    inv_ok = np.isclose(Pv1[0, 0], Pv2[0, 0], rtol=1e-4, equal_nan=True) and np.sign(Be1[0, 0]) == np.sign(Be2[0, 0])
    print(f"[self-check] affine [0,2] invariance: chi2-p {Pv1[0,0]:.4g} vs {Pv2[0,0]:.4g}  -> {'PASS' if inv_ok else 'FAIL'}")


def main():
    t0 = time.time()
    iids, iidx = S.load_order(); N = len(iids)
    X, _ = S.load_X(iidx); Xcov = X[:, 1:]
    sv_order = []; G = S.load_G(sv_order); M = len(sv_order); G2 = G * G
    self_checks(G, X)

    carr_frac = (G > CARRIER).mean(axis=1)
    spread = np.argsort(carr_frac)[np.linspace(0, M - 1, N_SV).astype(int)]
    known = [i for i, s in enumerate(sv_order) if s.startswith("chr11_sv")][:6]
    sv_idx = sorted(set(spread.tolist()) | set(known))
    sub_ids = [sv_order[i] for i in sv_idx]; Gs = G[sv_idx]; G2s = G2[sv_idx]
    print(f"grid {len(sv_idx)} SVs; carrier-frac pct {np.round(np.percentile(carr_frac[sv_idx],[0,50,100]),4)}",
          flush=True)

    hdr = open(f"{INP}/disease.pheno").readline().rstrip("\n").split("\t")
    phenos = pick_phenos(hdr[2:])
    colof = {ph: hdr.index(ph) for ph in phenos if ph in hdr}
    print(f"phenos: {phenos}", flush=True)
    Yv = {ph: np.full(N, -1, np.int8) for ph in colof}
    with open(f"{INP}/disease.pheno") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[1])
            if kk is None:
                continue
            for ph, c in colof.items():
                if p[c] == "0":
                    Yv[ph][kk] = 0
                elif p[c] == "1":
                    Yv[ph][kk] = 1

    rows = []
    for ph in colof:
        y = Yv[ph]; present = y >= 0; nca = int((y == 1).sum())
        res = S.score_test_pheno(Gs, G2s, X, y)
        if res is None:
            print(f"  {ph}: score-test skipped (cases={nca})", flush=True); continue
        beta_s, _, p_s, npr, _ = res
        case_mask = (y == 1)
        for a in range(len(sv_idx)):
            carr_cases = int((Gs[a][case_mask] > CARRIER).sum())
            bf, pf, status = firth_fit(Xcov, Gs[a], y, present)
            rows.append((sub_ids[a], ph, nca, carr_cases,
                         float(p_s[a]), float(beta_s[a]),
                         (float(pf) if pf is not None else np.nan),
                         (float(bf) if bf is not None else np.nan), status))
        nfail = sum(1 for r in rows if r[1] == ph and r[8] != "ok")
        print(f"  {ph}: cases={nca} firth ok={len(sub_ids)-nfail}/{len(sub_ids)} "
              f"(sep/fail={nfail}) ({time.time()-t0:.0f}s)", flush=True)

    out = f"{OUTA}/validation_firth_v4.tsv"
    with open(out, "w") as f:
        f.write("sv_locus\tphecode\tcases\tcarrier_cases\tp_score\tbeta_score\tp_firth\tbeta_firth\tfirth_status\n")
        for r in rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.3g}\t{r[5]:.4g}\t"
                    f"{'' if not np.isfinite(r[6]) else f'{r[6]:.3g}'}\t"
                    f"{'' if not np.isfinite(r[7]) else f'{r[7]:.4g}'}\t{r[8]}\n")

    both = [r for r in rows if r[8] == "ok" and r[4] > 0 and np.isfinite(r[6]) and r[6] > 0]
    print(f"\n=== VALIDATION v4: score test vs firthmodels Firth (LRT), {len(rows)} cells, "
          f"{len(both)} with a valid Firth fit ===", flush=True)
    nsep = sum(1 for r in rows if r[8] != "ok")
    print(f"Firth non-fit (separation/fail): {nsep}/{len(rows)} — by carrier_cases bin:", flush=True)
    for lo, hi in [(0, 5), (5, 20), (20, 100), (100, 10**9)]:
        s = [r for r in rows if lo <= r[3] < hi]
        f_ = sum(1 for r in s if r[8] != "ok")
        print(f"    carrier_cases [{lo},{hi}): {f_}/{len(s)} non-fit", flush=True)

    def ratio_report(name, sel):
        s = [r for r in both if sel(r)]
        if len(s) < 4:
            print(f"{name:42s} n={len(s):4d} (too few)"); return
        lr = np.array([np.log10(r[4] / r[6]) for r in s])
        nd = sum(1 for r in s if np.sign(r[5]) == np.sign(r[7]))
        antic = sum(1 for r in s if r[4] < 0.5 * r[6])
        print(f"{name:42s} n={len(s):4d} dir={100*nd/len(s):5.1f}% "
              f"med log10(ps/pf)={np.median(lr):+.2f} frac(ps<0.5pf)={100*antic/len(s):4.1f}%")

    print("\n-- p-ratio by carrier_cases-among-cases (log10(ps/pf)<0 = score anti-conservative) --", flush=True)
    ratio_report("carrier_cases [0,5)", lambda r: r[3] < 5)
    ratio_report("carrier_cases [5,20)", lambda r: 5 <= r[3] < 20)
    ratio_report("carrier_cases [20,100)", lambda r: 20 <= r[3] < 100)
    ratio_report("carrier_cases >=100", lambda r: r[3] >= 100)
    ratio_report("ALL", lambda r: True)

    print("\n-- SCREEN PRECISION: among score-test positives, does Firth confirm? --", flush=True)
    for thr in (1e-3, 1e-4, 1e-5):
        pos = [r for r in both if r[4] < thr]
        if not pos:
            print(f"  score p<{thr:.0e}: 0 positives"); continue
        conf = sum(1 for r in pos if r[6] < 0.05)
        conf_dir = sum(1 for r in pos if r[6] < 0.05 and np.sign(r[5]) == np.sign(r[7]))
        print(f"  score p<{thr:.0e}: {len(pos):4d} positives, {100*conf/len(pos):.1f}% Firth-confirmed (p<0.05), "
              f"{100*conf_dir/max(len(pos),1):.1f}% confirmed AND concordant", flush=True)
    print(f"\nwrote {out} ({time.time()-t0:.0f}s)\nVALIDATION_DONE", flush=True)


if __name__ == "__main__":
    main()
