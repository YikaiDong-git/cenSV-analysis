#!/usr/bin/env python3
"""Step 30h2 (DISEASE) — approximate-Firth recalibration of the score-test tail + separation handling."""
import os, sys, csv, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4"); os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.stats import chi2 as chi2dist, fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S

MATDIR = f"{S.OUTA}/disease_matrices"
OUTDIR = f"{S.OUTA}/disease_matrices_firth"; os.makedirs(OUTDIR, exist_ok=True)

PTHRESH       = 0.05
REFIT_CUT     = 1e-3
FIRTH_MAXIT   = 50
TOL           = 1e-7
MAXSTEP       = 5.0
SEP_BETA_BOUND= 10.0
CARRIER_THRESH= 0.5
EPS           = 1e-12
NWORK         = 14
BAND_BATCH    = 150


def load_score_matrices():
    Pv = np.load(f"{MATDIR}/disease_Pv.npy")
    Be = np.load(f"{MATDIR}/disease_Be.npy")
    svloc, tiers = [], []
    for ln in open(f"{MATDIR}/disease_svorder.txt"):
        if not ln.strip():
            continue
        a = ln.rstrip("\n").split("\t"); svloc.append(a[0]); tiers.append(a[1] if len(a) > 1 else "")
    meta = list(csv.DictReader(open(f"{MATDIR}/disease_meta.tsv"), delimiter="\t"))
    return Pv, Be, np.array(svloc), np.array(tiers), meta


def _firth_intercept_null(yp, offp):
    d = 0.0
    for _ in range(FIRTH_MAXIT):
        mu = 1.0 / (1.0 + np.exp(-np.clip(offp + d, -30, 30))); w = mu * (1 - mu)
        I = max(float(w.sum()), EPS); h = w / I
        step = float(np.clip(float(np.sum((yp - mu) + h * (0.5 - mu))) / I, -MAXSTEP, MAXSTEP)); d += step
        if abs(step) < TOL:
            break
    mu = 1.0 / (1.0 + np.exp(-np.clip(offp + d, -30, 30))); w = mu * (1 - mu)
    ll = float(np.sum(yp * np.log(np.clip(mu, EPS, 1 - EPS)) + (1 - yp) * np.log(np.clip(1 - mu, EPS, 1 - EPS))))
    return ll, w


def firth_offset(Gb, yp, offp):
    nb = Gb.shape[0]
    delta = np.zeros(nb); beta = np.zeros(nb); done = np.zeros(nb, bool)
    for _ in range(FIRTH_MAXIT):
        eta = np.clip(offp[None, :] + delta[:, None] + beta[:, None] * Gb, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta)); w = mu * (1.0 - mu)
        Sw = w.sum(1); Swg = (w * Gb).sum(1); Swg2 = (w * Gb * Gb).sum(1)
        det = np.maximum(Sw * Swg2 - Swg * Swg, EPS)
        quad = (Swg2[:, None] - 2 * Swg[:, None] * Gb + Sw[:, None] * Gb * Gb) / det[:, None]
        h = w * quad
        ytil = (yp[None, :] - mu) + h * (0.5 - mu)
        Ud = ytil.sum(1); Ub = (ytil * Gb).sum(1)
        dd = np.clip((Swg2 * Ud - Swg * Ub) / det, -MAXSTEP, MAXSTEP)
        db = np.clip((-Swg * Ud + Sw * Ub) / det, -MAXSTEP, MAXSTEP)
        dd[done] = 0.0; db[done] = 0.0
        delta = delta + dd; beta = beta + db
        done = done | ((np.abs(dd) < TOL) & (np.abs(db) < TOL))
        if done.all():
            break
    eta = np.clip(offp[None, :] + delta[:, None] + beta[:, None] * Gb, -30, 30)
    mu = 1.0 / (1.0 + np.exp(-eta)); w = mu * (1.0 - mu)
    Sw = w.sum(1); Swg = (w * Gb).sum(1); Swg2 = (w * Gb * Gb).sum(1)
    det1 = np.maximum(Sw * Swg2 - Swg * Swg, EPS)
    ll1 = (np.log(np.clip(mu, EPS, 1 - EPS)) * yp[None, :]
           + np.log(np.clip(1 - mu, EPS, 1 - EPS)) * (1 - yp[None, :])).sum(1)
    pl1 = ll1 + 0.5 * np.log(det1)
    ll0, w0 = _firth_intercept_null(yp, offp)
    Sw0 = float(w0.sum()); Swg0 = Gb @ w0; Swg2_0 = (Gb * Gb) @ w0
    det0 = np.maximum(Sw0 * Swg2_0 - Swg0 * Swg0, EPS)
    pl0 = ll0 + 0.5 * np.log(det0)
    stat = np.maximum(2.0 * (pl1 - pl0), 0.0)
    return beta, chi2dist.sf(stat, 1), done


def fisher_cell(g, yp_full, present):
    g_p = g[present]; y = yp_full[present]
    carr = g_p >= CARRIER_THRESH
    a = int(np.sum(carr & (y == 1))); b = int(np.sum(carr & (y == 0)))
    c = int(np.sum(~carr & (y == 1))); d = int(np.sum(~carr & (y == 0)))
    try:
        orr, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    except Exception:
        orr, p = float("nan"), float("nan")
    return a, c, float(orr), float(p)


_G = _X = _Y = _Pv = None


def process_phecode(j):
    col = _Pv[:, j]
    band = np.nonzero(np.isfinite(col) & (col < REFIT_CUT))[0]
    if band.size == 0:
        return None
    yb = _Y[:, j]; present = yb >= 0
    y = np.where(present, yb, 0).astype(np.float64)
    phat, _w, _r, _c = S.fit_null(_X, y, present)
    off = np.log(np.clip(phat, 1e-6, 1 - 1e-6) / (1 - np.clip(phat, 1e-6, 1 - 1e-6)))[present]
    yp = y[present]
    beta = np.empty(band.size); pfir = np.empty(band.size); cvg = np.zeros(band.size, bool)
    for s in range(0, band.size, BAND_BATCH):
        sub = band[s:s + BAND_BATCH]
        Gb = _G[sub][:, present].astype(np.float64)
        b, p, c = firth_offset(Gb, yp, off)
        beta[s:s + len(sub)] = b; pfir[s:s + len(sub)] = p; cvg[s:s + len(sub)] = c
    recs = []
    for k, m_i in enumerate(band):
        sep = bool((abs(beta[k]) >= SEP_BETA_BOUND) or (not bool(cvg[k])))
        fa = fc = forr = fp = np.nan
        if sep:
            fa, fc, forr, fp = fisher_cell(_G[m_i], y, present)
        recs.append((int(m_i), float(beta[k]), float(pfir[k]), int(bool(cvg[k])), int(sep), fa, fc, forr, fp))
    return int(j), recs


def main(validate=False):
    t0 = time.time()
    iids, iidx = S.load_order(); N = len(iids)
    X, _ = S.load_X(iidx)
    sv_order = []; G = S.load_G(sv_order); M = len(sv_order)
    phecodes = []; Y = S.load_disease_matrix(iidx, phecodes); J = len(phecodes)
    Pv, Be, svloc_mat, tiers_mat, meta = load_score_matrices()
    assert Pv.shape == (M, J), f"matrix shape {Pv.shape} != loaded G/Y ({M},{J})"
    assert list(svloc_mat) == sv_order, "sv_order mismatch between G and saved matrices — abort"
    ph_label = {m["phecode"]: m.get("label", "") for m in meta}
    ph_cases = {m["phecode"]: int(m["cases"]) for m in meta}
    ph_n = {m["phecode"]: int(m["n"]) for m in meta}
    ph_tested = {m["phecode"]: (m["tested"] == "1") for m in meta}

    global _G, _X, _Y, _Pv
    _G, _X, _Y, _Pv = G, X, Y, Pv
    Pv_f = Pv.copy(); Be_f = Be.copy()
    cells = []
    n_refit = 0; n_sep = 0
    tested_idx = [j for j, ph in enumerate(phecodes) if ph_tested.get(ph, False)]
    print(f"recalibrating: {len(tested_idx)} tested phecodes x {M} SVs, refit band = score-p<{REFIT_CUT} "
          f"(nominal {PTHRESH}; equivalent hit set), {NWORK} workers ...", flush=True)
    from multiprocessing import Pool
    t1 = time.time(); done_ph = 0
    with Pool(NWORK) as pool:
        for res in pool.imap_unordered(process_phecode, tested_idx, chunksize=8):
            done_ph += 1
            if res is not None:
                j, recs = res; ph = phecodes[j]
                for (m_i, beta_k, pfir_k, cvg_k, sep_k, fa, fc, forr, fp) in recs:
                    Pv_f[m_i, j] = pfir_k; Be_f[m_i, j] = beta_k
                    if sep_k:
                        n_sep += 1
                    cells.append((sv_order[m_i], tiers_mat[m_i], ph, ph_label.get(ph, ""),
                                  ph_cases.get(ph, 0), ph_n.get(ph, 0), float(Pv[m_i, j]), pfir_k,
                                  beta_k, float(np.exp(np.clip(beta_k, -50, 50))), cvg_k,
                                  sep_k, fa, fc, forr, fp))
                    n_refit += 1
            if done_ph % 25 == 0:
                print(f"  {done_ph}/{len(tested_idx)} phecodes, {n_refit} cells, {n_sep} sep "
                      f"({time.time()-t1:.0f}s)", flush=True)

    sv_of = np.repeat(np.arange(M), J); ph_of = np.tile(np.arange(J), M)
    tested_arr = np.array([ph_tested.get(p, False) for p in phecodes])
    Pf_flat = Pv_f.ravel(); Bf_flat = Be_f.ravel()
    valid = np.isfinite(Pf_flat) & tested_arr[ph_of]
    q = np.full(M * J, np.nan)
    for tn in ("signal", "null_ext"):
        idx = np.nonzero(valid & (tiers_mat[sv_of] == tn))[0]
        if idx.size:
            q[idx] = S.bh(Pf_flat[idx])
    n_sig = int(np.nansum(q < 0.05)); n_ind = int(np.nansum((q >= 0.05) & (q < 0.10)))

    np.save(f"{OUTDIR}/disease_Pv_firth.npy", Pv_f.astype(np.float32))
    np.save(f"{OUTDIR}/disease_Be_firth.npy", Be_f.astype(np.float32))
    with open(f"{OUTDIR}/disease_firth_cells.tsv", "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tp_score\tp_firth\tbeta_firth\tOR_firth\t"
                "firth_converged\tseparation\tfisher_carrier_cases\tfisher_noncarrier_cases\t"
                "fisher_OR\tfisher_p\n")
        for r in sorted(cells, key=lambda x: x[7]):
            f.write("\t".join("" if (isinstance(v, float) and np.isnan(v)) else
                              (f"{v:.4g}" if isinstance(v, float) else str(v)) for v in r) + "\n")
    sepmap = {(c[0], c[2]): (int(c[11]), c[15]) for c in cells}
    keep = np.nonzero(valid & ((np.nan_to_num(q, nan=1.0) < 0.10) | (Pf_flat < 1e-4)))[0]
    keep = keep[np.argsort(Pf_flat[keep])]
    n_hdr_sig = 0; n_demoted = 0
    with open(f"{OUTDIR}/disease_hits_firthrecal.tsv", "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tbeta_firth_logOR\tp_firth\tq_bh_tier\t"
                "separation\tfisher_p\ttier_label\n")
        for i in keep:
            a = i // J; b = i % J; qi = q[i]
            sepflag, fp = sepmap.get((sv_order[a], phecodes[b]), (0, float("nan")))
            sig = bool(np.isfinite(qi) and qi < 0.05)
            base = ("significant" if sig else "indicative" if (np.isfinite(qi) and qi < 0.10) else "sub-threshold")
            lab = "indicative (separation)" if sepflag else base
            if sig and sepflag:
                n_demoted += 1
            elif sig:
                n_hdr_sig += 1
            f.write(f"{sv_order[a]}\t{tiers_mat[a]}\t{phecodes[b]}\t{ph_label.get(phecodes[b],'')}\t"
                    f"{ph_cases.get(phecodes[b],0)}\t{ph_n.get(phecodes[b],0)}\t{Bf_flat[i]:.4g}\t"
                    f"{Pf_flat[i]:.3g}\t{'' if not np.isfinite(qi) else f'{qi:.4g}'}\t{sepflag}\t"
                    f"{'' if not np.isfinite(fp) else f'{fp:.3g}'}\t{lab}\n")
    print(f"FIRTH_RECAL_DONE: {n_refit} cells refit, {n_sep} separation-flagged; "
          f"headline: {n_hdr_sig} significant (q<0.05, non-separation), {n_demoted} separation-demoted, "
          f"{n_ind} indicative (q0.05-0.1 matrix) -> {OUTDIR} ({time.time()-t0:.0f}s)", flush=True)

    if validate:
        validate_against_exact(cells)


def validate_against_exact(cells):
    ref = f"{S.OUTA}/phewide_disease_full_firthrefit.tsv"
    if not os.path.exists(ref):
        print("VALIDATE: no exact-Firth refit found (run step30h first) — skipping parity check", flush=True)
        return
    exact = {}
    for r in csv.DictReader(open(ref), delimiter="\t"):
        try:
            exact[(r["sv_locus"], r["phecode"])] = (float(r["p_firth"]), int(r["firth_converged"]))
        except (ValueError, KeyError):
            continue
    xs, ys = [], []
    for c in cells:
        key = (c[0], c[2])
        if key in exact and exact[key][1] == 1 and np.isfinite(exact[key][0]) and exact[key][0] > 0:
            xs.append(-np.log10(max(exact[key][0], 1e-300))); ys.append(-np.log10(max(c[7], 1e-300)))
    if len(xs) >= 5:
        r = float(np.corrcoef(xs, ys)[0, 1])
        print(f"VALIDATE: approx-Firth vs exact-Firth on {len(xs)} converged hit cells: "
              f"corr(-log10p)={r:.4f}  (want >~0.95)", flush=True)
    else:
        print(f"VALIDATE: only {len(xs)} comparable converged cells — insufficient for parity", flush=True)


if __name__ == "__main__":
    main(validate="--validate" in sys.argv)
