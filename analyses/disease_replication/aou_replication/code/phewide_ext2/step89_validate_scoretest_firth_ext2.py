#!/usr/bin/env python3
"""Step 89 (ext2 arm) — SCORE-TEST CALIBRATION GATE (<- step30g_validate4_firth.py)."""
import os
import sys
import gc
import csv
import time
import argparse

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import warnings                                       # noqa: E402
import numpy as np                                    # noqa: E402

warnings.filterwarnings("ignore", message=".*Maximum number of iterations.*")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import step51_disease_scoretest_ext2 as S             # noqa: E402
from frozen_tree_sig import tree_sig                  # noqa: E402
from firthmodels import FirthLogisticRegression       # noqa: E402

HOME = os.path.expanduser("~")
RV = (f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/"
      f"aou_replication/results_vds")
OUT = f"{S.OUTA}/gates"
CELLS = f"{OUT}/validation_firth_v4_ext2.tsv"
REPORT = f"{OUT}/validation_firth_v4_ext2_report.tsv"
TIERF = f"{RV}/step28_phewide_ext2/sv_tier.tsv"
FROZEN = ["step30_phewide", "step28_phewide", "phewide_carryover", "step30_phewide_mcc1"]

FIRTH_MAXIT = 25
FIRTH_PEN = 0.5
CARRIER = 0.5
CASE_TARGETS = (60, 120, 300, 800, 2000, 8000, 30000, 90000)
HEADLINE = ("chr16_sv1029", "chr16_sv37", "chr16_sv73", "chr16_sv465",
            "chr10_sv697", "chr10_sv698", "chr11_sv200", "chr11_sv400")
CARRIER_BINS = ((0, 5), (5, 20), (20, 100), (100, 10 ** 9))

_CTX = {}


def say(m=""):
    print(m, flush=True)


def read_tier(path):
    t = {}
    with open(path) as f:
        f.readline()
        for line in f:
            q = line.rstrip("\n").split("\t")
            t[q[0]] = q[1]
    return t


def pick_phenos(colnames):
    rows = [r for r in csv.DictReader(open(S.PSUM), delimiter="\t") if r["phecode"] in set(colnames)]
    by = {r["phecode"]: (int(r["cases_mcc2"]), r.get("sex", "Both")) for r in rows}
    have = [(ph, by[ph][0]) for ph in colnames if ph in by]
    if not have:
        raise SystemExit("no phecode in disease.pheno matched the phenotype summary")
    have.sort(key=lambda x: x[1])
    chosen = []
    for target in CASE_TARGETS:
        cand = min(have, key=lambda x: abs(x[1] - target))
        if cand[0] not in [c[0] for c in chosen]:
            chosen.append(cand)
    for ph, (nc, sx) in by.items():
        if sx != "Both" and ph in colnames and ph not in [c[0] for c in chosen] and nc >= 100:
            chosen.append((ph, nc))
            break
    return [c[0] for c in chosen]


def build_grid(sv_order, carr_frac, tiers, n_per_tier):
    by_tier = {}
    for i, s in enumerate(sv_order):
        by_tier.setdefault(tiers.get(s, "UNTIERED"), []).append(i)
    idx = set()
    for t, ii in sorted(by_tier.items()):
        ii = np.array(ii)
        order = ii[np.argsort(carr_frac[ii])]
        k = min(n_per_tier, len(order))
        if k:
            idx |= set(order[np.linspace(0, len(order) - 1, k).astype(int)].tolist())
        say(f"  tier {t:<10} {len(ii):>6} SVs -> {k} sampled; "
            f"carrier-frac p0/p50/p100 = {np.round(np.percentile(carr_frac[ii], [0, 50, 100]), 5)}")
    add = [i for i, s in enumerate(sv_order) if any(s.startswith(h + "_") for h in HEADLINE)]
    idx |= set(add)
    say(f"  + {len(add)} headline SVs appended -> grid {len(idx)} SVs")
    return sorted(idx)


def self_checks(G, X):
    ok = True
    rng = np.random.default_rng(0)
    N = G.shape[1]
    y = (rng.random(N) < 0.05).astype(np.int8)
    g = G[0].astype(np.float32)

    def one(gv):
        Be, Pv, *_ = S.score_test_batch(gv[None], (gv * gv)[None], X, y[:, None].astype(np.int8), min_case=5)
        return float(Pv[0, 0]), float(Be[0, 0])

    p1, b1 = one(g)
    g2 = (2.0 * (g - g.min()) / (np.ptp(g) + 1e-12)).astype(np.float32)
    p2, b2 = one(g2)
    inv = np.isclose(p1, p2, rtol=1e-4, equal_nan=True) and np.sign(b1) == np.sign(b2)
    say(f"[self-check 1] affine [0,2] invariance : p {p1:.6g} vs {p2:.6g}, sign {np.sign(b1):+.0f}/"
        f"{np.sign(b2):+.0f}  -> {'PASS' if inv else 'FAIL'}")
    ok &= bool(inv)

    _, w, _, _ = S.fit_null(X, y.astype(np.float64), np.ones(N, bool))
    Gk = G[:32].astype(np.float32)

    def corr_of(Xu):
        WX = Xu * w[:, None]
        A = (Gk @ WX).astype(np.float64)
        return np.einsum("mk,kl,ml->m", A, np.linalg.pinv(Xu.T @ WX), A)

    c1 = corr_of(X)
    c2 = corr_of(np.column_stack([X, X[:, 1]]))
    rel = float(np.max(np.abs(c2 - c1) / np.maximum(np.abs(c1), 1e-300)))
    same = rel < 1e-6
    say(f"[self-check 2] pinv == column-drop     : max |rel diff| in the score correction = {rel:.3e}"
        f"  -> {'PASS' if same else 'FAIL'}")
    ok &= bool(same)

    Xd = np.column_stack([X, X[:, 1]])
    _, _, rr, conv = S.fit_null(Xd, y.astype(np.float64), np.ones(N, bool))
    _, Pd, *_ = S.score_test_batch(g[None], (g * g)[None], Xd, y[:, None].astype(np.int8), min_case=5)
    trips = abs(float(rr.sum())) > S.IDENT_TOL
    say(f"[observation ] rank-deficient design   : null IRLS converged={conv}, "
        f"|sum r|={abs(float(rr.sum())):.4g} vs IDENT_TOL={S.IDENT_TOL} -> degeneracy guard "
        f"{'TRIPS' if trips else 'does not trip'}; engine returns p={float(Pd[0,0]):.6g}")
    if trips and not np.isfinite(float(Pd[0, 0])):
        say("               ⟹ FAIL-CLOSED: the engine refuses the cell rather than emitting a number.")
        say("                 So the pinv is defensive coding for a case the |sum r| guard already")
        say("                 catches at production scale (it does NOT trip on a small synthetic design,")
        say("                 so this is scale/conditioning dependent, not universal).")
    elif np.isfinite(float(Pd[0, 0])):
        say("               ⟹ the engine emitted a p on a rank-deficient design; pinv carried it (see "
            "self-check 2 for whether that value is the column-drop value).")
    if not ok:
        raise SystemExit("*** SELF-CHECK FAILED — the engine's asserted invariances do not hold; "
                         "abort before generating any calibration number")
    return ok


def firth_fit(Xcov, g, y, present, max_iter=FIRTH_MAXIT):
    Xp = Xcov[present]
    gp = g[present]
    yp = y[present].astype(int)
    keep = np.array([np.ptp(Xp[:, j]) > 0 for j in range(Xp.shape[1])])
    Xf = np.column_stack([Xp[:, keep], gp])
    var_index = int(keep.sum())
    m = FirthLogisticRegression(fit_intercept=True, max_iter=max_iter, penalty_weight=FIRTH_PEN)
    try:
        m.fit(Xf, yp)
    except Exception as e:
        msg = str(e).lower()
        return None, None, ("separation" if ("separat" in msg or "singular" in msg or "converg" in msg)
                            else "fail")
    beta = float(m.coef_[var_index])
    try:
        m.lrt(var_index)
        p = float(m.lrt_pvalues_[var_index])
    except Exception:
        p = float(m.pvalues_[var_index])
    status = "ok" if getattr(m, "converged_", True) else "separation"
    return beta, p, status


def _work(job):
    a, ph = job
    Gs = _CTX["Gs"]
    Xcov = _CTX["Xcov"]
    y = _CTX["Y"][ph]
    present = y >= 0
    case_mask = y == 1
    carr_cases = int((Gs[a][case_mask] > CARRIER).sum())
    carr_all = int((Gs[a] > CARRIER).sum())
    bf, pf, status = firth_fit(Xcov, Gs[a], y, present)
    bf2 = pf2 = None
    st2 = ""
    if status != "ok" and _CTX["retry_iter"]:
        bf2, pf2, st2 = firth_fit(Xcov, Gs[a], y, present, max_iter=_CTX["retry_iter"])
    return dict(a=a, phecode=ph, carrier_cases=carr_cases, carrier_all=carr_all,
                p_firth=pf, beta_firth=bf, firth_status=status,
                p_firth_retry=pf2, beta_firth_retry=bf2, firth_status_retry=st2)


def fmt(v, spec):
    return "" if v is None or not np.isfinite(v) else format(v, spec)


def report(rows, out_rows, label, sel):
    s = [r for r in rows if sel(r)]
    both = [r for r in s if r["firth_status"] == "ok" and r["p_score"] > 0
            and r["p_firth"] is not None and np.isfinite(r["p_firth"]) and r["p_firth"] > 0]
    nsep = sum(1 for r in s if r["firth_status"] != "ok")
    rec = dict(group=label, n_cells=len(s), n_firth_ok=len(both), n_firth_nonfit=nsep)
    if len(both) >= 4:
        lr = np.array([np.log10(r["p_score"] / r["p_firth"]) for r in both])
        nd = sum(1 for r in both if np.sign(r["beta_score"]) == np.sign(r["beta_firth"]))
        antic = sum(1 for r in both if r["p_score"] < 0.5 * r["p_firth"])
        rec.update(dir_concordant_pct=f"{100*nd/len(both):.1f}",
                   median_log10_ps_over_pf=f"{np.median(lr):+.3f}",
                   frac_ps_lt_half_pf_pct=f"{100*antic/len(both):.1f}",
                   frac_null_pscore_gt_005_pct=f"{100*sum(1 for r in s if r['p_score']>0.05)/max(len(s),1):.1f}")
        say(f"  {label:<34} n={len(s):>5}  ok={len(both):>5}  nonfit={nsep:>4}  "
            f"dir={100*nd/len(both):5.1f}%  med log10(ps/pf)={np.median(lr):+.3f}  "
            f"frac(ps<0.5pf)={100*antic/len(both):5.1f}%")
    else:
        say(f"  {label:<34} n={len(s):>5}  ok={len(both):>5}  nonfit={nsep:>4}  (too few for a ratio)")
    out_rows.append(rec)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-tier", type=int, default=12)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="smoke test: cap the number of cells")
    ap.add_argument("--resume", action="store_true", help="skip cells already in the output file")
    ap.add_argument("--retry-iter", type=int, default=400,
                    help="re-fit non-converged Firth cells at this max_iter (0 = off) to separate genuine "
                         "quasi-separation from an exhausted iteration budget")
    a = ap.parse_args()
    t0 = time.time()
    say("=" * 100)
    say("step89 — score-test vs exact-Firth CALIBRATION GATE on ext2 predictors")
    say("=" * 100)
    if os.path.islink(OUT):
        raise SystemExit(f"{OUT} is a symlink — refusing to write through it")
    os.makedirs(OUT, exist_ok=True)
    base = {n: tree_sig(f"{RV}/{n}") for n in FROZEN}
    for n, (c, _m, k) in base.items():
        say(f"[V-pre] {n:<22} n={k:<5} content={c}")

    tiers = read_tier(TIERF)
    iids, iidx = S.load_order()
    N = len(iids)
    X, _ = S.load_X(iidx)
    Xcov = X[:, 1:]
    sv_order = []
    G = S.load_G(sv_order)
    M = len(sv_order)
    miss = [s for s in sv_order if s not in tiers]
    if miss:
        raise SystemExit(f"{len(miss)} SVs have no tier (fail-open risk), e.g. {miss[:3]}")

    carr_frac = (G > CARRIER).mean(axis=1)
    say("\nGRID (selection is by carrier frequency and tier ONLY — never by significance):")
    sv_idx = build_grid(sv_order, carr_frac, tiers, a.n_per_tier)
    sub_ids = [sv_order[i] for i in sv_idx]
    Gs = np.ascontiguousarray(G[sv_idx])
    G2s = Gs * Gs
    self_checks(G, X)
    del G
    gc.collect()
    say(f"  freed the full dosage matrix; grid array {Gs.nbytes/1e6:.0f} MB")

    hdr = open(f"{S.INP}/disease.pheno").readline().rstrip("\n").split("\t")
    phenos = pick_phenos(hdr[2:])
    colof = {ph: hdr.index(ph) for ph in phenos if ph in hdr}
    say(f"\nPHENOTYPES ({len(colof)}): {list(colof)}")
    Y = {ph: np.full(N, -1, np.int8) for ph in colof}
    with open(f"{S.INP}/disease.pheno") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            kk = iidx.get(p[1])
            if kk is None:
                continue
            for ph, c in colof.items():
                v = p[c]
                if v == "0":
                    Y[ph][kk] = 0
                elif v == "1":
                    Y[ph][kk] = 1

    score = {}
    for ph in colof:
        yb = Y[ph][:, None]
        Be, Pv, ncase, npres, tested, nonconv, degen = S.score_test_batch(Gs, G2s, X, yb, min_case=50)
        score[ph] = (Be[:, 0], Pv[:, 0], int(ncase[0]), bool(tested[0]))
        say(f"  score test {ph:<12} cases={int(ncase[0]):>7}  tested={bool(tested[0])}  "
            f"degenerate={bool(degen[0])}  null-nonconv={bool(nonconv[0])}")

    jobs = [(i, ph) for ph in colof if score[ph][3] for i in range(len(sv_idx))]
    if a.limit:
        jobs = jobs[:a.limit]
    done = set()
    if a.resume and os.path.exists(CELLS):
        with open(CELLS) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                done.add((r["sv_locus"], r["phecode"]))
        jobs = [(i, ph) for i, ph in jobs if (sub_ids[i], ph) not in done]
        say(f"\n--resume: {len(done)} cells already present, {len(jobs)} to run")
    say(f"\nFIRTH REFITS: {len(jobs)} cells over {a.workers} workers "
        f"({len(sv_idx)} SVs x {sum(1 for ph in colof if score[ph][3])} testable phenotypes)")

    _CTX["Gs"] = Gs
    _CTX["Xcov"] = Xcov
    _CTX["Y"] = Y
    _CTX["retry_iter"] = a.retry_iter
    HEAD = ("sv_locus\ttier\tphecode\tcases\tcarrier_cases\tcarrier_all\tcarrier_frac\t"
            "p_score\tbeta_score\tp_firth\tbeta_firth\tfirth_status\t"
            "p_firth_retry\tbeta_firth_retry\tfirth_status_retry\n")
    new_file = not (a.resume and os.path.exists(CELLS))
    fh = open(CELLS, "w" if new_file else "a")
    if new_file:
        fh.write(HEAD)
        fh.flush()

    def emit(r):
        i = r["a"]
        ph = r["phecode"]
        bs, ps, nca, _ = score[ph]
        row = (f"{sub_ids[i]}\t{tiers.get(sub_ids[i],'')}\t{ph}\t{nca}\t{r['carrier_cases']}\t"
               f"{r['carrier_all']}\t{carr_frac[sv_idx[i]]:.6g}\t{fmt(float(ps[i]),'.3g')}\t"
               f"{fmt(float(bs[i]),'.4g')}\t{fmt(r['p_firth'],'.3g')}\t{fmt(r['beta_firth'],'.4g')}\t"
               f"{r['firth_status']}\t{fmt(r.get('p_firth_retry'),'.3g')}\t"
               f"{fmt(r.get('beta_firth_retry'),'.4g')}\t{r.get('firth_status_retry','')}\n")
        fh.write(row)
        fh.flush()

    t1 = time.time()
    ndone = 0
    if a.workers > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(a.workers) as pool:
            for r in pool.imap_unordered(_work, jobs, chunksize=1):
                emit(r)
                ndone += 1
                if ndone % 25 == 0:
                    el = time.time() - t1
                    say(f"    {ndone}/{len(jobs)} cells ({el:.0f}s, {el/ndone:.2f} s/cell, "
                        f"eta {(len(jobs)-ndone)*el/ndone/60:.1f} min)")
    else:
        for j in jobs:
            emit(_work(j))
            ndone += 1
    fh.close()
    say(f"  {ndone} cells in {time.time()-t1:.0f}s ({(time.time()-t1)/max(ndone,1):.2f} s/cell)")

    rows = []
    with open(CELLS) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            def fl(k):
                try:
                    return float(r[k])
                except (ValueError, KeyError):
                    return float("nan")
            rows.append(dict(sv_locus=r["sv_locus"], tier=r["tier"], phecode=r["phecode"],
                             cases=int(r["cases"]), carrier_cases=int(r["carrier_cases"]),
                             p_score=fl("p_score"), beta_score=fl("beta_score"),
                             p_firth=fl("p_firth"), beta_firth=fl("beta_firth"),
                             firth_status=r["firth_status"],
                             p_firth_retry=fl("p_firth_retry"), beta_firth_retry=fl("beta_firth_retry"),
                             firth_status_retry=r.get("firth_status_retry", "")))
    rows = [r for r in rows if np.isfinite(r["p_score"])]

    out_rows = []
    say("\n" + "-" * 100)
    say(f"CALIBRATION REPORT — {len(rows)} cells with a score-test p "
        f"(log10(ps/pf) < 0 ⟹ the score test is ANTI-CONSERVATIVE)")
    say(f"  null support: {100*sum(1 for r in rows if r['p_score']>0.05)/max(len(rows),1):.1f}% of cells at "
        f"p_score>0.05, {100*sum(1 for r in rows if r['p_score']>0.5)/max(len(rows),1):.1f}% at p_score>0.5 "
        f"(published gate: 91.8% / 55.4%)")
    say("\n-- by CARRIER COUNT AMONG CASES --")
    for lo, hi in CARRIER_BINS:
        report(rows, out_rows, f"carrier_cases [{lo},{hi})",
               lambda r, lo=lo, hi=hi: lo <= r["carrier_cases"] < hi)
    say("\n-- by TIER (the reason this step exists: ext2_new and ext2_chrX were never covered) --")
    for t in sorted({r["tier"] for r in rows}):
        report(rows, out_rows, f"tier {t}", lambda r, t=t: r["tier"] == t)
    say("\n-- by CASE COUNT --")
    for lo, hi in ((0, 200), (200, 2000), (2000, 20000), (20000, 10 ** 9)):
        report(rows, out_rows, f"cases [{lo},{hi})", lambda r, lo=lo, hi=hi: lo <= r["cases"] < hi)
    say("\n-- ALL --")
    allrec = report(rows, out_rows, "ALL", lambda r: True)

    nf = [r for r in rows if r["firth_status"] != "ok"]
    if nf:
        say("\n-- NON-CONVERGED FIRTH: genuine quasi-separation, or just an exhausted iteration budget? --")
        byph = {}
        for r in nf:
            byph.setdefault(r["phecode"], []).append(r)
        say(f"  {len(nf)}/{len(rows)} cells non-converged at max_iter={FIRTH_MAXIT}, "
            f"concentrated in: " + ", ".join(f"{k} ({len(v)})" for k, v in sorted(byph.items())))
        rec = sum(1 for r in nf if r.get("firth_status_retry") == "ok")
        if any(r.get("firth_status_retry") for r in nf):
            say(f"  re-fit at the larger budget: {rec}/{len(nf)} converge ⟹ "
                + ("they were ITERATION-LIMITED, not separated. `firth_status` overstates separation, in "
                   "this gate and in the published one." if rec > len(nf) / 2 else
                   "most are genuinely separated, so `firth_status` is reporting real quasi-separation."))
            agree = [r for r in nf if r.get("firth_status_retry") == "ok"
                     and np.isfinite(r.get("p_firth_retry", np.nan)) and np.isfinite(r["p_firth"])
                     and r["p_firth"] > 0 and r["p_firth_retry"] > 0]
            if agree:
                d = np.abs(np.log10([r["p_firth"] / r["p_firth_retry"] for r in agree]))
                say(f"  where both are available, |log10(p25/p_retry)| median {np.median(d):.3f}, "
                    f"max {np.max(d):.3f} over {len(agree)} cells")
            out_rows.append(dict(group="firth nonconv -> converged at retry", n_cells=len(nf),
                                 n_firth_ok=rec, n_firth_nonfit=len(nf) - rec))
        merged = []
        for r in rows:
            rr = dict(r)
            if rr["firth_status"] != "ok" and rr.get("firth_status_retry") == "ok":
                rr["p_firth"] = rr["p_firth_retry"]
                rr["beta_firth"] = rr["beta_firth_retry"]
                rr["firth_status"] = "ok"
            merged.append(rr)
        say("\n-- ALL, using the retry fit wherever it converged (the fuller-coverage verdict) --")
        report(merged, out_rows, "ALL (retry-merged)", lambda r: True)
        for t in sorted({r["tier"] for r in merged}):
            report(merged, out_rows, f"tier {t} (retry-merged)", lambda r, t=t: r["tier"] == t)

    say("\n-- SCREEN PRECISION: among score-test positives, does exact Firth confirm? --")
    both = [r for r in rows if r["firth_status"] == "ok" and r["p_score"] > 0
            and np.isfinite(r["p_firth"]) and r["p_firth"] > 0]
    for thr in (1e-3, 1e-4, 1e-5):
        pos = [r for r in both if r["p_score"] < thr]
        if not pos:
            say(f"  score p<{thr:.0e}: 0 positives")
            out_rows.append(dict(group=f"precision p<{thr:.0e}", n_cells=0))
            continue
        conf = sum(1 for r in pos if r["p_firth"] < 0.05)
        cd = sum(1 for r in pos if r["p_firth"] < 0.05 and np.sign(r["beta_score"]) == np.sign(r["beta_firth"]))
        say(f"  score p<{thr:.0e}: {len(pos):>4} positives, {100*conf/len(pos):.1f}% Firth-confirmed "
            f"(p<0.05), {100*cd/len(pos):.1f}% confirmed AND concordant")
        out_rows.append(dict(group=f"precision p<{thr:.0e}", n_cells=len(pos),
                             dir_concordant_pct=f"{100*cd/len(pos):.1f}",
                             frac_ps_lt_half_pf_pct=f"{100*conf/len(pos):.1f}"))

    hdr2 = ["group", "n_cells", "n_firth_ok", "n_firth_nonfit", "dir_concordant_pct",
            "median_log10_ps_over_pf", "frac_ps_lt_half_pf_pct", "frac_null_pscore_gt_005_pct"]
    with open(REPORT, "w") as f:
        f.write("\t".join(hdr2) + "\n")
        for r in out_rows:
            f.write("\t".join(str(r.get(h, "")) for h in hdr2) + "\n")
    say(f"\nwrote {CELLS}\nwrote {REPORT}")

    changed = [n for n in FROZEN if tree_sig(f"{RV}/{n}")[0] != base[n][0]]
    if changed:
        raise SystemExit(f"*** frozen tree(s) modified: {changed}")
    say("[V-post] all 4 frozen trees byte-identical — PASS")
    say(f"\nstep89 done in {time.time()-t0:.0f}s")
    say("VALIDATION_EXT2_DONE")


if __name__ == "__main__":
    main()
