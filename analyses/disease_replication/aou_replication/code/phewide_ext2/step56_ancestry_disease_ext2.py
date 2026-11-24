#!/usr/bin/env python3
"""Step 56 (ext2 arm) — PER-ANCESTRY disease score test over all 10,799 ext2 predictors."""
import os, sys, csv, time, gc, argparse, resource
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8"); os.environ.setdefault("OMP_NUM_THREADS", "8")
import numpy as np
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S

ANC_DUMMIES = ["anc_amr", "anc_eas", "anc_eur", "anc_mid", "anc_sas"]
DUMMY_TO_ANC = {"anc_amr": "AMR", "anc_eas": "EAS", "anc_eur": "EUR", "anc_mid": "MID", "anc_sas": "SAS"}
ANCS = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
WELL = {"EUR", "AFR", "AMR"}
ORDER = ["MID", "SAS", "EAS", "AMR", "AFR", "EUR"]
MATDIR = f"{S.MATDIR}/per_ancestry"
SVTIER = f"{S.SRC}/sv_tier.tsv"
BATCH = 64
FLUSH_EVERY = 6
MEM_FRAC = 0.85



def score_test_batch_uv(G, G2, X, Yb, min_case=50):
    N, B = Yb.shape; M = G.shape[0]; K = X.shape[1]
    RHS = np.zeros((N, B + B * K), np.float32)
    Wm = np.zeros((N, B), np.float32)
    XtWXi = [None] * B
    ncase = np.zeros(B, np.int64); npres = np.zeros(B, np.int64); tested = np.zeros(B, bool)
    nonconv = np.zeros(B, bool); degen = np.zeros(B, bool)
    for b in range(B):
        yb = Yb[:, b]; present = yb >= 0
        nca = int((yb == 1).sum()); npr = int(present.sum())
        ncase[b] = nca; npres[b] = npr
        if nca < min_case or (npr - nca) < min_case:
            continue
        y = np.where(present, yb, 0).astype(np.float64)
        _, w, r, conv = S.fit_null(X, y, present)
        nonconv[b] = not conv
        if abs(float(r.sum())) > S.IDENT_TOL:
            degen[b] = True
            continue
        WX = X * w[:, None]
        XtWXi[b] = np.linalg.pinv(X.T @ WX)
        RHS[:, b] = r
        RHS[:, B + b * K:B + (b + 1) * K] = WX
        Wm[:, b] = w
        tested[b] = True
    GX = (G @ RHS).astype(np.float64)
    GU = GX[:, :B]; GA = GX[:, B:]
    Sww = (G2 @ Wm).astype(np.float64)
    Be = np.zeros((M, B), np.float32); Pv = np.full((M, B), np.nan, np.float32)
    Umat = np.full((M, B), np.nan, np.float32); Vmat = np.full((M, B), np.nan, np.float32)
    for b in range(B):
        if not tested[b]:
            continue
        A = GA[:, b * K:(b + 1) * K]
        corr = np.einsum("mk,kl,ml->m", A, XtWXi[b], A)
        Vraw = Sww[:, b]
        V = Vraw - corr
        U = GU[:, b]
        with np.errstate(divide="ignore", invalid="ignore"):
            stat = (U * U) / V; beta = U / V
        good = np.isfinite(stat) & (V > S.VREL * Vraw) & (Vraw > 0)
        stat = np.where(good, stat, np.nan)
        Pv[:, b] = chi2dist.sf(stat, 1); Be[:, b] = np.where(good, beta, np.nan).astype(np.float32)
        Umat[:, b] = np.where(good, U, np.nan).astype(np.float32)
        Vmat[:, b] = np.where(good, V, np.nan).astype(np.float32)
    return Be, Pv, ncase, npres, tested, nonconv, degen, Umat, Vmat


def save_matrices_tag(Pv, Be, phecodes, label, ncase, npres, tested, sv_order, tiers, out_tag, progress):
    def _save(path, arr):
        np.save(path + ".tmp.npy", arr); os.replace(path + ".tmp.npy", path)
    _save(f"{MATDIR}/disease_Pv{out_tag}.npy", Pv)
    _save(f"{MATDIR}/disease_Be{out_tag}.npy", Be)
    with open(f"{MATDIR}/disease_meta{out_tag}.tsv.tmp", "w") as f:
        f.write("phecode\tlabel\tcases\tn\ttested\n")
        for b in range(len(phecodes)):
            f.write(f"{phecodes[b]}\t{label.get(phecodes[b],'')}\t{ncase[b]}\t{npres[b]}\t{int(tested[b])}\n")
    os.replace(f"{MATDIR}/disease_meta{out_tag}.tsv.tmp", f"{MATDIR}/disease_meta{out_tag}.tsv")
    with open(f"{MATDIR}/disease_svorder{out_tag}.txt.tmp", "w") as f:
        f.write("\n".join(f"{sv_order[a]}\t{tiers[a]}" for a in range(len(sv_order))) + "\n")
    os.replace(f"{MATDIR}/disease_svorder{out_tag}.txt.tmp", f"{MATDIR}/disease_svorder{out_tag}.txt")
    with open(f"{MATDIR}/disease_progress{out_tag}.txt", "w") as f:
        f.write(progress + "\n")


def parity_gate(G, X, Y, keep_cols, n_sv=128, n_ph=24, n_samp=40000):
    t0 = time.time()
    Gs = np.ascontiguousarray(G[:n_sv, :n_samp]); G2s = Gs * Gs
    Xs = np.ascontiguousarray(X[:n_samp][:, keep_cols])
    Ys_all = Y[:n_samp]
    order = np.argsort(-(Ys_all == 1).sum(0))[:n_ph]
    Ys = np.ascontiguousarray(Ys_all[:, np.sort(order)])
    ref = S.score_test_batch(Gs, G2s, Xs, Ys, min_case=50)
    got = score_test_batch_uv(Gs, G2s, Xs, Ys, min_case=50)
    names = ["Be", "Pv", "ncase", "npres", "tested", "nonconv", "degen"]
    for nm, a, b in zip(names, ref, got[:7]):
        a = np.asarray(a); b = np.asarray(b)
        if a.dtype.kind == "f":
            ok = np.array_equal(a, b, equal_nan=True)
        else:
            ok = np.array_equal(a, b)
        if not ok:
            d = np.nanmax(np.abs(a.astype(float) - b.astype(float)))
            raise SystemExit(f"ABORT parity_gate: '{nm}' differs between step51.score_test_batch and "
                             f"score_test_batch_uv (max|delta|={d:.3e}). The transcription drifted.")
    U, V = got[7], got[8]
    tested = got[4]
    fin = np.isfinite(U) & np.isfinite(V)
    print(f"PARITY GATE PASS: 7/7 return values bitwise identical on a real "
          f"{n_sv}sv x {n_ph}ph x {n_samp}samp slice ({int(tested.sum())} cols tested); "
          f"U,V finite in {fin.sum()}/{U.size} cells ({time.time()-t0:.1f}s)", flush=True)
    del Gs, G2s, Xs, Ys; gc.collect()



def ancestry_labels(X, names, N):
    lab = np.array(["AFR"] * N, dtype=object)
    for d in ANC_DUMMIES:
        col = names.index(d)
        lab[X[:, col] == 1.0] = DUMMY_TO_ANC[d]
    return lab


def load_tiers(sv_order):
    if not os.path.exists(SVTIER):
        raise SystemExit(f"ABORT: ext2 sv_tier.tsv not found at {SVTIER}")
    tmap = {r["sv_locus"]: r["tier"] for r in csv.DictReader(open(SVTIER), delimiter="\t")}
    missing = [s for s in sv_order if not tmap.get(s)]
    if missing:
        raise SystemExit(f"ABORT: {len(missing)} scored SVs have no tier (e.g. {missing[:5]}) — "
                         f"the fail-open would drop them from FDR entirely.")
    tiers = np.array([tmap[s] for s in sv_order])
    from collections import Counter
    print(f"tiers (data-derived, fail-closed): {dict(Counter(tiers.tolist()))}", flush=True)
    return tiers


def mem_available_gb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / 1048576.0
    return 0.0


def peak_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0


def cur_rss_gb():
    for line in open("/proc/self/status"):
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) / 1048576.0
    return 0.0


def mem_preflight(M, N, J, K, n_max, batch):
    gb = lambda b: b / 2 ** 30
    parts = {
        "G [MxN] f32": gb(M * N * 4), "Y [NxJ] i8": gb(N * J),
        "Ga f32": gb(M * n_max * 4), "G2a f32": gb(M * n_max * 4), "Ya i8": gb(n_max * J),
        "RHS f32": gb(n_max * (batch + batch * K) * 4), "Pv/Be/U/V f32": gb(4 * M * J * 4),
    }
    proj = sum(parts.values())
    avail = mem_available_gb(); rss = cur_rss_gb(); headroom = avail + rss
    print("--- memory preflight (largest stratum n=%d, K=%d, batch=%d) ---" % (n_max, K, batch), flush=True)
    for k, v in parts.items():
        print(f"    {k:<16} {v:6.2f} GB", flush=True)
    print(f"    {'PROJECTED PEAK':<16} {proj:6.2f} GB   vs headroom {headroom:.2f} GB "
          f"(MemAvailable {avail:.2f} + own RSS {rss:.2f}); limit {MEM_FRAC:.0%} = {MEM_FRAC*headroom:.2f} GB",
          flush=True)
    if proj > MEM_FRAC * headroom:
        raise SystemExit(f"ABORT: projected peak {proj:.1f} GB exceeds {MEM_FRAC:.0%} of headroom "
                         f"({headroom:.1f} GB = MemAvailable {avail:.1f} + own RSS {rss:.1f}). "
                         f"Run strata separately (--only) or split the SV set.")
    print(f"    -> fits monolithically; no step40 batch pooling required.\n", flush=True)
    return proj



def scan_subset(Ga, G2a, Xa, Ya, phecodes, tiers, label, sv_order, out_tag, min_case, t0):
    M = Ga.shape[0]; J = len(phecodes); na = Ga.shape[1]
    Pv = np.full((M, J), np.nan, np.float32); Be = np.full((M, J), np.nan, np.float32)
    Um = np.full((M, J), np.nan, np.float32); Vm = np.full((M, J), np.nan, np.float32)
    ncase = np.zeros(J, np.int32); npres = np.zeros(J, np.int32); tested = np.zeros(J, bool)
    n_nonconv = 0; n_degen = 0

    def _save_npy(path, arr):
        np.save(path + ".tmp.npy", arr); os.replace(path + ".tmp.npy", path)

    pv_path = f"{MATDIR}/disease_Pv{out_tag}.npy"; ord_path = f"{MATDIR}/disease_svorder{out_tag}.txt"
    meta_path = f"{MATDIR}/disease_meta{out_tag}.tsv"; u_path = f"{MATDIR}/disease_U{out_tag}.npy"
    if all(os.path.exists(p) for p in (pv_path, ord_path, meta_path, u_path)):
        prev_order = [ln.split("\t")[0] for ln in open(ord_path) if ln.strip()]
        prevPv = np.load(pv_path)
        if prev_order == sv_order and prevPv.shape == (M, J):
            Pv[:] = prevPv; Be[:] = np.load(f"{MATDIR}/disease_Be{out_tag}.npy")
            Um[:] = np.load(u_path); Vm[:] = np.load(f"{MATDIR}/disease_V{out_tag}.npy")
            for b, r in enumerate(csv.DictReader(open(meta_path), delimiter="\t")):
                tested[b] = r["tested"] == "1"; ncase[b] = int(r["cases"]); npres[b] = int(r["n"])
            print(f"    [{out_tag}] RESUME: {int(tested.sum())}/{J} cols already done — skipping", flush=True)

    def flush(progress):
        save_matrices_tag(Pv, Be, phecodes, label, ncase, npres, tested, sv_order, tiers, out_tag, progress)
        _save_npy(f"{MATDIR}/disease_U{out_tag}.npy", Um)
        _save_npy(f"{MATDIR}/disease_V{out_tag}.npy", Vm)

    nbatch = 0
    for j0 in range(0, J, BATCH):
        cols = [j for j in range(j0, min(j0 + BATCH, J)) if not tested[j]]
        if not cols:
            continue
        Be_b, Pv_b, nca_b, npr_b, tt_b, nc_b, dg_b, U_b, V_b = score_test_batch_uv(
            Ga, G2a, Xa, Ya[:, cols], min_case=min_case)
        n_nonconv += int(nc_b.sum()); n_degen += int(dg_b.sum())
        for a, j in enumerate(cols):
            ncase[j] = nca_b[a]; npres[j] = npr_b[a]
            if tt_b[a]:
                Be[:, j] = Be_b[:, a]; Pv[:, j] = Pv_b[:, a]
                Um[:, j] = U_b[:, a]; Vm[:, j] = V_b[:, a]; tested[j] = True
        nbatch += 1
        if nbatch % 10 == 0:
            print(f"    [{out_tag}] {min(j0+BATCH,J)}/{J} scanned, {int(tested.sum())} tested, "
                  f"rss {peak_rss_gb():.1f} GB ({time.time()-t0:.0f}s)", flush=True)
        if nbatch % FLUSH_EVERY == 0:
            flush(f"in_progress {out_tag}: {min(j0+BATCH,J)}/{J} cols ({time.time()-t0:.0f}s)")

    sv_of = np.repeat(np.arange(M), J); ph_of = np.tile(np.arange(J), M)
    Pv_flat = Pv.ravel(); Be_flat = Be.ravel()
    valid = np.isfinite(Pv_flat) & tested[ph_of]
    q = np.full(M * J, np.nan)
    tier_names = sorted({t for t in tiers.tolist() if t})
    for tn in tier_names:
        idx = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if idx.size:
            q[idx] = S.bh(Pv_flat[idx])
    n_q05 = int(np.nansum(q < 0.05))
    keep = np.nonzero(valid & ((np.nan_to_num(q, nan=1.0) < 0.05) | (Pv_flat < 1e-4)))[0]
    keep = keep[np.argsort(Pv_flat[keep])]
    out_c = f"{S.OUTA}/per_ancestry/phewide_disease_full{out_tag}.tsv"
    lab = {r["phecode"]: r["label"] for r in csv.DictReader(open(S.PSUM), delimiter="\t")}
    with open(out_c, "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tbeta_logOR\tchi2_p\tq_bh_tier\n")
        for i in keep:
            aa = i // J; bb = i % J
            f.write(f"{sv_order[aa]}\t{tiers[aa]}\t{phecodes[bb]}\t{lab.get(phecodes[bb],'')}\t{ncase[bb]}\t"
                    f"{npres[bb]}\t{Be_flat[i]:.4g}\t{Pv_flat[i]:.3g}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    flush(f"DONE {out_tag}: {int(tested.sum())}/{J} tested, n={na} ({time.time()-t0:.0f}s)")
    print(f"  [{out_tag}] DONE: n={na}, {int(tested.sum())} tested, {n_degen} degenerate, "
          f"{n_nonconv} non-conv, {n_q05} q<0.05, {len(keep)} written -> {out_c}", flush=True)
    return int(tested.sum()), n_q05, len(keep)


def main():
    global BATCH
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="comma-separated ancestries (default: all 6)")
    ap.add_argument("--batch", type=int, default=64)
    a = ap.parse_args()
    BATCH = a.batch
    only = set(a.only.split(",")) if a.only else None
    t0 = time.time()
    os.makedirs(MATDIR, exist_ok=True); os.makedirs(f"{S.OUTA}/per_ancestry", exist_ok=True)

    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    keep_cols = [j for j, c in enumerate(names) if c not in ANC_DUMMIES]
    print(f"covar cols kept per-ancestry ({len(keep_cols)}): {[names[j] for j in keep_cols]}", flush=True)
    lab_anc = ancestry_labels(X, names, N)
    uniq, cnt = np.unique(lab_anc.astype(str), return_counts=True)
    ancn = dict(zip(uniq.tolist(), cnt.tolist()))
    print("ancestry N:", ancn, flush=True)

    sv_order = []; G = S.load_G(sv_order); M = len(sv_order)
    phecodes = []; Y = S.load_disease_matrix(iidx, phecodes); J = len(phecodes)
    tiers = load_tiers(sv_order)
    label = {r["phecode"]: r["label"] for r in csv.DictReader(open(S.PSUM), delimiter="\t")}
    print(f"loaded G{G.shape} Y{Y.shape}, {M} SVs x {J} phecodes ({time.time()-t0:.0f}s)", flush=True)

    cand = [an for an in ORDER if (not only or an in only)]
    n_max = max((ancn.get(an, 0) for an in cand), default=0)
    mem_preflight(M, N, J, len(keep_cols), n_max, BATCH)
    parity_gate(G, X, Y, keep_cols)

    summary = []
    for anc in cand:
        out_tag = f"_{anc}"
        donef = f"{MATDIR}/disease_progress{out_tag}.txt"
        if os.path.exists(donef) and open(donef).read().startswith("DONE"):
            print(f"[{anc}] already DONE — skipping", flush=True); continue
        mask = (lab_anc.astype(str) == anc)
        na = int(mask.sum())
        min_case = 100 if na < 10000 else 50
        print(f"\n=== {anc}: n={na}  (min_case={min_case}"
              f"{', CONFIRMATION-ONLY' if anc not in WELL else ''}) ===", flush=True)
        if na < 400:
            print(f"  [{anc}] n<400 — skipping (nothing testable at min_case=100)", flush=True)
            summary.append((anc, na, 0, 0, 0)); continue
        Ga = np.ascontiguousarray(G[:, mask]); G2a = Ga * Ga
        Xa = np.ascontiguousarray(X[np.ix_(mask, keep_cols)])
        Ya = np.ascontiguousarray(Y[mask])
        nt, nq, nk = scan_subset(Ga, G2a, Xa, Ya, phecodes, tiers, label, sv_order, out_tag, min_case, t0)
        summary.append((anc, na, nt, nq, nk))
        del Ga, G2a, Xa, Ya; gc.collect()
        print(f"  [{anc}] peak RSS so far {peak_rss_gb():.1f} GB", flush=True)

    print("\n=== PER-ANCESTRY SCORE SUMMARY (ext2, 10,799 predictors) ===", flush=True)
    print(f"{'anc':<6}{'N':>9}{'tested':>8}{'q<0.05':>9}{'written':>9}")
    for anc, na, nt, nq, nk in summary:
        print(f"{anc:<6}{na:>9}{nt:>8}{nq:>9}{nk:>9}")
    print(f"peak RSS {peak_rss_gb():.1f} GB")
    print(f"EXT2_ALL_ANC_DONE ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
