#!/usr/bin/env python3
import os, csv, sys, time, glob
import numpy as np
import pgenlib
from scipy.stats import chi2 as chi2dist

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
SC = "/tmp/phewide_mcc1"
INP = f"{SC}/inputs"
PGDIR = f"{SC}/pgen"
RV = f"{BASE}/aou_replication/results_vds"
OUTA = f"{RV}/step30_phewide_mcc1"; os.makedirs(OUTA, exist_ok=True)
assert "mcc1" in SC and "mcc1" in OUTA, (
    f"MCC1 GUARD TRIPPED: SC={SC!r} OUTA={OUTA!r} — refusing to run against the non-mcc1 tree")
MATDIR = f"{OUTA}/disease_matrices"; os.makedirs(MATDIR, exist_ok=True)
PSUM = f"{OUTA}/phecodeX_pheno_summary.tsv"
FLUSH_EVERY = 6
VREL = 1e-6
IDENT_TOL = 1.0


def chunk_prefixes():
    fs = sorted(glob.glob(f"{PGDIR}/svpgen300_*.pgen"),
                key=lambda p: int(p.split("svpgen300_")[1].split(".pgen")[0]))
    return [f[:-len(".pgen")] for f in fs]


def load_order():
    p0 = chunk_prefixes()[0]
    iids = []
    for line in open(f"{p0}.psam"):
        if line.startswith("#"):
            continue
        iids.append(line.split("\t")[1] if "\t" in line else line.split()[1])
    return iids, {i: k for k, i in enumerate(iids)}


def load_X(iidx):
    N = len(iidx)
    with open(f"{INP}/covar.txt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = [j for j, c in enumerate(hdr) if c not in ("FID", "IID")]
        names = ["intercept"] + [hdr[j] for j in ci]
        iidcol = hdr.index("IID")
        X = np.zeros((N, len(ci) + 1), np.float64); X[:, 0] = 1.0
        seen = 0
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[iidcol])
            if kk is None:
                continue
            for jj, j in enumerate(ci):
                X[kk, jj + 1] = float(p[j])
            seen += 1
    assert seen == N, (f"covar.txt covered {seen}/{N} cohort samples — any missing sample would keep an "
                       f"all-zero covariate row (intercept only) and silently bias the null fit")
    print(f"X {X.shape} filled ({seen} covar rows)", flush=True)
    return X, names


def load_G(sv_order_out):
    prefixes = chunk_prefixes()
    import hashlib
    def _psam_sig(pfx):
        h = hashlib.md5()
        with open(f"{pfx}.psam", "rb") as fh:
            for line in fh:
                if not line.startswith(b"#"):
                    h.update(line.split(b"\t")[1] if b"\t" in line else line.split()[1]); h.update(b"\n")
        return h.hexdigest()
    sig0 = _psam_sig(prefixes[0])
    for pfx in prefixes[1:]:
        if _psam_sig(pfx) != sig0:
            raise SystemExit(f"ABORT: psam sample order mismatch at {pfx}.psam vs {prefixes[0]}.psam "
                             f"— chunk pgens are not sample-aligned; cannot concatenate dosages.")
    idlists = []
    for pfx in prefixes:
        ids = [l.split("\t")[2] for l in open(f"{pfx}.pvar") if not l.startswith("#")]
        idlists.append(ids)
    M = sum(len(x) for x in idlists)
    r0 = pgenlib.PgenReader(f"{prefixes[0]}.pgen".encode()); N = r0.get_raw_sample_ct(); r0.close()
    G = np.empty((M, N), np.float32); buf = np.empty(N, np.float64)
    row = 0
    for pfx, ids in zip(prefixes, idlists):
        r = pgenlib.PgenReader(f"{pfx}.pgen".encode())
        for v in range(len(ids)):
            r.read_dosages(v, buf)
            G[row] = buf; row += 1
        r.close()
        sv_order_out.extend(ids)
    assert row == M
    print(f"G {G.shape} loaded from {len(prefixes)} chunk pgens ({G.nbytes/1e9:.1f} GB)", flush=True)
    return G


def load_disease_matrix(iidx, phecodes_out):
    N = len(iidx)
    with open(f"{INP}/disease.pheno") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        phecodes = hdr[2:]; J = len(phecodes)
        Y = np.full((N, J), -1, np.int8)
        t0 = time.time(); n = 0
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[1])
            if kk is None:
                continue
            row = p[2:]
            for a in range(J):
                v = row[a]
                if v == "0":
                    Y[kk, a] = 0
                elif v == "1":
                    Y[kk, a] = 1
            n += 1
            if n % 50000 == 0:
                print(f"  disease parse {n} rows ({time.time()-t0:.0f}s)", flush=True)
    phecodes_out.extend(phecodes)
    print(f"Y {Y.shape} loaded ({Y.nbytes/1e6:.0f} MB)", flush=True)
    return Y


def fit_null(Xk, y, present):
    kk = Xk.shape[1]
    beta = np.zeros(kk)
    ybar = float(y[present].mean())
    ybar = min(max(ybar, 1e-6), 1 - 1e-6)
    beta[0] = np.log(ybar / (1 - ybar))
    converged = False
    for _ in range(25):
        eta = np.clip(Xk @ beta, -30, 30)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = p * (1 - p); w[~present] = 0.0
        r = y - p; r[~present] = 0.0
        WX = Xk * w[:, None]
        XtWX = Xk.T @ WX
        grad = Xk.T @ r
        try:
            delta = np.linalg.solve(XtWX, grad)
        except np.linalg.LinAlgError:
            delta = np.linalg.lstsq(XtWX, grad, rcond=None)[0]
        beta += delta
        if np.max(np.abs(delta)) < 1e-6:
            converged = True
            break
    eta = np.clip(Xk @ beta, -30, 30); p = 1.0 / (1.0 + np.exp(-eta))
    w = p * (1 - p); w[~present] = 0.0
    r = y - p; r[~present] = 0.0
    return p, w, r, converged


def score_test_batch(G, G2, X, Yb, min_case=50, return_uv=False):
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
        _, w, r, conv = fit_null(X, y, present)
        nonconv[b] = not conv
        if abs(float(r.sum())) > IDENT_TOL:
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
    Umat = np.full((M, B), np.nan, np.float32) if return_uv else None
    Vmat = np.full((M, B), np.nan, np.float32) if return_uv else None
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
        good = np.isfinite(stat) & (V > VREL * Vraw) & (Vraw > 0)
        stat = np.where(good, stat, np.nan)
        Pv[:, b] = chi2dist.sf(stat, 1); Be[:, b] = np.where(good, beta, np.nan).astype(np.float32)
        if return_uv:
            Umat[:, b] = np.where(good, U, np.nan).astype(np.float32)
            Vmat[:, b] = np.where(good, V, np.nan).astype(np.float32)
    if return_uv:
        return Be, Pv, ncase, npres, tested, nonconv, degen, Umat, Vmat
    return Be, Pv, ncase, npres, tested, nonconv, degen


def score_test_pheno(G, G2, X, y_int, min_case=50):
    Be, Pv, ncase, npres, tested, _nonconv, _degen = score_test_batch(G, G2, X, y_int[:, None].astype(np.int8), min_case)
    if not tested[0]:
        return None
    return Be[:, 0], np.full(G.shape[0], np.nan, np.float32), Pv[:, 0], int(npres[0]), int(ncase[0])


def bh(pvals):
    ps = np.asarray(pvals, float); ok = np.isfinite(ps); n = int(ok.sum())
    q = np.full(len(ps), np.nan)
    if n == 0:
        return q
    ok_idx = np.nonzero(ok)[0]
    order = ok_idx[np.argsort(ps[ok_idx])]
    prev = 1.0
    for k, i in enumerate(order[::-1]):
        prev = min(prev, ps[i] * n / (n - k)); q[i] = prev
    return q


def save_matrices(Pv, Be, phecodes, label, ncase, npres, tested, sv_order, tiers, out_tag, progress):
    def _save_npy(path, arr):
        np.save(path + ".tmp.npy", arr); os.replace(path + ".tmp.npy", path)
    _save_npy(f"{MATDIR}/disease_Pv{out_tag}.npy", Pv)
    _save_npy(f"{MATDIR}/disease_Be{out_tag}.npy", Be)
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


def main(pheno_subset=None, out_tag=""):
    t0 = time.time()
    iids, iidx = load_order(); N = len(iids)
    X, _ = load_X(iidx)
    sv_order = []; G = load_G(sv_order); M = len(sv_order)
    print("computing G2 = G*G ...", flush=True); G2 = G * G
    phecodes = []; Y = load_disease_matrix(iidx, phecodes); J = len(phecodes)
    svtier_path = next((p for p in (f"{SC}/sv_tier.tsv", f"{INP}/sv_tier.tsv") if os.path.exists(p)), None)
    if svtier_path is None:
        raise SystemExit(f"ABORT: sv_tier.tsv not found at {SC}/ or {INP}/ — cannot assign SV tiers.")
    tier = {r["sv_locus"]: r["tier"] for r in csv.DictReader(open(svtier_path), delimiter="\t")}
    label = {r["phecode"]: r["label"] for r in csv.DictReader(open(PSUM), delimiter="\t")}
    tiers = np.array([tier.get(s, "") for s in sv_order])

    sub = set(pheno_subset) if pheno_subset else None
    Pv = np.full((M, J), np.nan, np.float32); Be = np.full((M, J), np.nan, np.float32)
    ncase = np.zeros(J, np.int32); npres = np.zeros(J, np.int32); tested = np.zeros(J, bool)
    n_nonconv = 0; n_degen = 0

    if sub is None:
        pv_path = f"{MATDIR}/disease_Pv{out_tag}.npy"; ord_path = f"{MATDIR}/disease_svorder{out_tag}.txt"
        meta_path = f"{MATDIR}/disease_meta{out_tag}.tsv"
        if os.path.exists(pv_path) and os.path.exists(ord_path) and os.path.exists(meta_path):
            prev_order = [ln.split("\t")[0] for ln in open(ord_path) if ln.strip()]
            prevPv = np.load(pv_path)
            if prev_order == sv_order and prevPv.shape == (M, J):
                Pv[:] = prevPv; Be[:] = np.load(f"{MATDIR}/disease_Be{out_tag}.npy")
                mrows = list(csv.DictReader(open(meta_path), delimiter="\t"))
                for b, r in enumerate(mrows):
                    tested[b] = r["tested"] == "1"; ncase[b] = int(r["cases"]); npres[b] = int(r["n"])
                print(f"RESUME: loaded partial, {int(tested.sum())}/{J} phecodes already done — skipping them",
                      flush=True)
            else:
                print("RESUME: existing partial is incompatible (sv_order/shape) — ignoring, fresh run", flush=True)

    BATCH = 64
    print(f"scoring {J if sub is None else len(sub)} phecodes x {M} SVs (batch={BATCH}) ...", flush=True)
    for bi, j0 in enumerate(range(0, J, BATCH)):
        cols = list(range(j0, min(j0 + BATCH, J)))
        if sub is not None:
            cols = [j for j in cols if phecodes[j] in sub]
        else:
            cols = [j for j in cols if not tested[j]]
        if not cols:
            continue
        Be_b, Pv_b, nca_b, npr_b, tt_b, nc_b, dg_b = score_test_batch(G, G2, X, Y[:, cols])
        n_nonconv += int(nc_b.sum()); n_degen += int(dg_b.sum())
        for a, j in enumerate(cols):
            if not tt_b[a]:
                ncase[j] = nca_b[a]; npres[j] = npr_b[a]
                continue
            Be[:, j] = Be_b[:, a]; Pv[:, j] = Pv_b[:, a]; ncase[j] = nca_b[a]; npres[j] = npr_b[a]
            tested[j] = True
        print(f"  {min(j0+BATCH,J)}/{J} scanned, {int(tested.sum())} tested ({time.time()-t0:.0f}s)", flush=True)
        if (bi + 1) % FLUSH_EVERY == 0:
            save_matrices(Pv, Be, phecodes, label, ncase, npres, tested, sv_order, tiers, out_tag,
                          f"in_progress: {min(j0+BATCH,J)}/{J} phecode-cols done ({time.time()-t0:.0f}s)")
            print(f"    [flush] partial matrices -> {MATDIR} ({time.time()-t0:.0f}s)", flush=True)

    sv_of = np.repeat(np.arange(M), J); ph_of = np.tile(np.arange(J), M)
    Pv_flat = Pv.ravel(); Be_flat = Be.ravel()
    valid = np.isfinite(Pv_flat) & tested[ph_of]
    q = np.full(M * J, np.nan)
    for tn in ("signal", "null_ext"):
        idx = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if idx.size:
            q[idx] = bh(Pv_flat[idx])
    n_q05 = int(np.nansum(q < 0.05))
    keep = np.nonzero(valid & ((np.nan_to_num(q, nan=1.0) < 0.05) | (Pv_flat < 1e-4)))[0]
    keep = keep[np.argsort(Pv_flat[keep])]
    out_c = f"{OUTA}/phewide_disease_full{out_tag}.tsv"
    with open(out_c, "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tbeta_logOR\tchi2_p\tq_bh_tier\n")
        for i in keep:
            a = i // J; b = i % J
            f.write(f"{sv_order[a]}\t{tiers[a]}\t{phecodes[b]}\t{label.get(phecodes[b],'')}\t{ncase[b]}\t"
                    f"{npres[b]}\t{Be_flat[i]:.4g}\t{Pv_flat[i]:.3g}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    save_matrices(Pv, Be, phecodes, label, ncase, npres, tested, sv_order, tiers, out_tag,
                  f"DONE: {int(tested.sum())}/{J} phecodes tested ({time.time()-t0:.0f}s)")
    print(f"DISEASE_DONE: {int(tested.sum())} phecodes tested, {int(valid.sum())} valid tests, "
          f"{n_nonconv} null-fits non-converged, {n_degen} degenerate phenos MASKED, "
          f"{n_q05} q<0.05 (tiered), {len(keep)} written -> {out_c}; "
          f"matrices -> {MATDIR}/disease_*.npy ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
