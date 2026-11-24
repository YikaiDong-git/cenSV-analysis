#!/usr/bin/env python3
import os, sys, csv, glob, json, time, argparse
import numpy as np
import pgenlib
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phewide_ext"))
import step42_gate2_headline_rescore as GG

RV = GG.RV
SRC = f"{RV}/step28_phewide_ext2"
PGDIR = f"{SRC}/pgen"
CARRY = GG.CARRY
INP = f"{CARRY}/inputs"
OUTA = f"{RV}/step30_phewide_ext2"
MATDIR = f"{OUTA}/disease_matrices"
PSUM = f"{OUTA}/reference/phecodeX_pheno_summary.tsv"
FLUSH_EVERY = 6
VREL = 1e-6
IDENT_TOL = 1.0


def manifest_prefixes():
    rows = list(csv.DictReader(open(f"{SRC}/pgen_manifest.tsv"), delimiter="\t"))
    rows.sort(key=lambda r: int(r["order"]))
    out = [(f"{PGDIR}/{r['prefix']}", int(r["n_variants"])) for r in rows]
    for p, _ in out:
        if not os.path.exists(p + ".pgen"):
            raise SystemExit(f"ABORT: manifest lists a missing pgen: {p}.pgen")
    ng = len(glob.glob(f"{PGDIR}/*.pgen"))
    if ng != len(out):
        raise SystemExit(f"ABORT: {ng} pgens on disk but {len(out)} in the manifest — stale or partial run")
    return out


def chunk_prefixes():
    return [p for p, _ in manifest_prefixes()]


def psam_iids(pfx):
    out = []
    for line in open(f"{pfx}.psam"):
        if line.startswith("#"):
            continue
        out.append(line.split("\t")[1] if "\t" in line else line.split()[1])
    return out


def load_order():
    pfx = manifest_prefixes()[0][0]
    iids = psam_iids(pfx)
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
    assert seen == N, f"covar.txt covered {seen}/{N} cohort samples"
    print(f"X {X.shape} filled ({seen} covar rows)", flush=True)
    return X, names


def load_G(sv_order_out):
    import hashlib
    pre = manifest_prefixes()
    sig0 = None
    for pfx, _ in pre:
        h = hashlib.md5()
        with open(f"{pfx}.psam", "rb") as fh:
            for line in fh:
                if not line.startswith(b"#"):
                    h.update(line.split(b"\t")[1] if b"\t" in line else line.split()[1]); h.update(b"\n")
        s = h.hexdigest()
        if sig0 is None:
            sig0 = s
        elif s != sig0:
            raise SystemExit(f"ABORT: psam sample order mismatch at {pfx}.psam — cannot concatenate")
    M = sum(n for _, n in pre)
    N = len(psam_iids(pre[0][0]))
    print(f"allocating G [{M} x {N}] float32 = {M*N*4/1e9:.1f} GB", flush=True)
    G = np.empty((M, N), np.float32)
    buf = np.empty(N, np.float64)
    off = 0; t0 = time.time()
    for k, (pfx, nv) in enumerate(pre):
        ids = [l.split("\t")[2] for l in open(f"{pfx}.pvar") if not l.startswith("#")]
        assert len(ids) == nv, f"{pfx}: pvar has {len(ids)} vars, manifest says {nv}"
        r = pgenlib.PgenReader(f"{pfx}.pgen".encode())
        for v in range(nv):
            r.read_dosages(v, buf); G[off + v, :] = buf
        r.close()
        sv_order_out.extend(ids); off += nv
        if (k + 1) % 10 == 0:
            print(f"  loaded {off}/{M} SVs ({time.time()-t0:.0f}s)", flush=True)
    assert off == M
    return G


def load_disease_matrix(iidx, phecodes_out):
    t0 = time.time()
    hdr = open(f"{INP}/disease.pheno").readline().rstrip("\n").split("\t")
    phecodes = hdr[2:]
    J = len(phecodes); N = len(iidx)
    Y = np.full((N, J), -1, np.int8)
    n = 0
    with open(f"{INP}/disease.pheno") as f:
        next(f)
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


def score_test_batch(G, G2, X, Yb, min_case=50):
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
    return Be, Pv, ncase, npres, tested, nonconv, degen


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


def save_matrices(Pv, Be, phecodes, label, ncase, npres, tested, sv_order, tiers, progress):
    def _save(path, arr):
        np.save(path + ".tmp.npy", arr); os.replace(path + ".tmp.npy", path)
    _save(f"{MATDIR}/disease_Pv.npy", Pv)
    _save(f"{MATDIR}/disease_Be.npy", Be)
    with open(f"{MATDIR}/disease_meta.tsv.tmp", "w") as f:
        f.write("phecode\tlabel\tcases\tn\ttested\n")
        for b in range(len(phecodes)):
            f.write(f"{phecodes[b]}\t{label.get(phecodes[b],'')}\t{ncase[b]}\t{npres[b]}\t{int(tested[b])}\n")
    os.replace(f"{MATDIR}/disease_meta.tsv.tmp", f"{MATDIR}/disease_meta.tsv")
    with open(f"{MATDIR}/disease_svorder.txt.tmp", "w") as f:
        f.write("\n".join(f"{sv_order[a]}\t{tiers[a]}" for a in range(len(sv_order))) + "\n")
    os.replace(f"{MATDIR}/disease_svorder.txt.tmp", f"{MATDIR}/disease_svorder.txt")
    with open(f"{MATDIR}/disease_progress.txt", "w") as f:
        f.write(progress + "\n")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--batch", type=int, default=64)
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(OUTA, exist_ok=True); os.makedirs(MATDIR, exist_ok=True)
    GG.guard(OUTA); GG.guard(MATDIR)

    iids, iidx = load_order(); N = len(iids)
    print(f"cohort N={N:,}", flush=True)
    X, _ = load_X(iidx)
    sv_order = []; G = load_G(sv_order); M = len(sv_order)
    print(f"G loaded: {M:,} SVs. computing G2 ...", flush=True)
    G2 = G * G
    phecodes = []; Y = load_disease_matrix(iidx, phecodes); J = len(phecodes)

    tier = {r["sv_locus"]: r["tier"] for r in csv.DictReader(open(f"{SRC}/sv_tier.tsv"), delimiter="\t")}
    missing = [s for s in sv_order if not tier.get(s, "")]
    if missing:
        raise SystemExit(f"ABORT: {len(missing)} scored SVs have no tier (e.g. {missing[:3]}). "
                         f"step30g would have silently dropped them from FDR — refusing to repeat that.")
    tiers = np.array([tier[s] for s in sv_order])
    from collections import Counter
    tier_counts = Counter(tiers.tolist())
    print(f"tiers present: {dict(tier_counts)}", flush=True)
    label = {r["phecode"]: r["label"] for r in csv.DictReader(open(PSUM), delimiter="\t")}

    Pv = np.full((M, J), np.nan, np.float32); Be = np.full((M, J), np.nan, np.float32)
    ncase = np.zeros(J, np.int32); npres = np.zeros(J, np.int32); tested = np.zeros(J, bool)
    n_nonconv = 0; n_degen = 0

    pv_path = f"{MATDIR}/disease_Pv.npy"; ord_path = f"{MATDIR}/disease_svorder.txt"
    meta_path = f"{MATDIR}/disease_meta.tsv"
    if all(os.path.exists(p) for p in (pv_path, ord_path, meta_path)):
        prev_order = [ln.split("\t")[0] for ln in open(ord_path) if ln.strip()]
        prevPv = np.load(pv_path)
        if prev_order == sv_order and prevPv.shape == (M, J):
            Pv[:] = prevPv; Be[:] = np.load(f"{MATDIR}/disease_Be.npy")
            for b, r in enumerate(csv.DictReader(open(meta_path), delimiter="\t")):
                tested[b] = r["tested"] == "1"; ncase[b] = int(r["cases"]); npres[b] = int(r["n"])
            print(f"RESUME: {int(tested.sum())}/{J} phecodes already done", flush=True)
        else:
            print("RESUME: partial incompatible — fresh run", flush=True)

    BATCH = a.batch
    print(f"scoring {J} phecodes x {M} SVs (batch={BATCH}) ...", flush=True)
    for bi, j0 in enumerate(range(0, J, BATCH)):
        cols = [j for j in range(j0, min(j0 + BATCH, J)) if not tested[j]]
        if not cols:
            continue
        Be_b, Pv_b, nca_b, npr_b, tt_b, nc_b, dg_b = score_test_batch(G, G2, X, Y[:, cols])
        n_nonconv += int(nc_b.sum()); n_degen += int(dg_b.sum())
        for k, j in enumerate(cols):
            ncase[j] = nca_b[k]; npres[j] = npr_b[k]
            if not tt_b[k]:
                continue
            Be[:, j] = Be_b[:, k]; Pv[:, j] = Pv_b[:, k]; tested[j] = True
        print(f"  {min(j0+BATCH,J)}/{J} scanned, {int(tested.sum())} tested ({time.time()-t0:.0f}s)",
              flush=True)
        if (bi + 1) % FLUSH_EVERY == 0:
            save_matrices(Pv, Be, phecodes, label, ncase, npres, tested, sv_order, tiers,
                          f"in_progress: {min(j0+BATCH,J)}/{J} ({time.time()-t0:.0f}s)")
            print(f"    [flush] -> {MATDIR}", flush=True)

    save_matrices(Pv, Be, phecodes, label, ncase, npres, tested, sv_order, tiers,
                  f"DONE {J}/{J} ({time.time()-t0:.0f}s)")

    sv_of = np.repeat(np.arange(M), J); ph_of = np.tile(np.arange(J), M)
    Pv_flat = Pv.ravel(); Be_flat = Be.ravel()
    valid = np.isfinite(Pv_flat) & tested[ph_of]
    q = np.full(M * J, np.nan)
    per_tier = {}
    for tn in sorted(tier_counts):
        idx = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if idx.size:
            q[idx] = bh(Pv_flat[idx])
            per_tier[tn] = dict(n_tests=int(idx.size), n_q05=int(np.nansum(q[idx] < 0.05)))
    n_q05 = int(np.nansum(q < 0.05))
    print(f"\ntiered BH: {n_q05} cells q<0.05", flush=True)
    for tn, d in per_tier.items():
        print(f"  {tn:12s} {d['n_tests']:>10,} tests  {d['n_q05']:>5} q<0.05", flush=True)

    keep = np.nonzero(valid & ((np.nan_to_num(q, nan=1.0) < 0.05) | (Pv_flat < 1e-4)))[0]
    keep = keep[np.argsort(Pv_flat[keep])]
    out_c = f"{OUTA}/phewide_disease_full.tsv"
    with open(GG.guard(out_c), "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tbeta_logOR\tchi2_p\tq_bh_tier\n")
        for i in keep:
            s = sv_of[i]; j = ph_of[i]
            f.write(f"{sv_order[s]}\t{tiers[s]}\t{phecodes[j]}\t{label.get(phecodes[j],'')}\t"
                    f"{ncase[j]}\t{npres[j]}\t{Be_flat[i]:.4g}\t{Pv_flat[i]:.4g}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")

    lam = float(np.nanmedian(chi2dist.isf(Pv_flat[valid], 1)) / chi2dist.ppf(0.5, 1))
    summ = dict(M=M, J=J, n_tested_phecodes=int(tested.sum()), n_valid_tests=int(valid.sum()),
                n_q05=n_q05, per_tier=per_tier, lambda_GC=lam,
                n_nonconv=n_nonconv, n_degen=n_degen, n_hits_written=int(keep.size),
                secs=time.time() - t0)
    with open(GG.guard(f"{OUTA}/step51_summary.json"), "w") as f:
        json.dump(summ, f, indent=1)
    print(json.dumps(summ, indent=1), flush=True)
    print(f"\nSTEP51_DONE in {(time.time()-t0)/60:.1f} min -> {out_c}", flush=True)


if __name__ == "__main__":
    main()
