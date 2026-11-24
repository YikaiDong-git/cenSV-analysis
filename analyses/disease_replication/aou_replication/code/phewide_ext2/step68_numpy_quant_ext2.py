#!/usr/bin/env python3
"""Step 68 (ext2 arm) — all-ancestry QUANTITATIVE-trait PheWAS over the 10,799 ext2 predictors (← step30d)."""
import os, sys, csv, time, gc
import numpy as np
import pgenlib
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S

RV = S.RV
OUTA = S.OUTA
INP = S.INP
SVTIER = f"{S.SRC}/sv_tier.tsv"
XW = f"{OUTA}/reference/trait_crosswalk.tsv"
QSUM = f"{OUTA}/reference/quant_pheno_summary.tsv"
SCRATCH = "/tmp/phewide_ext2"
RSV = f"{SCRATCH}/R_sv.f32"
QDIR = f"{OUTA}/quant"
SVCHUNK = 100
TOP_N = 50


def build_Q(iids, iidx):
    N = len(iids)
    with open(f"{INP}/covar.txt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = [j for j, c in enumerate(hdr) if c not in ("FID", "IID")]; iidcol = hdr.index("IID")
        C = np.zeros((N, len(ci) + 1), np.float64); C[:, 0] = 1.0
        seen = 0
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[iidcol])
            if kk is None:
                continue
            for jj, j in enumerate(ci):
                C[kk, jj + 1] = float(p[j])
            seen += 1
    assert seen == N, (f"covar.txt covered {seen}/{N} cohort samples — a missing sample would leave a "
                       f"zero covariate row and bias the QR basis")
    Q, _ = np.linalg.qr(C)
    print(f"covariate basis Q {Q.shape}", flush=True)
    return Q.astype(np.float32)


def build_Rsv(N, Q):
    pre = S.manifest_prefixes()
    sv_order = []
    for pfx, _ in pre:
        for line in open(pfx + ".pvar"):
            if not line.startswith("#"):
                sv_order.append(line.split("\t")[2])
    M = len(sv_order)
    os.makedirs(SCRATCH, exist_ok=True)
    if os.path.exists(RSV) and os.path.getsize(RSV) == M * N * 4:
        print(f"R_sv exists ({M} rows) — skip rebuild", flush=True)
        return sv_order
    need = M * N * 4
    free = os.statvfs(SCRATCH).f_bavail * os.statvfs(SCRATCH).f_frsize
    print(f"building R_sv [{M} x {N}] = {need/2**30:.1f} GB (free {free/2**30:.1f} GB)", flush=True)
    if free < need * 1.05:
        raise SystemExit(f"ABORT: need {need/2**30:.1f} GB for R_sv, only {free/2**30:.1f} GB free")
    t0 = time.time()
    R = np.memmap(RSV, dtype=np.float32, mode="w+", shape=(M, N))
    buf = np.empty(N, np.float64); row = 0
    Qd = Q.astype(np.float64)
    for pfx, _ in pre:
        r = pgenlib.PgenReader((pfx + ".pgen").encode())
        for v in range(r.get_variant_ct()):
            r.read_dosages(v, buf)
            R[row] = (buf - Qd @ (Qd.T @ buf)).astype(np.float32); row += 1
            if row % 2000 == 0:
                print(f"  R_sv {row}/{M} ({time.time()-t0:.0f}s)", flush=True)
        r.close()
    R.flush(); del R; gc.collect()
    print(f"R_sv built ({row} rows, {time.time()-t0:.0f}s)", flush=True)
    return sv_order


def select_panel():
    mapped = {r["aou_id"] for r in csv.DictReader(open(XW), delimiter="\t")
              if r["aou_type"] == "loinc" and r["aou_id"]}
    top = [(r["loinc"], int(r["n_people"])) for r in csv.DictReader(open(QSUM), delimiter="\t")]
    top.sort(key=lambda x: -x[1])
    panel = list(dict.fromkeys([lo for lo, _ in top if lo in mapped] + [lo for lo, _ in top[:TOP_N]]))
    print(f"panel: {len(panel)} LOINCs ({len(mapped)} discovery-mapped, top-{TOP_N} by coverage)", flush=True)
    return panel


def build_panel_ry(iids, iidx, Q, panel):
    hdr = open(f"{INP}/quant.pheno").readline().rstrip("\n").split("\t")
    loincs = hdr[2:]; colidx = {lo: j for j, lo in enumerate(loincs)}
    keep = [(lo, colidx[lo]) for lo in panel if lo in colidx]
    P = len(keep); N = len(iids)
    Y = np.full((N, P), np.nan, np.float32)
    t0 = time.time(); n = 0
    with open(f"{INP}/quant.pheno") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[0])
            if kk is None:
                continue
            for a, (lo, j) in enumerate(keep):
                v = p[2 + j]
                if v != "NA":
                    Y[kk, a] = float(v)
            n += 1
            if n % 100000 == 0:
                print(f"  quant parse {n} rows ({time.time()-t0:.0f}s)", flush=True)
    Msk = np.isfinite(Y).astype(np.float32)
    np.nan_to_num(Y, copy=False)
    Y -= Q @ (Q.T @ Y)
    Y *= Msk
    ss_y = np.einsum("ij,ij->j", Y, Y).astype(np.float64); n_eff = Msk.sum(0).astype(np.float64)
    RyT = np.ascontiguousarray(Y.T); MkT = np.ascontiguousarray(Msk.T)
    del Y, Msk; gc.collect()
    print(f"panel residuals [{P} x {N}] in RAM ({RyT.nbytes/1e6:.0f} MB x2)", flush=True)
    return [lo for lo, _ in keep], RyT, MkT, ss_y, n_eff


def bh(pvals):
    ps = np.asarray(pvals); ok = np.isfinite(ps); n = int(ok.sum())
    q = np.full(len(ps), np.nan); order = np.argsort(np.where(ok, ps, 2)); prev = 1.0
    for rank, i in enumerate(order[::-1]):
        if not ok[i]:
            continue
        prev = min(prev, ps[i] * n / (n - rank)); q[i] = prev
    return q


def main():
    t0 = time.time()
    os.makedirs(QDIR, exist_ok=True)
    iids, iidx = S.load_order(); N = len(iids)
    Q = build_Q(iids, iidx); k = Q.shape[1]
    sv_order = build_Rsv(N, Q); M = len(sv_order)

    tmap = {r["sv_locus"]: r["tier"] for r in csv.DictReader(open(SVTIER), delimiter="\t")}
    missing = [s for s in sv_order if not tmap.get(s)]
    if missing:
        raise SystemExit(f"ABORT: {len(missing)} scored SVs have no tier (e.g. {missing[:3]}) — "
                         f"fail-CLOSED (step30d's tier.get(s,'') would have dropped them from BH)")
    tiers = np.array([tmap[s] for s in sv_order])
    from collections import Counter
    print(f"tiers (data-derived, fail-closed): {dict(Counter(tiers.tolist()))}", flush=True)

    panel = select_panel()
    ploincs, RyT, MkT, ss_y, n_eff = build_panel_ry(iids, iidx, Q, panel)
    P = len(ploincs)
    R = np.memmap(RSV, dtype=np.float32, mode="r", shape=(M, N))
    svs = np.array(sv_order)

    T = M * P
    Pv = np.full(T, np.nan, np.float64); Be = np.zeros(T, np.float32); Rc = np.zeros(T, np.float32)
    for si in range(0, M, SVCHUNK):
        Rs = np.array(R[si:si + SVCHUNK], np.float32)
        num = (Rs @ RyT.T).astype(np.float64)
        denx = ((Rs * Rs) @ MkT.T).astype(np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            rcorr = num / np.sqrt(denx * ss_y[None, :]); beta = num / denx
        df = n_eff[None, :] - k - 1
        tstat = rcorr * np.sqrt(np.maximum(df, 1) / np.maximum(1 - rcorr ** 2, 1e-12))
        pval = 2 * norm.sf(np.abs(tstat))
        srows = np.arange(si, si + Rs.shape[0])
        flat = (srows[:, None] * P + np.arange(P)[None, :]).ravel()
        Pv[flat] = pval.ravel(); Be[flat] = beta.ravel().astype(np.float32)
        Rc[flat] = rcorr.ravel().astype(np.float32)
        if (si // SVCHUNK) % 20 == 0:
            print(f"  SV {min(si+SVCHUNK,M)}/{M} ({time.time()-t0:.0f}s)", flush=True)

    Ne_row = np.broadcast_to(n_eff, (M, P)).ravel()
    valid = (Ne_row > k + 5) & np.isfinite(Pv)
    sv_of = np.repeat(np.arange(M), P)
    q = np.full(T, np.nan)
    for tn in sorted({t for t in tiers.tolist() if t}):
        idx = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if idx.size:
            q[idx] = bh(Pv[idx])
    out = np.nonzero(valid & ((Pv < 1e-3) | (np.nan_to_num(q, nan=1.0) < 0.1)))[0]
    out = out[np.argsort(Pv[out])]
    hits = f"{QDIR}/phewide_quant_ext2.tsv"
    with open(hits, "w") as f:
        f.write("sv_locus\ttier\tloinc\tn\tbeta\tpartial_r\tp\tq_bh_tier\n")
        for i in out:
            a = i // P; b = i % P
            f.write(f"{svs[a]}\t{tiers[a]}\t{ploincs[b]}\t{int(n_eff[b])}\t{Be[i]:.4g}\t{Rc[i]:.4g}\t"
                    f"{Pv[i]:.3g}\t{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    dump = f"{QDIR}/phewide_quant_allpanel_ext2.tsv"
    with open(dump, "w") as f:
        f.write("sv_locus\ttier\tloinc\tn\tbeta\tp\tq_bh_tier\n")
        for i in np.nonzero(valid)[0]:
            a = i // P; b = i % P
            f.write(f"{svs[a]}\t{tiers[a]}\t{ploincs[b]}\t{int(n_eff[b])}\t{Be[i]:.4g}\t{Pv[i]:.3g}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    sig = int(np.nansum(q < 0.05))
    print(f"\nQUANT_EXT2_DONE: panel={P} LOINCs x {M} SVs, {int(valid.sum())} valid tests, "
          f"{sig} q<0.05 (tiered), {len(out)} written", flush=True)
    print(f"  -> {hits}")
    print(f"  -> {dump}  ({os.path.getsize(dump)/1e6:.1f} MB, bucket only)")
    print(f"  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
