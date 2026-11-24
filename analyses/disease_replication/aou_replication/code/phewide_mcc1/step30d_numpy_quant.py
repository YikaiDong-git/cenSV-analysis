#!/usr/bin/env python3
"""Step 30d (phewide P2 QUANT, numpy engine — FOCUSED panel): all-ancestry quant PheWAS for 6,513 SVs."""
import os, csv, glob
import numpy as np
import pgenlib
from scipy.stats import norm

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
SC = "/tmp/phewide"
INP = f"{SC}/inputs"
PGD = f"{BASE}/aou_replication/results_vds/step28_phewide/pgen"
RV = f"{BASE}/aou_replication/results_vds"
OUTA = f"{RV}/step30_phewide"; os.makedirs(OUTA, exist_ok=True)
XW = f"{RV}/step31_phewide/trait_crosswalk.tsv"
QSUM = f"{OUTA}/quant_pheno_summary.tsv"
RSV = f"{SC}/R_sv.f32"
SVCHUNK = 100; TOP_N = 50


def load_order():
    iids = [l.split()[1] for l in open(f"{SC}/phewide_pgen/sv.fam")]
    return iids, {i: k for k, i in enumerate(iids)}


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


def build_Rsv(iids, Q):
    pgens = sorted(glob.glob(f"{PGD}/svpgen300_*.pgen"),
                   key=lambda p: int(p.split("svpgen300_")[1].split(".pgen")[0]))
    N = len(iids); sv_order = []
    for pg in pgens:
        for line in open(pg[:-5] + ".pvar"):
            if not line.startswith("#"):
                sv_order.append(line.split("\t")[2])
    M = len(sv_order)
    if os.path.exists(RSV) and os.path.getsize(RSV) == M * N * 4:
        print(f"R_sv exists ({M} rows) — skip rebuild", flush=True); return sv_order
    print(f"building R_sv [{M} x {N}]", flush=True)
    R = np.memmap(RSV, dtype=np.float32, mode="w+", shape=(M, N)); buf = np.empty(N, np.float64); row = 0
    Qd = Q.astype(np.float64)
    for pg in pgens:
        r = pgenlib.PgenReader(pg.encode())
        for v in range(r.get_variant_ct()):
            r.read_dosages(v, buf)
            R[row] = (buf - Qd @ (Qd.T @ buf)).astype(np.float32); row += 1
        r.close()
    R.flush(); print(f"R_sv built ({row})", flush=True); return sv_order


def select_panel():
    mapped = {r["aou_id"] for r in csv.DictReader(open(XW), delimiter="\t") if r["aou_type"] == "loinc" and r["aou_id"]}
    top = []
    for r in csv.DictReader(open(QSUM), delimiter="\t"):
        top.append((r["loinc"], int(r["n_people"])))
    top.sort(key=lambda x: -x[1])
    panel = list(dict.fromkeys([lo for lo, _ in top if lo in mapped] + [lo for lo, _ in top[:TOP_N]]))
    return panel


def build_panel_ry(iids, iidx, Q, panel):
    hdr = open(f"{INP}/quant.pheno").readline().rstrip("\n").split("\t")
    loincs = hdr[2:]; colidx = {lo: j for j, lo in enumerate(loincs)}
    keep = [(lo, colidx[lo]) for lo in panel if lo in colidx]
    P = len(keep); N = len(iids)
    Y = np.full((N, P), np.nan, np.float32)
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
    import gc
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
    iids, iidx = load_order(); N = len(iids)
    Q = build_Q(iids, iidx); k = Q.shape[1]
    sv_order = build_Rsv(iids, Q); M = len(sv_order)
    tier = {r["sv_locus"]: r["tier"] for r in csv.DictReader(open(f"{SC}/phewide_pgen/sv_tier.tsv"), delimiter="\t")}
    panel = select_panel()
    ploincs, RyT, MkT, ss_y, n_eff = build_panel_ry(iids, iidx, Q, panel)
    P = len(ploincs)
    R = np.memmap(RSV, dtype=np.float32, mode="r", shape=(M, N))
    svs = np.array(sv_order); tiers = np.array([tier.get(s, "") for s in sv_order])

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
        Pv[flat] = pval.ravel(); Be[flat] = beta.ravel().astype(np.float32); Rc[flat] = rcorr.ravel().astype(np.float32)
        print(f"  SV {min(si+SVCHUNK,M)}/{M}", flush=True)

    Ne_row = np.broadcast_to(n_eff, (M, P)).ravel()
    valid = (Ne_row > k + 5) & np.isfinite(Pv)
    sv_of = np.repeat(np.arange(M), P)
    q = np.full(T, np.nan)
    for tn in ("signal", "null_ext"):
        idx = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if idx.size:
            q[idx] = bh(Pv[idx])
    out = np.nonzero(valid & ((Pv < 1e-3) | (np.nan_to_num(q, nan=1.0) < 0.1)))[0]
    out = out[np.argsort(Pv[out])]
    with open(f"{OUTA}/phewide_quant.tsv", "w") as f:
        f.write("sv_locus\ttier\tloinc\tn\tbeta\tpartial_r\tp\tq_bh_tier\n")
        for i in out:
            a = i // P; b = i % P
            f.write(f"{svs[a]}\t{tiers[a]}\t{ploincs[b]}\t{int(n_eff[b])}\t{Be[i]:.4g}\t{Rc[i]:.4g}\t{Pv[i]:.3g}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    with open(f"{SC}/phewide_quant_allpanel.tsv", "w") as f:
        f.write("sv_locus\ttier\tloinc\tn\tbeta\tp\tq_bh_tier\n")
        vi = np.nonzero(valid)[0]
        for i in vi:
            a = i // P; b = i % P
            f.write(f"{svs[a]}\t{tiers[a]}\t{ploincs[b]}\t{int(n_eff[b])}\t{Be[i]:.4g}\t{Pv[i]:.3g}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    sig = int(np.nansum(q < 0.05))
    print(f"QUANT_DONE: panel={P} LOINCs, {int(valid.sum())} valid tests, {sig} q<0.05 (tiered), "
          f"{len(out)} written -> {OUTA}/phewide_quant.tsv (+ full panel dump in scratch)", flush=True)


if __name__ == "__main__":
    main()
