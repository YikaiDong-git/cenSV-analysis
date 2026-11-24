#!/usr/bin/env python3
"""finemap_core — numpy engine for SuSiE fine-mapping of AoU centromeric-SV disease hits."""
from __future__ import annotations

import csv
import fcntl
import glob
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
import pgenlib

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
FM = f"{BASE}/aou_replication/finemap"
REPVDS = f"{BASE}/aou_replication/results_vds"
SVPGEN_DIR = f"{REPVDS}/step28_phewide_ext2/pgen"
CENTRES = f"{HOME}/workspace/working/jwein22/centromere_phewas/results"
PERI_BFILE = CENTRES + "/vds_pericentromeric_snps_min_ac1000/{chrom}/{chrom}.pericentromeric_snps.min_ac1000"
PRED_IDS_DIR = f"{FM}/inputs_ext2/pred_ids"


def load_pred_ids(chrom):
    p = f"{PRED_IDS_DIR}/{chrom}.pred_ids.txt"
    assert os.path.exists(p), f"missing ext2 EN-predictor ids: {p} (run m6_build_pred_ids_ext2.py)"
    with open(p) as fh:
        ids = {ln.strip() for ln in fh if ln.strip()}
    assert ids, f"empty ext2 EN-predictor id set: {p}"
    return ids
PLINK2 = os.environ.get("PLINK2", f"{HOME}/bin/plink2")

KEEP_DIR = f"{FM}/inputs_ext2/keep"
COVAR_DIR = f"{FM}/inputs_ext2/covar"
PHENO_DIR = f"{FM}/inputs_ext2/pheno"

VREL = 1e-6
IDENT_TOL = 1.0


def log(m):
    print(m, flush=True)


def _psam_iids(psam):
    out = []
    with open(psam) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            out.append(line.split("\t")[1] if "\t" in line else line.split()[1])
    return out


def load_keep(stratum):
    path = f"{KEEP_DIR}/{stratum}.keep.fid_iid.tsv"
    iids = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1] != "IID":
                iids.append(p[1])
    return iids, {i: k for k, i in enumerate(iids)}


def load_covar(stratum, iidx):
    N = len(iidx)
    path = f"{COVAR_DIR}/{stratum}.covar.txt"
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = [j for j, c in enumerate(hdr) if c not in ("FID", "IID")]
        names = ["intercept"] + [hdr[j] for j in ci]
        iidcol = hdr.index("IID")
        X = np.zeros((N, len(ci) + 1), np.float64)
        X[:, 0] = 1.0
        seen = 0
        for line in f:
            p = line.rstrip("\n").split("\t")
            kk = iidx.get(p[iidcol])
            if kk is None:
                continue
            for jj, j in enumerate(ci):
                X[kk, jj + 1] = float(p[j])
            seen += 1
    assert seen == N, f"covar covered {seen}/{N} keep samples — missing rows would bias the null fit"
    return X, names


def load_pheno(mcc, phecode, iidx):
    N = len(iidx)
    path = f"{PHENO_DIR}/disease.{mcc}.pheno"
    y = np.full(N, -1, np.int8)
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        if phecode not in hdr:
            raise KeyError(f"{phecode} not in {path}")
        col = hdr.index(phecode)
        iidcol = hdr.index("IID")
        for line in f:
            p = line.rstrip("\n").split("\t")
            kk = iidx.get(p[iidcol])
            if kk is None:
                continue
            v = p[col]
            if v == "1":
                y[kk] = 1
            elif v == "0":
                y[kk] = 0
    return y


def _svpgen_prefixes():
    fs = sorted(glob.glob(f"{SVPGEN_DIR}/svpgen_*.pgen"), key=os.path.basename)
    assert fs, f"no ext2 svpgen chunks under {SVPGEN_DIR}"
    return [f[:-len(".pgen")] for f in fs]


def _svpgen_index():
    prefixes = _svpgen_prefixes()

    def sig(pfx):
        h = hashlib.md5()
        with open(f"{pfx}.psam", "rb") as fh:
            for line in fh:
                if not line.startswith(b"#"):
                    h.update(line.split(b"\t")[1] if b"\t" in line else line.split()[1])
                    h.update(b"\n")
        return h.hexdigest()

    s0 = sig(prefixes[0])
    for pfx in prefixes[1:]:
        if sig(pfx) != s0:
            raise SystemExit(f"ABORT: svpgen300 psam order mismatch at {pfx} — cannot align dosages")
    idx = {}
    for pfx in prefixes:
        with open(f"{pfx}.pvar") as fh:
            vi = 0
            for line in fh:
                if line.startswith("#"):
                    continue
                sid = line.split("\t")[2]
                idx[sid] = (pfx, vi)
                vi += 1
    src_iids = _psam_iids(f"{prefixes[0]}.psam")
    return idx, src_iids


def read_sv_dosages(sv_ids, iids, iidx, sv_index=None, src_iids=None):
    if sv_index is None:
        sv_index, src_iids = _svpgen_index()
    scatter = np.array([iidx.get(i, -1) for i in src_iids], dtype=np.int64)
    keep_src = scatter >= 0
    N = len(iids)
    missing = [s for s in sv_ids if s not in sv_index]
    if missing:
        raise KeyError(f"{len(missing)} SV(s) absent from svpgen300 store, e.g. {missing[:3]}")
    by_chunk = {}
    for k, s in enumerate(sv_ids):
        pfx, vi = sv_index[s]
        by_chunk.setdefault(pfx, []).append((vi, k))
    G = np.empty((len(sv_ids), N), np.float32)
    Nsrc = len(src_iids)
    buf = np.empty(Nsrc, np.float64)
    for pfx, items in by_chunk.items():
        r = pgenlib.PgenReader(f"{pfx}.pgen".encode())
        for vi, k in items:
            r.read_dosages(vi, buf)
            row = np.full(N, np.nan, np.float32)
            row[scatter[keep_src]] = buf[keep_src].astype(np.float32)
            G[k] = row
        r.close()
    return G


def list_window_svs(chrom, sv_index=None):
    if sv_index is None:
        sv_index, _ = _svpgen_index()
    pre = f"{chrom}_"
    return sorted(s for s in sv_index if s.startswith(pre))


EXTRACT_THREADS = int(os.environ.get("FINEMAP_EXTRACT_THREADS", "4"))


def extract_snp_pgen(chrom, keep_fid_iid, out_prefix, bfile=None, force=False):
    if bfile is None:
        bfile = PERI_BFILE.format(chrom=chrom)
    if os.path.exists(out_prefix + ".pgen") and not force:
        return out_prefix
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    with open(out_prefix + ".build.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if os.path.exists(out_prefix + ".pgen") and not force:
            return out_prefix
        tmp = out_prefix + ".tmp"
        cmd = [PLINK2, "--bfile", bfile, "--keep", keep_fid_iid, "--make-pgen",
               "--out", tmp, "--threads", str(EXTRACT_THREADS)]
        log("+ " + " ".join(cmd))
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        for ext in (".pvar", ".psam", ".pgen"):
            os.replace(tmp + ext, out_prefix + ext)
    return out_prefix


def snp_ids_pos(pgen_prefix):
    ids, pos, ref, alt = [], [], [], []
    with open(pgen_prefix + ".pvar") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            ids.append(p[2]); pos.append(int(p[1])); ref.append(p[3]); alt.append(p[4])
    return ids, np.array(pos, np.int64), ref, alt


def iter_snp_blocks(pgen_prefix, iids, iidx, block=2048):
    ids, pos, ref, alt = snp_ids_pos(pgen_prefix)
    src_iids = _psam_iids(pgen_prefix + ".psam")
    scatter = np.array([iidx.get(i, -1) for i in src_iids], dtype=np.int64)
    keep_src = scatter >= 0
    Nsrc = len(src_iids); N = len(iids)
    r = pgenlib.PgenReader((pgen_prefix + ".pgen").encode())
    buf = np.empty(Nsrc, np.float64)
    M = len(ids)
    for b0 in range(0, M, block):
        b1 = min(b0 + block, M)
        D = np.empty((b1 - b0, N), np.float32)
        for j in range(b0, b1):
            r.read_dosages(j, buf)
            row = np.full(N, np.nan, np.float32)
            row[scatter[keep_src]] = buf[keep_src].astype(np.float32)
            D[j - b0] = row
        yield ids[b0:b1], pos[b0:b1], D
    r.close()


def fit_null(X, y, present):
    K = X.shape[1]
    beta = np.zeros(K)
    ybar = float(np.clip(y[present].mean(), 1e-6, 1 - 1e-6))
    beta[0] = np.log(ybar / (1 - ybar))
    conv = False
    for _ in range(25):
        eta = np.clip(X @ beta, -30, 30)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = p * (1 - p); w[~present] = 0.0
        r = y - p; r[~present] = 0.0
        XtWX = X.T @ (X * w[:, None])
        grad = X.T @ r
        try:
            delta = np.linalg.solve(XtWX, grad)
        except np.linalg.LinAlgError:
            delta = np.linalg.lstsq(XtWX, grad, rcond=None)[0]
        beta += delta
        if np.max(np.abs(delta)) < 1e-6:
            conv = True
            break
    eta = np.clip(X @ beta, -30, 30); p = 1.0 / (1.0 + np.exp(-eta))
    w = p * (1 - p); w[~present] = 0.0
    r = y - p; r[~present] = 0.0
    return w, r, conv


def missing_mask(D):
    A = np.asarray(D)
    return ~np.isfinite(A) | (A < -0.0001) | (A > 2.0001)


def masked_meanimpute(D, present):
    Dm = np.array(D, np.float64)
    miss = missing_mask(Dm)
    valid = present[None, :] & ~miss
    cnt = valid.sum(1)
    ssum = np.where(valid, Dm, 0.0).sum(1)
    mu = np.where(cnt > 0, ssum / np.maximum(cnt, 1), 0.0)
    if miss.any():
        mi = np.where(miss)
        Dm[mi] = mu[mi[0]]
    Dm[:, ~present] = 0.0
    npres = max(1, int(present.sum()))
    return Dm, cnt / npres


def std_present(Dm, present):
    S = Dm[:, present].copy()
    S -= S.mean(1, keepdims=True)
    sd = S.std(1, keepdims=True); sd[sd == 0] = 1.0
    return S / sd


class ScoreNull:

    def __init__(self, X, y, min_case=50):
        self.present = y >= 0
        yy = np.where(self.present, y, 0).astype(np.float64)
        self.ncase = int((y == 1).sum()); self.ncontrol = int((y == 0).sum())
        self.tested = self.ncase >= min_case and self.ncontrol >= min_case
        self.X = X
        if not self.tested:
            return
        w, r, conv = fit_null(X, yy, self.present)
        self.converged = conv
        self.degenerate = abs(float(r.sum())) > IDENT_TOL
        self.w = w; self.r = r
        self.XtWXi = np.linalg.pinv(X.T @ (X * w[:, None]))
        self.WX = X * w[:, None]

    def score_pre(self, Dm):
        U = Dm @ self.r
        Vraw = (Dm * Dm) @ self.w
        A = Dm @ self.WX
        corr = np.einsum("bk,kl,bl->b", A, self.XtWXi, A)
        V = Vraw - corr
        good = np.isfinite(V) & (V > VREL * Vraw) & (Vraw > 0)
        z = np.where(good, U / np.sqrt(np.where(good, V, 1.0)), np.nan)
        beta = np.where(good, U / np.where(good, V, 1.0), np.nan)
        from scipy.stats import chi2 as _chi2
        p = np.where(good, _chi2.sf(np.where(good, U * U / np.where(good, V, 1.0), 0.0), 1), np.nan)
        return z.astype(np.float64), beta.astype(np.float64), p.astype(np.float64)

    def score(self, D, return_callrate=False):
        Dm, callrate = masked_meanimpute(D, self.present)
        z, beta, p = self.score_pre(Dm)
        if return_callrate:
            return z, beta, p, callrate
        return z, beta, p


def corr_to_svs(D, Gsv, present):
    def std(M):
        Mm, _ = masked_meanimpute(M, present)
        Ms = Mm[:, present]
        Ms -= Ms.mean(1, keepdims=True)
        sd = Ms.std(1, keepdims=True); sd[sd == 0] = 1.0
        return Ms / sd
    Ds = std(D); Gs = std(Gsv); n = Ds.shape[1]
    C = (Ds @ Gs.T) / n
    return np.max(np.abs(C), axis=1)


def residualize(M, X, present):
    Xp = X[present]
    beta = np.linalg.lstsq(Xp, M[present], rcond=None)[0]
    R = np.zeros_like(M)
    R[present] = M[present] - Xp @ beta
    return R


def _finalize_fit(out_dir, *, chrom, phecode, stratum, mcc, panel_mode, snp_select, sv_ids, hit_sv_ids,
                  Gsv, zsv, psv, null, X, y, present, iids, iidx, wpref, pred_ids, snp_recs, nblk,
                  min_callrate, L_snp_cap, top_assoc, top_tag, t0):
    N = len(iids)
    n_lowcr = sum(1 for r in snp_recs if r[7] < min_callrate)
    log(f"  scored {len(snp_recs)} SNPs in {nblk} blocks; {n_lowcr} below call-rate {min_callrate} "
        f"(dropped from panel) ({time.time()-t0:.0f}s)")

    cand = [r for r in snp_recs if r[7] >= min_callrate]
    by_assoc = sorted([r for r in cand if np.isfinite(r[2])], key=lambda r: -abs(r[2]))[:top_assoc]
    by_tag = sorted([r for r in cand if np.isfinite(r[5])], key=lambda r: -r[5])[:top_tag]
    if snp_select == "assoc_only":
        by_tag = []
    chosen, seen = [], set()
    for r in by_assoc + by_tag:
        if r[0] not in seen:
            seen.add(r[0]); chosen.append(r)
    chosen = chosen[:L_snp_cap]

    win_ids, _, _, _ = snp_ids_pos(wpref)
    idpos = {i: k for k, i in enumerate(win_ids)}
    src_iids = _psam_iids(wpref + ".psam")
    scatter = np.array([iidx.get(i, -1) for i in src_iids], dtype=np.int64)
    keep_src = scatter >= 0
    snp_order = [r[0] for r in chosen]
    Gsnp = np.empty((len(snp_order), N), np.float32)
    r = pgenlib.PgenReader((wpref + ".pgen").encode())
    buf = np.empty(len(src_iids), np.float64)
    for j, sid in enumerate(snp_order):
        r.read_dosages(idpos[sid], buf)
        row = np.full(N, np.nan, np.float32)
        row[scatter[keep_src]] = buf[keep_src].astype(np.float32)
        Gsnp[j] = row
    r.close()

    var_ids = list(sv_ids) + snp_order
    var_type = ["sv"] * len(sv_ids) + ["snp"] * len(snp_order)
    is_pred = [False] * len(sv_ids) + [r[6] for r in chosen]
    zrec = {r[0]: r for r in chosen}
    zvec = list(zsv) + [zrec[i][2] for i in snp_order]
    prec = list(psv) + [zrec[i][4] for i in snp_order]
    Gall = np.vstack([Gsv, Gsnp])
    Gall, _ = masked_meanimpute(Gall, present)
    del Gsv, Gsnp
    Mfull = Gall.T

    yb = np.where(present, y, 0).astype(np.float64)
    Mres = residualize(Mfull, X, present)
    yb_res = np.zeros(N)
    yb_res[present] = yb[present] - X[present] @ np.linalg.lstsq(X[present], yb[present], rcond=None)[0]

    def _corr(Mp):
        Mp -= Mp.mean(0); sd = Mp.std(0); sd[sd == 0] = 1.0; Mp /= sd
        C = (Mp.T @ Mp) / Mp.shape[0]; np.fill_diagonal(C, 1.0); return (C + C.T) / 2
    R = _corr(Mres[present])
    R_raw = _corr(Mfull[present])
    r_raw_res_maxdelta = float(np.max(np.abs(R - R_raw)))

    n_eff = 4.0 / (1.0 / null.ncase + 1.0 / null.ncontrol)

    np.savetxt(f"{out_dir}/R.tsv", R.astype(np.float64), fmt="%.10g", delimiter="\t")
    (Mres if present.all() else Mres[present]).tofile(f"{out_dir}/Mres.bin")
    yb_res[present].astype("<f8").tofile(f"{out_dir}/yres.bin")
    np.save(f"{out_dir}/R.npy", R.astype(np.float64))
    with open(f"{out_dir}/shapes.tsv", "w") as f:
        f.write(f"n_present\t{int(present.sum())}\nn_variants\t{len(var_ids)}\nn_eff\t{n_eff:.6f}\n")
    with open(f"{out_dir}/variant_meta.tsv", "w") as f:
        f.write("idx\tvariant_id\tvariant_type\tis_en_predictor\tz\tp\n")
        for j, vid in enumerate(var_ids):
            f.write(f"{j}\t{vid}\t{var_type[j]}\t{int(is_pred[j])}\t{zvec[j]:.6g}\t{prec[j]:.3g}\n")
    with open(f"{out_dir}/fit_meta.json", "w") as f:
        json.dump({
            "weight_set": "ext2_11129", "svpgen_dir": SVPGEN_DIR,
            "keep_dir": KEEP_DIR, "covar_dir": COVAR_DIR, "pheno_dir": PHENO_DIR,
            "chrom": chrom, "phecode": phecode, "stratum": stratum, "mcc": mcc,
            "panel_mode": panel_mode, "snp_select": snp_select,
            "n_sv": len(sv_ids), "n_snp": len(snp_order),
            "n_variants": len(var_ids), "n_present": int(present.sum()),
            "ncase": null.ncase, "ncontrol": null.ncontrol, "n_eff": n_eff,
            "null_converged": bool(null.converged), "n_snp_scanned": len(snp_recs),
            "n_snp_below_callrate": int(n_lowcr), "min_callrate": min_callrate,
            "r_raw_vs_residualized_maxdelta": r_raw_res_maxdelta,
            "n_en_predictor_in_panel": int(sum(is_pred)),
            "hit_sv_ids": list(hit_sv_ids),
            "elapsed_s": round(time.time() - t0, 1),
        }, f, indent=2)
    log(f"  BUILT {chrom}.{phecode}.{stratum} [{panel_mode}] V={len(var_ids)} "
        f"(sv={len(sv_ids)} snp={len(snp_order)}) n_present={int(present.sum())} "
        f"ncase={null.ncase} ({time.time()-t0:.0f}s) -> {out_dir}")
    return out_dir


def build_fit(chrom, phecode, stratum, mcc, hit_sv_ids, out_dir, panel_mode="hit_only",
              snp_select="assoc_tag", top_assoc=60, top_tag=60, L_snp_cap=200,
              min_callrate=0.90, force=False):
    t0 = time.time()
    os.makedirs(out_dir, exist_ok=True)
    iids, iidx = load_keep(stratum)
    X, covar_names = load_covar(stratum, iidx)
    y = load_pheno(mcc, phecode, iidx)
    N = len(iids)

    sv_index, _ = _svpgen_index()
    if panel_mode == "hit_only":
        sv_ids = list(dict.fromkeys(hit_sv_ids))
    elif panel_mode == "all_window":
        sv_ids = list_window_svs(chrom, sv_index)
    else:
        raise ValueError(panel_mode)
    Gsv = read_sv_dosages(sv_ids, iids, iidx, sv_index=sv_index, src_iids=_psam_iids(_svpgen_prefixes()[0] + ".psam"))

    null = ScoreNull(X, y)
    if not null.tested:
        raise SystemExit(f"SKIP {chrom}.{phecode}.{stratum}: below case/control floor "
                         f"(ncase={null.ncase}, ncontrol={null.ncontrol})")
    if null.degenerate:
        raise SystemExit(f"SKIP {chrom}.{phecode}.{stratum}: degenerate null fit (|sum r|>{IDENT_TOL})")

    present = null.present
    zsv, bsv, psv = null.score(Gsv)

    scratch = os.environ.get("FINEMAP_SCRATCH", "/tmp/finemap_scratch_ext2")
    wpref = f"{scratch}/{stratum}/{chrom}.peri.keep"
    extract_snp_pgen(chrom, f"{KEEP_DIR}/{stratum}.keep.fid_iid.tsv", wpref, force=force)
    pred_ids = load_pred_ids(chrom)

    Gsv_imp, _ = masked_meanimpute(Gsv, present)
    Gsv_std = std_present(Gsv_imp, present)
    npres = int(present.sum())

    blk = max(256, min(2048, int(1.5e9 / (max(N, 1) * 8))))
    snp_recs = []
    nblk = 0
    for ids_b, pos_b, D in iter_snp_blocks(wpref, iids, iidx, block=blk):
        Dm, cr = masked_meanimpute(D, present)
        zb, bb, pb = null.score_pre(Dm)
        Ds = std_present(Dm, present)
        cb = np.max(np.abs((Ds @ Gsv_std.T) / npres), axis=1)
        for k in range(len(ids_b)):
            snp_recs.append((ids_b[k], int(pos_b[k]), float(zb[k]), float(bb[k]),
                             float(pb[k]), float(cb[k]), ids_b[k] in pred_ids, float(cr[k])))
        nblk += 1

    return _finalize_fit(
        out_dir, chrom=chrom, phecode=phecode, stratum=stratum, mcc=mcc, panel_mode=panel_mode,
        snp_select=snp_select, sv_ids=sv_ids, hit_sv_ids=hit_sv_ids, Gsv=Gsv, zsv=zsv, psv=psv,
        null=null, X=X, y=y, present=present, iids=iids, iidx=iidx, wpref=wpref, pred_ids=pred_ids,
        snp_recs=snp_recs, nblk=nblk, min_callrate=min_callrate, L_snp_cap=L_snp_cap,
        top_assoc=top_assoc, top_tag=top_tag, t0=t0)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--phecode", required=True)
    ap.add_argument("--stratum", required=True)
    ap.add_argument("--mcc", required=True, choices=["mcc1", "mcc2"])
    ap.add_argument("--hit-svs", required=True, help="comma-separated sv_locus ids")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--panel-mode", default="hit_only", choices=["hit_only", "all_window"])
    ap.add_argument("--snp-select", default="assoc_tag", choices=["assoc_tag", "assoc_only"])
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    build_fit(a.chrom, a.phecode, a.stratum, a.mcc, a.hit_svs.split(","), a.out_dir,
              panel_mode=a.panel_mode, snp_select=a.snp_select, force=a.force)
