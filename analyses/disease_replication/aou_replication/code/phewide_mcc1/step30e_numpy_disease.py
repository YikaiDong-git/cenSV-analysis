#!/usr/bin/env python3
import os, csv, importlib.util
import numpy as np
from scipy.stats import norm

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
SC = "/tmp/phewide"
INP = f"{SC}/inputs"
RV = f"{BASE}/aou_replication/results_vds"
OUTA = f"{RV}/step30_phewide"
XW = f"{RV}/step31_phewide/trait_crosswalk.tsv"
DSUM = f"{OUTA}/phecodeX_pheno_summary.tsv"
RSV = f"{SC}/R_sv.f32"
SVCHUNK = 100; TOP_N = 60

_s = importlib.util.spec_from_file_location("s30d", os.path.join(os.path.dirname(__file__), "step30d_numpy_quant.py"))
d = importlib.util.module_from_spec(_s); _s.loader.exec_module(d)


def select_panel():
    mapped = {r["aou_id"] for r in csv.DictReader(open(XW), delimiter="\t")
              if r["aou_type"] == "phecode" and r["aou_id"]}
    top = [(r["phecode"], int(r["cases_mcc2"])) for r in csv.DictReader(open(DSUM), delimiter="\t")]
    top.sort(key=lambda x: -x[1])
    return list(dict.fromkeys([p for p, _ in top if p in mapped] + [p for p, _ in top[:TOP_N]]))


def build_panel_ry(iids, iidx, Q, panel):
    hdr = open(f"{INP}/disease.pheno").readline().rstrip("\n").split("\t")
    phs = hdr[2:]; colidx = {p: j for j, p in enumerate(phs)}
    keep = [(p, colidx[p]) for p in panel if p in colidx]
    P = len(keep); N = len(iids)
    Y = np.full((N, P), np.nan, np.float32)
    with open(f"{INP}/disease.pheno") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[0])
            if kk is None:
                continue
            for a, (ph, j) in enumerate(keep):
                v = p[2 + j]
                if v != "NA":
                    Y[kk, a] = float(v)
    import gc
    Msk = np.isfinite(Y).astype(np.float32)
    np.nan_to_num(Y, copy=False)
    Y -= Q @ (Q.T @ Y); Y *= Msk
    ss_y = np.einsum("ij,ij->j", Y, Y).astype(np.float64); n_eff = Msk.sum(0).astype(np.float64)
    RyT = np.ascontiguousarray(Y.T); MkT = np.ascontiguousarray(Msk.T)
    del Y, Msk; gc.collect()
    print(f"disease panel residuals [{P} x {N}] in RAM", flush=True)
    return [p for p, _ in keep], RyT, MkT, ss_y, n_eff


def main():
    iids, iidx = d.load_order(); N = len(iids)
    Q = d.build_Q(iids, iidx); k = Q.shape[1]
    sv_order = d.build_Rsv(iids, Q); M = len(sv_order)
    tier = {r["sv_locus"]: r["tier"] for r in csv.DictReader(open(f"{SC}/phewide_pgen/sv_tier.tsv"), delimiter="\t")}
    pphs, RyT, MkT, ss_y, n_eff = build_panel_ry(iids, iidx, Q, select_panel())
    P = len(pphs)
    R = np.memmap(RSV, dtype=np.float32, mode="r", shape=(M, N))
    svs = np.array(sv_order); tiers = np.array([tier.get(s, "") for s in sv_order])
    T = M * P
    Pv = np.full(T, np.nan, np.float64); Be = np.zeros(T, np.float32); Rc = np.zeros(T, np.float32)
    for si in range(0, M, SVCHUNK):
        Rs = np.array(R[si:si + SVCHUNK], np.float32)
        num = (Rs @ RyT.T).astype(np.float64); denx = ((Rs * Rs) @ MkT.T).astype(np.float64)
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
            q[idx] = d.bh(Pv[idx])
    lab = {r["phecode"]: r["label"] for r in csv.DictReader(open(DSUM), delimiter="\t")}
    with open(f"{SC}/phewide_disease_allpanel.tsv", "w") as f:
        f.write("sv_locus\ttier\tphecode\tn\tbeta\tp\tq_bh_tier\n")
        for i in np.nonzero(valid)[0]:
            a = i // P; b = i % P
            f.write(f"{svs[a]}\t{tiers[a]}\t{pphs[b]}\t{int(n_eff[b])}\t{Be[i]:.4g}\t{Pv[i]:.3g}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    out = np.nonzero(valid & ((Pv < 1e-3) | (np.nan_to_num(q, nan=1.0) < 0.1)))[0]
    out = out[np.argsort(Pv[out])]
    with open(f"{OUTA}/phewide_disease.tsv", "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tn\tbeta_lpm\tpartial_r\tp\tq_bh_tier\n")
        for i in out:
            a = i // P; b = i % P
            f.write(f"{svs[a]}\t{tiers[a]}\t{pphs[b]}\t{lab.get(pphs[b],'')}\t{int(n_eff[b])}\t"
                    f"{Be[i]:.4g}\t{Rc[i]:.4g}\t{Pv[i]:.3g}\t{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    sig = int(np.nansum(q < 0.05))
    print(f"DISEASE_DONE: panel={P} phecodes, {int(valid.sum())} valid, {sig} q<0.05 (tiered), "
          f"{len(out)} written (LPM screen) -> {OUTA}/phewide_disease.tsv", flush=True)


if __name__ == "__main__":
    main()
