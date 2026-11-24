#!/usr/bin/env python3
"""step30r — PER-ANCESTRY QUANT PheWAS (stress-test chr1->WBC, the last-standing 'real biology' claim)."""
import os, csv, sys, time, gc
import numpy as np
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S
import step30d_numpy_quant as D
import step30q_ancestry_meta as MQ

ANC_DUMMIES = ["anc_amr", "anc_eas", "anc_eur", "anc_mid", "anc_sas"]
DUMMY_TO_ANC = {"anc_amr": "AMR", "anc_eas": "EAS", "anc_eur": "EUR", "anc_mid": "MID", "anc_sas": "SAS"}
ANCS = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
WELL = ["EUR", "AFR", "AMR"]
QMAT = f"{S.OUTA}/quant_matrices"; os.makedirs(QMAT, exist_ok=True)
SVCHUNK = 200


def ancestry_labels(X, names, N):
    lab = np.array(["AFR"] * N, dtype=object)
    for d in ANC_DUMMIES:
        lab[X[:, names.index(d)] == 1.0] = DUMMY_TO_ANC[d]
    return lab


def read_panel(iidx, panel, N):
    hdr = open(f"{S.INP}/quant.pheno").readline().rstrip("\n").split("\t")
    loincs = hdr[2:]; colidx = {lo: j for j, lo in enumerate(loincs)}
    keep = [(lo, colidx[lo]) for lo in panel if lo in colidx]
    P = len(keep); Y = np.full((N, P), np.nan, np.float32)
    with open(f"{S.INP}/quant.pheno") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[1])
            if kk is None:
                continue
            for a, (lo, j) in enumerate(keep):
                v = p[2 + j]
                if v != "NA":
                    Y[kk, a] = float(v)
    return [lo for lo, _ in keep], Y


def scan_anc(anc, G, Xa_full_cols, Y_full, mask, ploincs, sv_order):
    na = int(mask.sum()); M = G.shape[0]; P = len(ploincs)
    Qa, _ = np.linalg.qr(Xa_full_cols[mask].astype(np.float64)); Qa = Qa.astype(np.float32); k = Qa.shape[1]
    Ya = Y_full[mask].copy(); Msk = np.isfinite(Ya).astype(np.float32)
    np.nan_to_num(Ya, copy=False)
    Ya -= Qa @ (Qa.T @ Ya); Ya *= Msk
    ss_y = np.einsum("ij,ij->j", Ya, Ya).astype(np.float64); n_eff = Msk.sum(0).astype(np.float64)
    RyT = np.ascontiguousarray(Ya.T); MkT = np.ascontiguousarray(Msk.T)
    del Ya, Msk; gc.collect()
    beta = np.full((M, P), np.nan, np.float32); pv = np.full((M, P), np.nan, np.float32)
    invv = np.full((M, P), np.nan, np.float32)
    Ga = np.ascontiguousarray(G[:, mask])
    for si in range(0, M, SVCHUNK):
        blk = Ga[si:si + SVCHUNK]
        blk -= (blk @ Qa) @ Qa.T
    for si in range(0, M, SVCHUNK):
        Rs = Ga[si:si + SVCHUNK].astype(np.float32)
        num = (Rs @ RyT.T).astype(np.float64)
        denx = ((Rs * Rs) @ MkT.T).astype(np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            rcorr = num / np.sqrt(denx * ss_y[None, :]); b = num / denx
        df = n_eff[None, :] - k - 1
        good = (denx > 0) & np.isfinite(rcorr) & (np.abs(rcorr) < 1) & (df > 0)
        t = np.where(good, rcorr * np.sqrt(np.maximum(df, 1) / np.maximum(1 - rcorr ** 2, 1e-12)), np.nan)
        p = 2 * norm.sf(np.abs(t))
        iv = np.where(good, df * denx / np.maximum(ss_y[None, :] * (1 - rcorr ** 2), 1e-30), np.nan)
        sl = slice(si, si + Rs.shape[0])
        beta[sl] = np.where(good, b, np.nan).astype(np.float32)
        pv[sl] = p.astype(np.float32); invv[sl] = iv.astype(np.float32)
    del Ga, RyT, MkT; gc.collect()
    for nm, arr in (("beta", beta), ("invvar", invv), ("p", pv)):
        np.save(f"{QMAT}/quant_{nm}_{anc}.npy.tmp.npy", arr)
        os.replace(f"{QMAT}/quant_{nm}_{anc}.npy.tmp.npy", f"{QMAT}/quant_{nm}_{anc}.npy")
    with open(f"{QMAT}/quant_loinc_{anc}.txt", "w") as f:
        f.write("\n".join(ploincs) + "\n")
    with open(f"{QMAT}/quant_svorder_{anc}.txt", "w") as f:
        f.write("\n".join(sv_order) + "\n")
    nsig = int(np.nansum((pv < 0.05) & (np.abs(beta) < 3)))
    print(f"  [{anc}] n={na}, panel={P}, sane sig(p<0.05,|b|<3)={nsig}", flush=True)
    return na, nsig


def main():
    t0 = time.time()
    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    keep_cols = [j for j, c in enumerate(names) if c not in ANC_DUMMIES]
    Xa_cols = np.ascontiguousarray(X[:, keep_cols])
    lab = ancestry_labels(X, names, N)
    print("ancestry N:", {a: int((lab.astype(str) == a).sum()) for a in ANCS}, flush=True)
    sv_order = []; G = S.load_G(sv_order); M = len(sv_order)
    panel = D.select_panel(); ploincs, Y_full = read_panel(iidx, panel, N)
    print(f"panel={len(ploincs)} LOINCs; Y_full {Y_full.shape}", flush=True)

    summ = []
    for anc in ANCS:
        mask = (lab.astype(str) == anc); na = int(mask.sum())
        if na < 400:
            print(f"[{anc}] n={na} <400 skip", flush=True); summ.append((anc, na, 0)); continue
        need = [f"{QMAT}/quant_{x}_{anc}.npy" for x in ("beta", "invvar", "p")] + \
               [f"{QMAT}/quant_loinc_{anc}.txt", f"{QMAT}/quant_svorder_{anc}.txt"]
        if all(os.path.exists(f) for f in need):
            print(f"[{anc}] complete — skip", flush=True); summ.append((anc, na, -1)); continue
        print(f"=== {anc}: n={na} ({time.time()-t0:.0f}s) ===", flush=True)
        _, ns = scan_anc(anc, G, Xa_cols, Y_full, mask, ploincs, sv_order); summ.append((anc, na, ns))
        gc.collect()
    print("\nPER-ANC QUANT SUMMARY:", summ, flush=True)

    def load(anc):
        if not os.path.exists(f"{QMAT}/quant_beta_{anc}.npy"):
            return None
        b = np.load(f"{QMAT}/quant_beta_{anc}.npy"); iv = np.load(f"{QMAT}/quant_invvar_{anc}.npy")
        lo = [l.strip() for l in open(f"{QMAT}/quant_loinc_{anc}.txt") if l.strip()]
        sv = [l.strip() for l in open(f"{QMAT}/quant_svorder_{anc}.txt") if l.strip()]
        return dict(b=b, iv=iv, loc={l: j for j, l in enumerate(lo)}, svr={s: i for i, s in enumerate(sv)})
    A = {a: load(a) for a in ANCS}; have = [a for a in ANCS if A[a]]
    WBC = ["6690-2", "770-8", "751-8"]
    KNOWN = {"29463-7": "weight", "8302-2": "height", "39156-5": "BMI", "32623-1": "MPV", "776-5": "PDW"}
    out = open(f"{S.OUTA}/quant_ancestry_meta.tsv", "w")
    out.write("sv_locus\tloinc\t" + "\t".join(f"{a}_beta\t{a}_p" for a in WELL) +
              "\tmeta_beta\tmeta_p\tI2\tk\tverdict\n")
    highlight = []
    ref = A.get("EUR") or A[have[0]]
    svs = list(ref["svr"].keys()); locs = list(ref["loc"].keys())
    for sv in svs:
        for lo in locs:
            cells = {}
            for a in have:
                d = A[a]; i = d["svr"].get(sv); j = d["loc"].get(lo)
                if i is None or j is None:
                    continue
                bv = float(d["b"][i, j]); iv = float(d["iv"][i, j])
                if np.isfinite(bv) and np.isfinite(iv) and iv > 0:
                    cells[a] = (bv * iv, iv, 0)
            if len([a for a in cells if a in WELL]) < 1:
                continue
            res = MQ.analyze(cells)
            if not np.isfinite(res["p"]):
                continue
            row = [sv, lo]
            for a in WELL:
                st = res["strata"].get(a)
                row += [f"{st['beta']:.4g}" if st else "", f"{st['p']:.3g}" if st else ""]
            row += [f"{res['beta']:.4g}", f"{res['p']:.3g}", f"{res['I2']:.0f}", str(res["k"]), res["verdict"]]
            out.write("\t".join(row) + "\n")
            if lo in WBC or lo in KNOWN:
                highlight.append((sv, lo, res))
    out.close()
    print(f"\nwrote {S.OUTA}/quant_ancestry_meta.tsv", flush=True)

    from collections import Counter
    print("\n=== chr1 -> WBC/neutrophil (the 'real biology' test) ===", flush=True)
    print(f"  {'SV':<40}{'loinc':<8}{'EUR_b/p':>16}{'AFR_b/p':>16}{'AMR_b/p':>16}{'metaP':>9} verdict")
    for sv, lo, res in sorted([h for h in highlight if h[1] in WBC and h[0].startswith("chr1_")],
                              key=lambda x: x[2]["p"] if np.isfinite(x[2]["p"]) else 9):
        s = f"  {sv[:39]:<40}{lo:<8}"
        for a in WELL:
            st = res["strata"].get(a)
            cell = f"{st['beta']:+.3f}/{st['p']:.0e}" if st else "na"
            s += f"{cell:>16}"
        s += f"{res['p']:>9.0e} {res['verdict']}"
        print(s, flush=True)
    print("\n=== known quant signals (weight/height/BMI/MPV/PDW) ===", flush=True)
    for sv, lo, res in sorted([h for h in highlight if h[1] in KNOWN], key=lambda x: x[2]["p"]):
        if not np.isfinite(res["p"]) or res["p"] > 1e-4:
            continue
        s = f"  {sv[:36]:<38}{KNOWN[lo]:<8}"
        for a in WELL:
            st = res["strata"].get(a)
            cell = f"{st['beta']:+.3f}/{st['p']:.0e}" if st else "na"
            s += f"{cell:>16}"
        s += f"{res['p']:>9.0e} I2={res['I2']:.0f} {res['verdict']}"
        print(s, flush=True)
    print("QUANT_ANCESTRY_DONE", flush=True)


if __name__ == "__main__":
    main()
