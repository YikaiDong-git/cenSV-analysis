#!/usr/bin/env python3
"""step30h4 — PER-ANCESTRY Firth recalibration of the disease PheWAS (closes the per-ancestry Firth gap)."""
import os, sys, csv, time, gc
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S
import step30p_ancestry_disease as P
import step30h_firth_refit_hits as H
import step30h2_approx_firth_recalibrate as H2

SCORE_MATDIR = f"{S.OUTA}/disease_matrices/per_ancestry"
FIRTH_MATDIR = f"{S.OUTA}/disease_matrices_firth/per_ancestry"
os.makedirs(FIRTH_MATDIR, exist_ok=True)
ORDER = ["MID", "SAS", "EAS", "AMR", "AFR", "EUR"]
NWORK = 14


def exact_refit(hitrows, svidx, phidx, Ga, Xa, Ya):
    loci = sorted({r["sv_locus"] for r in hitrows if r["sv_locus"] in svidx})
    phes = sorted({r["phecode"] for r in hitrows if r["phecode"] in phidx})
    if not loci or not phes:
        return {}
    H._DOS = {sid: np.ascontiguousarray(Ga[svidx[sid]]) for sid in loci}
    H._YV = {ph: np.ascontiguousarray(Ya[:, phidx[ph]]) for ph in phes}
    H._XCOV = np.ascontiguousarray(Xa[:, 1:])
    tasks = []; keys = []
    for r in hitrows:
        if r["sv_locus"] in H._DOS and r["phecode"] in H._YV:
            b0 = float(r.get("beta_logOR") or 0.0)
            tasks.append((len(keys), r["sv_locus"], r["phecode"], b0)); keys.append((r["sv_locus"], r["phecode"]))
    out = {}
    if not tasks:
        return out
    with Pool(NWORK) as pool:
        for idx, res in pool.imap_unordered(H._refit_one, tasks, chunksize=16):
            bf, pf, conv, _conf = res
            out[keys[idx]] = (bf, pf, bool(conv))
    H._DOS = H._YV = H._XCOV = None
    return out


def approx_band(Ga, Xa, Ya, Pv_anc, tested_idx):
    H2._G, H2._X, H2._Y, H2._Pv = Ga, Xa, Ya, Pv_anc
    cells = {}
    with Pool(H2.NWORK) as pool:
        for res in pool.imap_unordered(H2.process_phecode, tested_idx, chunksize=8):
            if res is None:
                continue
            j, recs = res
            for (m_i, beta_k, pfir_k, cvg_k, sep_k, fa, fc, forr, fp) in recs:
                cells[(int(m_i), int(j))] = (float(pfir_k), float(beta_k), int(cvg_k), int(sep_k),
                                             float(fp) if fp == fp else np.nan)
    H2._G = H2._X = H2._Y = H2._Pv = None
    return cells


def run_stratum(anc, G, X, Y, lab, keep_cols, sv_order, tiers, phecodes, svidx, phidx, t0):
    final_pv = f"{FIRTH_MATDIR}/disease_Pv_final_{anc}.npy"
    if os.path.exists(final_pv):
        print(f"[{anc}] disease_Pv_final_{anc}.npy exists — skipping", flush=True); return None
    pva = f"{SCORE_MATDIR}/disease_Pv_{anc}.npy"
    if not os.path.exists(pva):
        print(f"[{anc}] no score matrix ({pva}) — skipping", flush=True); return None
    mask = (lab.astype(str) == anc); na = int(mask.sum())
    print(f"\n=== {anc}: n={na}{'  CONFIRMATION-ONLY' if anc not in P.WELL else ''} ===", flush=True)
    if na < 400:
        print(f"  [{anc}] n<400 — skipping", flush=True); return (anc, na, 0, 0, 0, 0, 0)

    Pv_anc = np.load(pva).astype(np.float64)
    Be_anc = np.load(f"{SCORE_MATDIR}/disease_Be_{anc}.npy").astype(np.float64)
    meta = list(csv.DictReader(open(f"{SCORE_MATDIR}/disease_meta_{anc}.tsv"), delimiter="\t"))
    M, J = Pv_anc.shape
    assert len(meta) == J and [m["phecode"] for m in meta] == phecodes, f"{anc}: meta/phecode misalignment"
    assert (M, J) == (len(sv_order), len(phecodes)), f"{anc}: matrix shape {Pv_anc.shape} != ({len(sv_order)},{len(phecodes)})"
    svo_anc = [ln.split("\t")[0] for ln in open(f"{SCORE_MATDIR}/disease_svorder_{anc}.txt") if ln.strip()]
    assert svo_anc == sv_order, f"{anc}: svorder mismatch vs loaded G — matrix rows would be misaligned"
    ph_label = {m["phecode"]: m.get("label", "") for m in meta}
    ph_cases = {m["phecode"]: int(m["cases"]) for m in meta}
    ph_n = {m["phecode"]: int(m["n"]) for m in meta}
    tested_arr = np.array([m["tested"] == "1" for m in meta])
    tested_idx = [j for j in range(J) if tested_arr[j]]

    Ga = np.ascontiguousarray(G[:, mask])
    Ya = np.ascontiguousarray(Y[mask])
    Xa = np.ascontiguousarray(X[np.ix_(mask, keep_cols)])
    print(f"  [{anc}] subset: Ga{Ga.shape} Xa{Xa.shape} ({Xa.shape[1]} covars, dummies dropped), "
          f"{len(tested_idx)} tested phecodes ({time.time()-t0:.0f}s)", flush=True)

    hits_path = f"{S.OUTA}/phewide_disease_full_{anc}.tsv"
    hitrows = list(csv.DictReader(open(hits_path), delimiter="\t")) if os.path.exists(hits_path) else []
    exact = exact_refit(hitrows, svidx, phidx, Ga, Xa, Ya)
    print(f"  [{anc}] exact Firth on {len(exact)} score-hit cells "
          f"({sum(1 for v in exact.values() if v[2])} converged) ({time.time()-t0:.0f}s)", flush=True)

    approx = approx_band(Ga, Xa, Ya, Pv_anc, tested_idx)
    print(f"  [{anc}] approx Firth band = {len(approx)} cells "
          f"({sum(1 for v in approx.values() if v[3])} separation-flagged) ({time.time()-t0:.0f}s)", flush=True)

    Pf = Pv_anc.copy(); Bf = Be_anc.copy(); source = {}; sepmap = {}
    n_ex = n_ap = 0
    for (sv, ph), (bf, pf, conv) in exact.items():
        i = svidx.get(sv); j = phidx.get(ph)
        if i is None or j is None:
            continue
        if conv and np.isfinite(pf):
            Pf[i, j] = pf
            if np.isfinite(bf):
                Bf[i, j] = bf
            source[(i, j)] = "exact"; n_ex += 1
        elif (i, j) in approx:
            pa, ba, _cv, sep, fp = approx[(i, j)]
            Pf[i, j] = pa; Bf[i, j] = ba; source[(i, j)] = "approx"; sepmap[(i, j)] = (sep, fp); n_ap += 1
    for (i, j), (pa, ba, _cv, sep, fp) in approx.items():
        if (i, j) in source:
            continue
        Pf[i, j] = pa; Bf[i, j] = ba; source[(i, j)] = "approx"
        if sep:
            sepmap[(i, j)] = (sep, fp)
    print(f"  [{anc}] merged {n_ex} exact + {n_ap} approx-fallback (+ {len(source)-n_ex-n_ap} approx-only band)",
          flush=True)

    sv_of = np.repeat(np.arange(M), J); ph_of = np.tile(np.arange(J), M)
    flat = Pf.ravel(); Bflat = Bf.ravel()
    valid = np.isfinite(flat) & tested_arr[ph_of]
    q = np.full(M * J, np.nan)
    for tn in ("signal", "null_ext"):
        ix = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if ix.size:
            q[ix] = S.bh(flat[ix])

    np.save(f"{FIRTH_MATDIR}/disease_Pv_final_{anc}.npy", Pf.astype(np.float32))
    np.save(f"{FIRTH_MATDIR}/disease_Be_final_{anc}.npy", Bf.astype(np.float32))

    with open(f"{FIRTH_MATDIR}/disease_firth_cells_{anc}.tsv", "w") as f:
        f.write("sv_locus\ttier\tphecode\tp_final\tp_source\tbeta_final_logOR\tseparation\tfisher_p\n")
        for (i, j), src in sorted(source.items(), key=lambda kv: flat[kv[0][0] * J + kv[0][1]]):
            sep, fp = sepmap.get((i, j), (0, np.nan))
            f.write(f"{sv_order[i]}\t{tiers[i]}\t{phecodes[j]}\t{Pf[i,j]:.3g}\t{src}\t{Bf[i,j]:.4g}\t"
                    f"{sep}\t{'' if not np.isfinite(fp) else f'{fp:.3g}'}\n")

    keep = np.nonzero(valid & ((np.nan_to_num(q, nan=1.0) < 0.10) | (flat < 1e-4)))[0]
    keep = keep[np.argsort(flat[keep])]
    n_sig = n_ind = n_sep = 0
    with open(f"{FIRTH_MATDIR}/disease_hits_final_{anc}.tsv", "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tbeta_final_logOR\tp_final\tp_source\t"
                "q_bh_tier\tseparation\tfisher_p\ttier_label\n")
        for idx in keep:
            i = idx // J; j = idx % J; qi = q[idx]
            src = source.get((i, j), "score")
            sep, fp = sepmap.get((i, j), (0, np.nan))
            sig = bool(np.isfinite(qi) and qi < 0.05)
            base = ("significant" if sig else "indicative" if (np.isfinite(qi) and qi < 0.10) else "sub-threshold")
            lab = "indicative (separation)" if sep else base
            if sep and (sig or base == "indicative"):
                n_sep += 1
            elif sig:
                n_sig += 1
            elif base == "indicative":
                n_ind += 1
            ph = phecodes[j]
            f.write(f"{sv_order[i]}\t{tiers[i]}\t{ph}\t{ph_label.get(ph,'')}\t{ph_cases.get(ph,0)}\t{ph_n.get(ph,0)}\t"
                    f"{Bflat[idx]:.4g}\t{flat[idx]:.3g}\t{src}\t{'' if not np.isfinite(qi) else f'{qi:.4g}'}\t"
                    f"{sep}\t{'' if not np.isfinite(fp) else f'{fp:.3g}'}\t{lab}\n")
    print(f"  [{anc}] FINAL: {n_sig} significant / {n_ind} indicative / {n_sep} separation-flagged "
          f"-> disease_{{Pv,Be}}_final_{anc}.npy ({time.time()-t0:.0f}s)", flush=True)
    del Ga, Ya, Xa, Pv_anc, Be_anc; gc.collect()
    return (anc, na, len(exact), n_sig, n_ind, n_sep, len(approx))


def main():
    t0 = time.time()
    only = None
    if "--only" in sys.argv:
        only = set(sys.argv[sys.argv.index("--only") + 1].split(","))
    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    keep_cols = [j for j, c in enumerate(names) if c not in P.ANC_DUMMIES]
    print(f"per-ancestry covars ({len(keep_cols)}): {[names[j] for j in keep_cols]}", flush=True)
    lab = P.ancestry_labels(X, names, N)
    uniq, cnt = np.unique(lab.astype(str), return_counts=True)
    print("ancestry N:", dict(zip(uniq.tolist(), cnt.tolist())), flush=True)

    sv_order = []; G = S.load_G(sv_order); M = len(sv_order)
    phecodes = []; Y = S.load_disease_matrix(iidx, phecodes); J = len(phecodes)
    svidx = {s: i for i, s in enumerate(sv_order)}; phidx = {p: j for j, p in enumerate(phecodes)}
    svtier_path = next((p for p in (f"{S.SC}/sv_tier.tsv", f"{S.INP}/sv_tier.tsv") if os.path.exists(p)), None)
    if svtier_path is None:
        raise SystemExit(f"ABORT: sv_tier.tsv not found at {S.SC}/ or {S.INP}/")
    tier = {r["sv_locus"]: r["tier"] for r in csv.DictReader(open(svtier_path), delimiter="\t")}
    tiers = np.array([tier.get(s, "") for s in sv_order])
    print(f"loaded G{G.shape} Y{Y.shape}, {M} SVs x {J} phecodes ({time.time()-t0:.0f}s)", flush=True)

    summary = []
    for anc in ORDER:
        if only and anc not in only:
            continue
        r = run_stratum(anc, G, X, Y, lab, keep_cols, sv_order, tiers, phecodes, svidx, phidx, t0)
        if r is not None:
            summary.append(r)

    print("\n=== PER-ANCESTRY FIRTH SUMMARY ===", flush=True)
    print(f"{'anc':<6}{'N':>9}{'exact':>8}{'sig':>6}{'ind':>6}{'sep':>6}{'band':>7}")
    for (anc, na, nex, nsig, nind, nsep, nband) in summary:
        print(f"{anc:<6}{na:>9}{nex:>8}{nsig:>6}{nind:>6}{nsep:>6}{nband:>7}")
    print(f"PERANC_FIRTH_DONE ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
