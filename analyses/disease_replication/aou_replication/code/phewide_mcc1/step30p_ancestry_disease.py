#!/usr/bin/env python3
"""step30p — PER-ANCESTRY disease PheWAS (the decisive artifact-vs-biology test for the pooled hits)."""
import os, csv, sys, time, gc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S

ANC_DUMMIES = ["anc_amr", "anc_eas", "anc_eur", "anc_mid", "anc_sas"]
DUMMY_TO_ANC = {"anc_amr": "AMR", "anc_eas": "EAS", "anc_eur": "EUR", "anc_mid": "MID", "anc_sas": "SAS"}
ANCS = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
WELL = {"EUR", "AFR", "AMR"}
BATCH = 64
FLUSH_EVERY = 6


def ancestry_labels(X, names, N):
    lab = np.array(["AFR"] * N, dtype=object)
    for d in ANC_DUMMIES:
        col = names.index(d)
        lab[X[:, col] == 1.0] = DUMMY_TO_ANC[d]
    return lab


def _save_npy(path, arr):
    np.save(path + ".tmp.npy", arr); os.replace(path + ".tmp.npy", path)


def scan_subset(Ga, G2a, Xa, Ya, phecodes, tiers, label, sv_order, out_tag, min_case, t0):
    M = Ga.shape[0]; J = len(phecodes); na = Ga.shape[1]
    Pv = np.full((M, J), np.nan, np.float32); Be = np.full((M, J), np.nan, np.float32)
    Um = np.full((M, J), np.nan, np.float32); Vm = np.full((M, J), np.nan, np.float32)
    ncase = np.zeros(J, np.int32); npres = np.zeros(J, np.int32); tested = np.zeros(J, bool)
    n_nonconv = 0; n_degen = 0

    pv_path = f"{S.MATDIR}/disease_Pv{out_tag}.npy"; ord_path = f"{S.MATDIR}/disease_svorder{out_tag}.txt"
    meta_path = f"{S.MATDIR}/disease_meta{out_tag}.tsv"; u_path = f"{S.MATDIR}/disease_U{out_tag}.npy"
    if all(os.path.exists(p) for p in (pv_path, ord_path, meta_path, u_path)):
        prev_order = [ln.split("\t")[0] for ln in open(ord_path) if ln.strip()]
        prevPv = np.load(pv_path)
        if prev_order == sv_order and prevPv.shape == (M, J):
            Pv[:] = prevPv; Be[:] = np.load(f"{S.MATDIR}/disease_Be{out_tag}.npy")
            Um[:] = np.load(u_path); Vm[:] = np.load(f"{S.MATDIR}/disease_V{out_tag}.npy")
            for b, r in enumerate(csv.DictReader(open(meta_path), delimiter="\t")):
                tested[b] = r["tested"] == "1"; ncase[b] = int(r["cases"]); npres[b] = int(r["n"])
            print(f"    [{out_tag}] RESUME: {int(tested.sum())}/{J} cols already done — skipping", flush=True)

    def flush(progress):
        S.save_matrices(Pv, Be, phecodes, label, ncase, npres, tested, sv_order, tiers, out_tag, progress)
        _save_npy(f"{S.MATDIR}/disease_U{out_tag}.npy", Um)
        _save_npy(f"{S.MATDIR}/disease_V{out_tag}.npy", Vm)

    nbatch = 0
    for j0 in range(0, J, BATCH):
        cols = [j for j in range(j0, min(j0 + BATCH, J)) if not tested[j]]
        if not cols:
            continue
        Be_b, Pv_b, nca_b, npr_b, tt_b, nc_b, dg_b, U_b, V_b = S.score_test_batch(
            Ga, G2a, Xa, Ya[:, cols], min_case=min_case, return_uv=True)
        n_nonconv += int(nc_b.sum()); n_degen += int(dg_b.sum())
        for a, j in enumerate(cols):
            ncase[j] = nca_b[a]; npres[j] = npr_b[a]
            if tt_b[a]:
                Be[:, j] = Be_b[:, a]; Pv[:, j] = Pv_b[:, a]
                Um[:, j] = U_b[:, a]; Vm[:, j] = V_b[:, a]; tested[j] = True
        nbatch += 1
        if nbatch % 10 == 0:
            print(f"    [{out_tag}] {min(j0+BATCH,J)}/{J} scanned, {int(tested.sum())} tested "
                  f"({time.time()-t0:.0f}s)", flush=True)
        if nbatch % FLUSH_EVERY == 0:
            flush(f"in_progress {out_tag}: {min(j0+BATCH,J)}/{J} cols ({time.time()-t0:.0f}s)")

    sv_of = np.repeat(np.arange(M), J); ph_of = np.tile(np.arange(J), M)
    Pv_flat = Pv.ravel(); Be_flat = Be.ravel()
    valid = np.isfinite(Pv_flat) & tested[ph_of]
    q = np.full(M * J, np.nan)
    for tn in ("signal", "null_ext"):
        idx = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if idx.size:
            q[idx] = S.bh(Pv_flat[idx])
    n_q05 = int(np.nansum(q < 0.05))
    keep = np.nonzero(valid & ((np.nan_to_num(q, nan=1.0) < 0.05) | (Pv_flat < 1e-4)))[0]
    keep = keep[np.argsort(Pv_flat[keep])]
    out_c = f"{S.OUTA}/phewide_disease_full{out_tag}.tsv"
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
    t0 = time.time()
    S.MATDIR = os.path.join(S.MATDIR, "per_ancestry"); os.makedirs(S.MATDIR, exist_ok=True)
    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    keep_cols = [j for j, c in enumerate(names) if c not in ANC_DUMMIES]
    print(f"covar cols kept per-ancestry ({len(keep_cols)}): {[names[j] for j in keep_cols]}", flush=True)
    lab = ancestry_labels(X, names, N)
    uniq, cnt = np.unique(lab.astype(str), return_counts=True)
    print("ancestry N:", dict(zip(uniq.tolist(), cnt.tolist())), flush=True)

    sv_order = []; G = S.load_G(sv_order); M = len(sv_order)
    phecodes = []; Y = S.load_disease_matrix(iidx, phecodes); J = len(phecodes)
    svtier_path = next((p for p in (f"{S.SC}/sv_tier.tsv", f"{S.INP}/sv_tier.tsv") if os.path.exists(p)), None)
    if svtier_path is None:
        raise SystemExit(f"ABORT: sv_tier.tsv not found at {S.SC}/ or {S.INP}/ — cannot assign SV tiers.")
    tier = {r["sv_locus"]: r["tier"] for r in csv.DictReader(open(svtier_path), delimiter="\t")}
    label = {r["phecode"]: r["label"] for r in csv.DictReader(open(S.PSUM), delimiter="\t")}
    tiers = np.array([tier.get(s, "") for s in sv_order])

    summary = []
    for anc in ANCS:
        out_tag = f"_{anc}"
        donef = f"{S.MATDIR}/disease_progress{out_tag}.txt"
        if os.path.exists(donef) and open(donef).read().startswith("DONE"):
            print(f"[{anc}] already DONE — skipping", flush=True); continue
        mask = (lab.astype(str) == anc)
        na = int(mask.sum())
        min_case = 100 if na < 10000 else 50
        print(f"\n=== {anc}: n={na}  (min_case={min_case}{', CONFIRMATION-ONLY' if anc not in WELL else ''}) ===",
              flush=True)
        if na < 400:
            print(f"  [{anc}] n<400 — skipping (nothing testable at min_case=100)", flush=True)
            summary.append((anc, na, 0, 0, 0)); continue
        Ga = np.ascontiguousarray(G[:, mask]); G2a = Ga * Ga
        Xa = np.ascontiguousarray(X[np.ix_(mask, keep_cols)])
        Ya = np.ascontiguousarray(Y[mask])
        nt, nq, nk = scan_subset(Ga, G2a, Xa, Ya, phecodes, tiers, label, sv_order, out_tag, min_case, t0)
        summary.append((anc, na, nt, nq, nk))
        del Ga, G2a, Xa, Ya; gc.collect()

    print("\n=== PER-ANCESTRY SUMMARY ===", flush=True)
    print(f"{'anc':<6}{'N':>8}{'tested':>8}{'q<0.05':>8}{'written':>8}")
    for anc, na, nt, nq, nk in summary:
        print(f"{anc:<6}{na:>8}{nt:>8}{nq:>8}{nk:>8}")
    print(f"ALL_ANC_DONE ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
