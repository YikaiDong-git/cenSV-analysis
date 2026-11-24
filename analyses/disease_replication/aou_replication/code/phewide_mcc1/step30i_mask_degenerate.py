#!/usr/bin/env python3
import os, sys, csv, time, shutil
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S

INP = S.INP; OUTA = S.OUTA; MATDIR = S.MATDIR; IDENT_TOL = S.IDENT_TOL
HITS = f"{OUTA}/phewide_disease_full.tsv"


FLOOD_FRAC = 0.02
FLOOD_P = 1e-6


def main():
    t0 = time.time()
    meta = list(csv.DictReader(open(f"{MATDIR}/disease_meta.tsv"), delimiter="\t"))
    phecodes = [m["phecode"] for m in meta]; J = len(phecodes)
    tested0 = np.array([m["tested"] == "1" for m in meta])

    Pv = np.load(f"{MATDIR}/disease_Pv.npy"); Be = np.load(f"{MATDIR}/disease_Be.npy")
    M = Pv.shape[0]
    with np.errstate(invalid="ignore"):
        frac_flood = np.array([(np.nan_to_num(Pv[:, j], nan=1.0) < FLOOD_P).mean() for j in range(J)])
    degen = tested0 & (frac_flood > FLOOD_FRAC)
    print(f"DEGENERATE phecodes masked: {int(degen.sum())} -> "
          f"{[(phecodes[j], round(float(frac_flood[j]),3)) for j in np.where(degen)[0]]}", flush=True)
    print(f"  (next-highest flood fraction among kept: {np.sort(frac_flood[tested0 & ~degen])[-1]:.4g})", flush=True)
    svrows = [ln.rstrip("\n").split("\t") for ln in open(f"{MATDIR}/disease_svorder.txt") if ln.strip()]
    sv_order = [r[0] for r in svrows]; tiers = np.array([r[1] if len(r) > 1 else "" for r in svrows]); M = len(sv_order)
    label = {m["phecode"]: m["label"] for m in meta}
    ncase = np.array([int(m["cases"]) for m in meta]); npres = np.array([int(m["n"]) for m in meta])
    Pv[:, degen] = np.nan; Be[:, degen] = np.nan
    tested = tested0 & ~degen

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

    arch = f"{MATDIR}/prefix_degenerate"; os.makedirs(arch, exist_ok=True)
    for fn in ("disease_Pv.npy", "disease_Be.npy", "disease_meta.tsv", "disease_svorder.txt"):
        if os.path.exists(f"{MATDIR}/{fn}"):
            shutil.copy2(f"{MATDIR}/{fn}", f"{arch}/{fn}")
    if os.path.exists(HITS):
        shutil.copy2(HITS, f"{arch}/phewide_disease_full.prefix.tsv")

    np.save(f"{MATDIR}/disease_Pv.npy", Pv); np.save(f"{MATDIR}/disease_Be.npy", Be)
    with open(f"{MATDIR}/disease_meta.tsv", "w") as f:
        f.write("phecode\tlabel\tcases\tn\ttested\tdegenerate\n")
        for j in range(J):
            f.write(f"{phecodes[j]}\t{label.get(phecodes[j],'')}\t{ncase[j]}\t{npres[j]}\t{int(tested[j])}\t{int(degen[j])}\n")
    with open(HITS, "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tbeta_logOR\tchi2_p\tq_bh_tier\n")
        for i in keep:
            a = i // J; b = i % J
            f.write(f"{sv_order[a]}\t{tiers[a]}\t{phecodes[b]}\t{label.get(phecodes[b],'')}\t{ncase[b]}\t"
                    f"{npres[b]}\t{Be_flat[i]:.4g}\t{Pv_flat[i]:.3g}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    print(f"MASK_DONE: {int(degen.sum())} degenerate phenos masked; {int(tested.sum())} phecodes now valid; "
          f"{n_q05} q<0.05 (tiered); {len(keep)} hits -> {HITS} (pre-fix archived in {arch}) ({time.time()-t0:.0f}s)",
          flush=True)


if __name__ == "__main__":
    main()
