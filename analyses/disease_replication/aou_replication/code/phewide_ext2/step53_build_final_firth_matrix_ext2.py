#!/usr/bin/env python3
"""Step 30h3 (DISEASE) — build the FINAL Firth-recalibrated matrix = EXACT Firth primary + approx-Firth tail."""

import os, csv, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S

D = f"{S.OUTA}"
MAT = f"{D}/disease_matrices"
OUT = f"{D}/disease_matrices_firth"; os.makedirs(OUT, exist_ok=True)


def main():
    Pv = np.load(f"{MAT}/disease_Pv.npy").astype(np.float64)
    Be = np.load(f"{MAT}/disease_Be.npy").astype(np.float64)
    svloc, tiers = [], []
    for ln in open(f"{MAT}/disease_svorder.txt"):
        if ln.strip():
            a = ln.rstrip("\n").split("\t"); svloc.append(a[0]); tiers.append(a[1] if len(a) > 1 else "")
    tiers = np.array(tiers); svidx = {s: i for i, s in enumerate(svloc)}
    meta = list(csv.DictReader(open(f"{MAT}/disease_meta.tsv"), delimiter="\t"))
    phidx = {m["phecode"]: j for j, m in enumerate(meta)}
    ph_label = {m["phecode"]: m.get("label", "") for m in meta}
    ph_cases = {m["phecode"]: int(m["cases"]) for m in meta}
    ph_n = {m["phecode"]: int(m["n"]) for m in meta}
    tested = np.array([m["tested"] == "1" for m in meta])
    M, J = Pv.shape

    Pf = Pv.copy(); Bf = Be.copy()
    source = {}
    sepmap = {}

    approx = {}
    for r in csv.DictReader(open(f"{OUT}/disease_firth_cells.tsv"), delimiter="\t"):
        try:
            approx[(r["sv_locus"], r["phecode"])] = (
                float(r["p_firth"]), float(r["beta_firth"]), int(r["separation"]),
                float(r["fisher_p"]) if r["fisher_p"] not in ("", "nan") else np.nan)
        except (ValueError, KeyError):
            continue

    n_ex = n_ap = 0
    for r in csv.DictReader(open(f"{D}/phewide_disease_full_firthrefit.tsv"), delimiter="\t"):
        k = (r["sv_locus"], r["phecode"]); i = svidx.get(k[0]); j = phidx.get(k[1])
        if i is None or j is None:
            continue
        conv = r["firth_converged"] == "1"
        if conv and r["p_firth"] not in ("", "nan"):
            Pf[i, j] = float(r["p_firth"])
            if r["beta_firth"] not in ("", "nan"):
                Bf[i, j] = float(r["beta_firth"])
            source[(i, j)] = "exact"; n_ex += 1
        elif k in approx:
            pa, ba, sep, fp = approx[k]
            Pf[i, j] = pa; Bf[i, j] = ba; source[(i, j)] = "approx"
            sepmap[(i, j)] = (sep, fp); n_ap += 1
    for (svl, ph), (pa, ba, sep, fp) in approx.items():
        i = svidx.get(svl); j = phidx.get(ph)
        if i is None or j is None or (i, j) in source:
            continue
        Pf[i, j] = pa; Bf[i, j] = ba; source[(i, j)] = "approx"
        if sep:
            sepmap[(i, j)] = (sep, fp)
    print(f"overlaid {n_ex} exact-Firth + {n_ap} approx-fallback (+ {len(source)-n_ex-n_ap} approx-only band) cells",
          flush=True)

    sv_of = np.repeat(np.arange(M), J); ph_of = np.tile(np.arange(J), M)
    flat = Pf.ravel(); Bflat = Bf.ravel()
    valid = np.isfinite(flat) & tested[ph_of]
    q = np.full(M * J, np.nan)
    for tn in sorted({t for t in tiers.tolist() if t}):
        ix = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if ix.size:
            q[ix] = S.bh(flat[ix])

    np.save(f"{OUT}/disease_Pv_final.npy", Pf.astype(np.float32))
    np.save(f"{OUT}/disease_Be_final.npy", Bf.astype(np.float32))

    keep = np.nonzero(valid & ((np.nan_to_num(q, nan=1.0) < 0.10) | (flat < 1e-4)))[0]
    keep = keep[np.argsort(flat[keep])]
    n_sig = n_ind = n_sep_hits = 0
    with open(f"{OUT}/disease_hits_final.tsv", "w") as f:
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
                n_sep_hits += 1
            elif sig:
                n_sig += 1
            elif base == "indicative":
                n_ind += 1
            ph = meta[j]["phecode"]
            f.write(f"{svloc[i]}\t{tiers[i]}\t{ph}\t{ph_label.get(ph,'')}\t{ph_cases.get(ph,0)}\t{ph_n.get(ph,0)}\t"
                    f"{Bflat[idx]:.4g}\t{flat[idx]:.3g}\t{src}\t{'' if not np.isfinite(qi) else f'{qi:.4g}'}\t"
                    f"{sep}\t{'' if not np.isfinite(fp) else f'{fp:.3g}'}\t{lab}\n")
    print(f"FINAL: {n_sig} significant (q<0.05), {n_ind} indicative (q0.05-0.1), {n_sep_hits} separation-flagged "
          f"in the hit band -> {OUT}/disease_hits_final.tsv ; matrices -> disease_{{Pv,Be}}_final.npy", flush=True)


if __name__ == "__main__":
    main()
