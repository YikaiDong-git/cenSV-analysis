#!/usr/bin/env python3
"""Step 62 (ext2 arm) — degenerate-null-fit VERIFICATION (the ext2 counterpart of step30i)."""
import os, sys, csv, time, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S

MATDIR = S.MATDIR
OUTA = S.OUTA
HITS = f"{OUTA}/phewide_disease_full.tsv"
FLOOD_FRAC = 0.02
FLOOD_P = 1e-6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="write masked outputs to a NEW tree (never overwrites the published ext2 files)")
    a = ap.parse_args()
    t0 = time.time()

    meta = list(csv.DictReader(open(f"{MATDIR}/disease_meta.tsv"), delimiter="\t"))
    phecodes = [m["phecode"] for m in meta]; J = len(phecodes)
    tested0 = np.array([m["tested"] == "1" for m in meta])
    Pv = np.load(f"{MATDIR}/disease_Pv.npy")
    M = Pv.shape[0]
    with np.errstate(invalid="ignore"):
        frac_flood = np.array([(np.nan_to_num(Pv[:, j], nan=1.0) < FLOOD_P).mean() for j in range(J)])
    degen = tested0 & (frac_flood > FLOOD_FRAC)

    print(f"ext2 matrices: {M} SVs x {J} phecodes, {int(tested0.sum())} tested "
          f"({time.time()-t0:.0f}s)", flush=True)
    print(f"flood detector: phecodes with >{FLOOD_FRAC:.0%} of ALL SVs at p<{FLOOD_P:g}")
    print(f"  DEGENERATE FOUND: {int(degen.sum())}", flush=True)
    if degen.sum():
        for j in np.where(degen)[0]:
            print(f"    {phecodes[j]:12s} flood_frac={frac_flood[j]:.4f} cases={meta[j]['cases']}")
    order = np.argsort(-np.where(tested0, frac_flood, -1))
    print("  worst 8 tested phecodes by flood fraction (for calibration):")
    for j in order[:8]:
        mark = "  <-- DEGENERATE" if degen[j] else ""
        print(f"    {phecodes[j]:12s} {frac_flood[j]:.6f}  cases={meta[j]['cases']}{mark}")

    if not degen.sum():
        print("\nVERIFIED: the engine's inline |Sigma r| > IDENT_TOL guard left NO phecode that floods the")
        print("result. The independent result-side detector agrees with the fit-side guard, so the ext2")
        print("headline needs no re-thresholding and nothing is rewritten.")
        print(f"DEGEN_CHECK_CLEAN ({time.time()-t0:.0f}s)", flush=True)
        return

    print("\n⚠️ The result-side detector found phecodes the inline guard did not. The ext2 hit set would")
    print("   need re-thresholding. Re-run with --apply to write masked outputs to a NEW tree.")
    if not a.apply:
        print("DEGEN_CHECK_DIRTY (report-only; pass --apply to write)", flush=True)
        return

    out = f"{OUTA}/degenmask"; os.makedirs(out, exist_ok=True)
    Be = np.load(f"{MATDIR}/disease_Be.npy")
    svrows = [ln.rstrip("\n").split("\t") for ln in open(f"{MATDIR}/disease_svorder.txt") if ln.strip()]
    sv_order = [r[0] for r in svrows]
    tiers = np.array([r[1] if len(r) > 1 else "" for r in svrows])
    label = {m["phecode"]: m["label"] for m in meta}
    ncase = np.array([int(m["cases"]) for m in meta]); npres = np.array([int(m["n"]) for m in meta])
    Pv[:, degen] = np.nan; Be[:, degen] = np.nan
    tested = tested0 & ~degen
    sv_of = np.repeat(np.arange(M), J); ph_of = np.tile(np.arange(J), M)
    Pv_flat = Pv.ravel(); Be_flat = Be.ravel()
    valid = np.isfinite(Pv_flat) & tested[ph_of]
    q = np.full(M * J, np.nan)
    for tn in sorted({t for t in tiers.tolist() if t}):
        idx = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if idx.size:
            q[idx] = S.bh(Pv_flat[idx])
    keep = np.nonzero(valid & ((np.nan_to_num(q, nan=1.0) < 0.05) | (Pv_flat < 1e-4)))[0]
    keep = keep[np.argsort(Pv_flat[keep])]
    np.save(f"{out}/disease_Pv.npy", Pv); np.save(f"{out}/disease_Be.npy", Be)
    with open(f"{out}/disease_meta.tsv", "w") as f:
        f.write("phecode\tlabel\tcases\tn\ttested\tdegenerate\n")
        for j in range(J):
            f.write(f"{phecodes[j]}\t{label.get(phecodes[j],'')}\t{ncase[j]}\t{npres[j]}\t"
                    f"{int(tested[j])}\t{int(degen[j])}\n")
    with open(f"{out}/phewide_disease_full.tsv", "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tbeta_logOR\tchi2_p\tq_bh_tier\n")
        for i in keep:
            aa = i // J; bb = i % J
            f.write(f"{sv_order[aa]}\t{tiers[aa]}\t{phecodes[bb]}\t{label.get(phecodes[bb],'')}\t{ncase[bb]}\t"
                    f"{npres[bb]}\t{Be_flat[i]:.4g}\t{Pv_flat[i]:.3g}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    print(f"MASK_APPLIED -> {out}/ ({len(keep)} hits). Published ext2 files untouched.", flush=True)


if __name__ == "__main__":
    main()
