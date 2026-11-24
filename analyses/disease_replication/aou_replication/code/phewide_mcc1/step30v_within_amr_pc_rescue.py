#!/usr/bin/env python3
"""step30v — within-AMR PC rescue test for the AMR chr16 -> substance-use cluster."""
import os, sys, csv
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S
import step30p_ancestry_disease as P
import step30h_firth_refit_hits as H

EIGENVEC = "/tmp/phewide/pca/amr_pca.eigenvec"
KPC = int(sys.argv[sys.argv.index("--kpc")+1]) if "--kpc" in sys.argv else 10
SP = S.OUTA
PA = f"{SP}/disease_matrices_firth/per_ancestry"


def load_within_amr_pcs(amr_iids):
    rows = {}
    with open(EIGENVEC) as f:
        hdr = f.readline().rstrip("\n").split("\t") if "\t" in open(EIGENVEC).readline() else None
    with open(EIGENVEC) as f:
        head = f.readline().split()
        pc_cols = [i for i, c in enumerate(head) if c.upper().startswith("PC")]
        iid_col = 0
        for ln in f:
            p = ln.split()
            rows[p[iid_col]] = [float(p[i]) for i in pc_cols[:KPC]]
    PC = np.full((len(amr_iids), KPC), np.nan)
    miss = 0
    for i, iid in enumerate(amr_iids):
        v = rows.get(iid) or rows.get(str(iid))
        if v is None:
            miss += 1
        else:
            PC[i, :len(v)] = v
    return PC, miss


def main():
    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    keep_cols = [j for j, c in enumerate(names) if c not in P.ANC_DUMMIES]
    lab = P.ancestry_labels(X, names, N)
    amr = (lab.astype(str) == "AMR"); na = int(amr.sum())
    amr_iids = [iids[k] for k in range(N) if amr[k]]
    Xa = np.ascontiguousarray(X[np.ix_(amr, keep_cols)])
    PC, miss = load_within_amr_pcs(amr_iids)
    print(f"AMR n={na}; within-AMR PCs loaded (K={KPC}); {miss} AMR samples missing a PC row", flush=True)
    ok = np.all(np.isfinite(PC), axis=1)
    print(f"samples with complete within-AMR PCs: {int(ok.sum())}/{na}", flush=True)

    hits = [h for h in csv.DictReader(open(f"{PA}/disease_hits_final_AMR.tsv"), delimiter="\t")
            if h["tier_label"] == "significant" and h["sv_locus"].startswith("chr16_") and h["phecode"].startswith("MB_28")]
    seen = {}
    for h in hits:
        k = (h["sv_locus"], h["phecode"])
        if k not in seen or float(h["p_final"]) < float(seen[k]["p_final"]): seen[k] = h
    hits = sorted(seen.values(), key=lambda h: float(h["p_final"]))
    loci = sorted({h["sv_locus"] for h in hits}); phes = sorted({h["phecode"] for h in hits})
    cvmap = {r["sv_locus"]: float(r["cv_r2"]) for r in csv.DictReader(open(f"{SP}/sv_cv_r2_all.tsv"), delimiter="\t")}
    print(f"target panel: {len(hits)} cells ({len(loci)} chr16 SVs x {len(phes)} MB_28x phecodes)", flush=True)

    DOS = H.load_sv_dosages(loci, iidx)
    YV = H.load_pheno_cols(phes, iidx)

    print("\n== SV dosage vs within-AMR PCs (R^2 = how much the SV is an admixture proxy) ==", flush=True)
    Zpc = np.column_stack([np.ones(int(ok.sum())), PC[ok]])
    for sv in loci:
        g = DOS[sv][amr][ok].astype(np.float64)
        beta, *_ = np.linalg.lstsq(Zpc, g, rcond=None); r2 = 1 - np.sum((g - Zpc@beta)**2)/np.sum((g-g.mean())**2)
        print(f"  {'_'.join(sv.split('_')[:2]):14s} cv_r2={cvmap.get(sv,-1):.2f}  R^2(dosage~PCs)={r2:.3f}", flush=True)

    Xcov0 = Xa[:, 1:]
    Xcov1 = np.column_stack([Xa[:, 1:], PC[:, :KPC]])
    print(f"\n== refit: Model0 (global gPCs) vs Model1 (+{KPC} within-AMR PCs) ==", flush=True)
    print(f"{'sv':14s} {'cv':>4} {'phecode':10s} {'cases':>5} | {'b0':>6} {'p0':>9} | {'b1':>6} {'p1':>9} | verdict", flush=True)
    out = []
    for h in hits:
        sv = h["sv_locus"]; ph = h["phecode"]
        g_full = DOS.get(sv); y_full = YV.get(ph)
        if g_full is None or y_full is None: continue
        g = g_full[amr].astype(np.float64); y = y_full[amr]
        present = (y >= 0) & ok
        b0, p0, c0 = H.firth(Xcov0, g, y, present)
        b1, p1, c1 = H.firth(Xcov1, g, y, present)
        surv = (np.isfinite(p1) and p1 < 1e-4 and np.sign(b1) == np.sign(b0))
        atten = (np.isfinite(p0) and np.isfinite(p1) and p1 > 0 and p0 > 0 and (np.log10(p1) - np.log10(p0)) > 2)
        verdict = "SURVIVES" if surv else ("attenuated" if atten else "weak/ns")
        cvv = cvmap.get(sv, -1)
        print(f"{'_'.join(sv.split('_')[:2]):14s} {cvv:>4.2f} {ph:10s} {h['cases']:>5} | "
              f"{b0:>6.2f} {p0:>9.1e} | {b1:>6.2f} {p1:>9.1e} | {verdict}", flush=True)
        out.append((sv, ph, cvv, h["cases"], b0, p0, b1, p1, verdict))

    with open(f"{PA}/chr16_within_amr_pc_rescue.tsv", "w") as f:
        f.write("sv_locus\tphecode\tcv_r2\tcases\tbeta_model0\tp_model0\tbeta_model1_withinPC\tp_model1_withinPC\tverdict\n")
        for r in out:
            f.write("\t".join(str(x) for x in r) + "\n")
    nsurv = sum(1 for r in out if r[8] == "SURVIVES"); natt = sum(1 for r in out if r[8] == "attenuated")
    print(f"\nSUMMARY: {nsurv}/{len(out)} cells SURVIVE within-AMR PC adjustment (p<1e-4, concordant); "
          f"{natt} attenuated -> {PA}/chr16_within_amr_pc_rescue.tsv", flush=True)


if __name__ == "__main__":
    main()
