#!/usr/bin/env python3
"""Step 61 (ext2 arm) — within-AMR PC rescue for the ext2 chr16 -> substance-use cluster."""
import os, sys, csv, argparse
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S
import step56_ancestry_disease_ext2 as P
import step51b_exact_firth_refit_ext2 as H

EIGENVEC = f"{S.OUTA}/pca/amr_pca.eigenvec"
PA = f"{S.OUTA}/disease_matrices_firth/per_ancestry"
CVR2 = f"{S.SRC}/sv_cv_r2_all.tsv"
OUT = f"{PA}/chr16_within_amr_pc_rescue_ext2.tsv"


def load_within_amr_pcs(amr_iids, kpc):
    if not os.path.exists(EIGENVEC):
        raise SystemExit(f"ABORT: {EIGENVEC} not found — run step60_build_within_amr_pcs.py first.")
    rows = {}
    with open(EIGENVEC) as f:
        head = f.readline().split()
        pc_cols = [i for i, c in enumerate(head) if c.upper().startswith("PC")]
        if not pc_cols:
            raise SystemExit(f"ABORT: no PC columns in {EIGENVEC} header: {head[:6]}")
        iid_col = 1 if len(head) > 1 and head[0].upper().lstrip("#") in ("FID", "IID") and \
            head[1].upper() == "IID" else 0
        for ln in f:
            p = ln.split()
            rows[p[iid_col]] = [float(p[i]) for i in pc_cols[:kpc]]
    PC = np.full((len(amr_iids), kpc), np.nan)
    miss = 0
    for i, iid in enumerate(amr_iids):
        v = rows.get(iid) or rows.get(str(iid))
        if v is None:
            miss += 1
        else:
            PC[i, :len(v)] = v
    return PC, miss


def r2_on_pcs(vec, Z):
    beta, *_ = np.linalg.lstsq(Z, vec, rcond=None)
    ss = np.sum((vec - vec.mean()) ** 2)
    return float(1 - np.sum((vec - Z @ beta) ** 2) / ss) if ss > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kpc", type=int, default=10)
    ap.add_argument("--sv-prefix", default="chr16_")
    ap.add_argument("--ph-prefix", default="MB_28")
    a = ap.parse_args()

    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    keep_cols = [j for j, c in enumerate(names) if c not in P.ANC_DUMMIES]
    lab = P.ancestry_labels(X, names, N)
    amr = (lab.astype(str) == "AMR"); na = int(amr.sum())
    amr_iids = [iids[k] for k in range(N) if amr[k]]
    Xa = np.ascontiguousarray(X[np.ix_(amr, keep_cols)])
    PC, miss = load_within_amr_pcs(amr_iids, a.kpc)
    ok = np.all(np.isfinite(PC), axis=1)
    print(f"AMR n={na}; within-AMR PCs K={a.kpc}; {miss} AMR samples missing a PC row; "
          f"{int(ok.sum())}/{na} with complete PCs", flush=True)
    if ok.sum() < 0.8 * na:
        print(f"  ⚠️ only {100*ok.sum()/na:.1f}% of AMR samples have PCs — interpret with care", flush=True)

    hitf = f"{PA}/disease_hits_final_AMR.tsv"
    if not os.path.exists(hitf):
        raise SystemExit(f"ABORT: {hitf} not found — run step57 first.")
    hits = [h for h in csv.DictReader(open(hitf), delimiter="\t")
            if h["tier_label"] == "significant" and h["sv_locus"].startswith(a.sv_prefix)
            and h["phecode"].startswith(a.ph_prefix)]
    seen = {}
    for h in hits:
        k = (h["sv_locus"], h["phecode"])
        if k not in seen or float(h["p_final"]) < float(seen[k]["p_final"]):
            seen[k] = h
    hits = sorted(seen.values(), key=lambda h: float(h["p_final"]))
    loci = sorted({h["sv_locus"] for h in hits}); phes = sorted({h["phecode"] for h in hits})
    cvmap = {r["sv_locus"]: float(r["cv_r2"]) for r in csv.DictReader(open(CVR2), delimiter="\t")}
    print(f"target panel: {len(hits)} AMR-significant cells "
          f"({len(loci)} {a.sv_prefix}* SVs x {len(phes)} {a.ph_prefix}* phecodes)", flush=True)
    if not hits:
        print("nothing to rescue — the cluster is not AMR-significant under ext2 weights."); return

    DOS = H.load_sv_dosages(loci, iidx)
    YV = H.load_pheno_cols(phes, iidx)
    Zpc = np.column_stack([np.ones(int(ok.sum())), PC[ok]])

    print("\n== SV dosage ~ within-AMR PCs (R^2 = how much the SV is itself an admixture proxy) ==")
    sv_r2 = {}
    for sv in loci:
        g = DOS[sv][amr][ok].astype(np.float64)
        sv_r2[sv] = r2_on_pcs(g, Zpc)
        print(f"  {'_'.join(sv.split('_')[:2]):16s} cv_r2={cvmap.get(sv,-1):5.2f}  "
              f"R^2(dosage~PCs)={sv_r2[sv]:.4f}", flush=True)

    print("\n== CONTROL: does the PHENOTYPE track the within-AMR PCs? (if not, 'survives' is vacuous) ==")
    ph_r2 = {}
    for ph in phes:
        y = YV[ph][amr]
        pres = (y >= 0) & ok
        if pres.sum() < 50:
            continue
        Zp = np.column_stack([np.ones(int(pres.sum())), PC[pres]])
        ph_r2[ph] = r2_on_pcs(y[pres].astype(np.float64), Zp)
    for ph in sorted(ph_r2, key=lambda k: -ph_r2[k])[:8]:
        print(f"  {ph:12s} R^2(phenotype~PCs)={ph_r2[ph]:.4f}", flush=True)
    print(f"  -> max phenotype R^2 {max(ph_r2.values()):.4f} over {len(ph_r2)} phecodes"
          if ph_r2 else "  -> no phecode testable", flush=True)

    Xcov0 = Xa[:, 1:]
    Xcov1 = np.column_stack([Xa[:, 1:], PC[:, :a.kpc]])
    print(f"\n== refit: Model0 (global gPCs) vs Model1 (+{a.kpc} within-AMR PCs) ==")
    print(f"{'sv':16s} {'cv':>5} {'phecode':11s} {'cases':>6} | {'b0':>7} {'p0':>10} | "
          f"{'b1':>7} {'p1':>10} | verdict", flush=True)
    out = []
    for h in hits:
        sv = h["sv_locus"]; ph = h["phecode"]
        g_full = DOS.get(sv); y_full = YV.get(ph)
        if g_full is None or y_full is None:
            continue
        g = g_full[amr].astype(np.float64); y = y_full[amr]
        present = (y >= 0) & ok
        b0, p0, _c0 = H.firth(Xcov0, g, y, present)
        b1, p1, _c1 = H.firth(Xcov1, g, y, present)
        surv = bool(np.isfinite(p1) and p1 < 1e-4 and np.sign(b1) == np.sign(b0))
        atten = bool(np.isfinite(p0) and np.isfinite(p1) and p1 > 0 and p0 > 0
                     and (np.log10(p1) - np.log10(p0)) > 2)
        verdict = "SURVIVES" if surv else ("attenuated" if atten else "weak/ns")
        cvv = cvmap.get(sv, -1)
        print(f"{'_'.join(sv.split('_')[:2]):16s} {cvv:>5.2f} {ph:11s} {h['cases']:>6} | "
              f"{b0:>7.2f} {p0:>10.1e} | {b1:>7.2f} {p1:>10.1e} | {verdict}", flush=True)
        out.append((sv, ph, cvv, h["cases"], sv_r2.get(sv, np.nan), ph_r2.get(ph, np.nan),
                    b0, p0, b1, p1, verdict))

    with open(OUT, "w") as f:
        f.write("sv_locus\tphecode\tcv_r2\tcases\tsv_dosage_r2_on_withinAMR_PCs\tphenotype_r2_on_withinAMR_PCs\t"
                "beta_model0\tp_model0\tbeta_model1_withinPC\tp_model1_withinPC\tverdict\n")
        for r in out:
            f.write("\t".join(str(x) for x in r) + "\n")
    nsurv = sum(1 for r in out if r[10] == "SURVIVES")
    natt = sum(1 for r in out if r[10] == "attenuated")
    print(f"\nSUMMARY: {nsurv}/{len(out)} cells SURVIVE within-AMR PC adjustment (p<1e-4, concordant sign); "
          f"{natt} attenuated")
    print(f"  max SV dosage R^2 on PCs {max(sv_r2.values()):.4f}; "
          f"max phenotype R^2 on PCs {max(ph_r2.values()):.4f}" if ph_r2 else "")
    print(f"  -> {OUT}")
    print("  Reminder: SURVIVES refutes within-AMR GENETIC stratification only; environmental confounding "
          "and an AMR-enriched tag-haplotype remain open.")
    print("EXT2_AMR_PC_RESCUE_DONE", flush=True)


if __name__ == "__main__":
    main()
