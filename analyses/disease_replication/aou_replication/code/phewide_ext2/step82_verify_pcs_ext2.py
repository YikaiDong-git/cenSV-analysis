#!/usr/bin/env python3
"""Step 82 (ext2 arm) — adversarial validity check for the within-AMR PC rescue (step61)."""
import os, sys, csv, argparse
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S
import step51b_exact_firth_refit_ext2 as H
import step56_ancestry_disease_ext2 as P
import step61_within_amr_pc_rescue_ext2 as V

RESCUE = f"{S.OUTA}/disease_matrices_firth/per_ancestry/chr16_within_amr_pc_rescue_ext2.tsv"
OUT = f"{S.OUTA}/per_ancestry/triage/within_amr_pc_verification_ext2.tsv"


def logit_fit(Xd, y, maxit=100, tol=1e-8):
    b = np.zeros(Xd.shape[1])
    for _ in range(maxit):
        eta = np.clip(Xd @ b, -30, 30); mu = 1 / (1 + np.exp(-eta)); w = np.clip(mu * (1 - mu), 1e-9, None)
        z = eta + (y - mu) / w
        WX = Xd * w[:, None]
        b_new = np.linalg.solve(Xd.T @ WX + 1e-8 * np.eye(Xd.shape[1]), Xd.T @ (w * z))
        if np.max(np.abs(b_new - b)) < tol:
            b = b_new; break
        b = b_new
    eta = np.clip(Xd @ b, -30, 30); mu = 1 / (1 + np.exp(-eta))
    ll = float(np.sum(y * np.log(np.clip(mu, 1e-12, 1)) + (1 - y) * np.log(np.clip(1 - mu, 1e-12, 1))))
    return ll, b


def r2_lin(vec, Z):
    beta, *_ = np.linalg.lstsq(Z, vec, rcond=None)
    ss = np.sum((vec - vec.mean()) ** 2)
    return float(1 - np.sum((vec - Z @ beta) ** 2) / ss) if ss > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kpc", type=int, default=10)
    a = ap.parse_args()
    KPC = a.kpc
    rows = []

    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    gpc_idx = [j for j, c in enumerate(names) if c.startswith("gPC")]
    lab = P.ancestry_labels(X, names, N)
    amr = (lab.astype(str) == "AMR")
    amr_iids = [iids[k] for k in range(N) if amr[k]]
    gPC = X[np.ix_(amr, gpc_idx)]
    PCm, miss = V.load_within_amr_pcs(amr_iids, KPC)
    ok = np.all(np.isfinite(PCm), axis=1)
    gPC = gPC[ok]; wPC = PCm[ok, :KPC]; na = int(ok.sum())
    print(f"AMR samples with complete within-AMR PCs: {na} (missing {miss})", flush=True)

    phes, svs = [], []
    if os.path.exists(RESCUE):
        rr = list(csv.DictReader(open(RESCUE), delimiter="\t"))
        phes = sorted({r["phecode"] for r in rr})
        svs = sorted({r["sv_locus"] for r in rr})
    if not phes:
        raise SystemExit(f"ABORT: cannot derive the rescue's phecode set from {RESCUE} — "
                         f"run step61 first (hardcoding the published list would test the wrong cells).")
    print(f"rescue covers {len(svs)} SVs x {len(phes)} phecodes: {phes}\n", flush=True)

    Zg = np.column_stack([np.ones(na), gPC])
    print("A. within-AMR PC variance explained by the 16 global gPCs (is this NEW structure?):", flush=True)
    for k in range(min(5, KPC)):
        r2 = r2_lin(wPC[:, k], Zg)
        print(f"   within-AMR PC{k + 1}: R^2 on global gPCs = {r2:.3f}  "
              f"({'largely NEW structure' if r2 < 0.5 else 'overlaps global'})", flush=True)
        rows.append(dict(section="A", item=f"withinAMR_PC{k + 1}", metric="R2_on_global_gPCs", value=r2))

    YV = H.load_pheno_cols(phes, iidx)

    def zscore(Mz):
        return (Mz - Mz.mean(0)) / (Mz.std(0) + 1e-9)
    wz = zscore(wPC); gz = zscore(gPC)
    print("\nB. phenotype ~ within-AMR PCs (logistic LRT + McFadden R^2):", flush=True)
    pheno_r2 = {}
    for ph in phes:
        yf = YV.get(ph)
        if yf is None:
            continue
        ya = yf[amr][ok]; present = ya >= 0
        y = ya[present].astype(float); n1 = int(y.sum())
        if n1 < 50:
            continue
        Wp = wz[present]
        ll0, _ = logit_fit(np.ones((int(present.sum()), 1)), y)
        ll1, _ = logit_fit(np.column_stack([np.ones(int(present.sum())), Wp]), y)
        lr = 2 * (ll1 - ll0); p = float(chi2dist.sf(lr, KPC)); mcf = 1 - ll1 / ll0
        r2l = r2_lin(y, np.column_stack([np.ones(int(present.sum())), Wp]))
        pheno_r2[ph] = r2l
        print(f"   {ph:11s} (cases={n1:>5}): LRT chi2={lr:8.1f} (df={KPC}) p={p:.2e}  "
              f"McFadden R^2={mcf:.4f}  linear R^2={r2l:.4f}", flush=True)
        rows.append(dict(section="B", item=ph, metric="LRT_chi2", value=lr, extra=f"p={p:.3e}",
                         n=n1, mcfadden_r2=mcf, linear_r2=r2l))

    print("\nC. phenotype ~ gPC1-16  vs  + within-AMR PCs (do they add structure the gPCs missed?):", flush=True)
    for ph in phes[:3]:
        yf = YV.get(ph)
        if yf is None:
            continue
        ya = yf[amr][ok]; present = ya >= 0; y = ya[present].astype(float)
        if int(y.sum()) < 50:
            continue
        llg, _ = logit_fit(np.column_stack([np.ones(int(present.sum())), gz[present]]), y)
        llgw, _ = logit_fit(np.column_stack([np.ones(int(present.sum())), gz[present], wz[present]]), y)
        lr = 2 * (llgw - llg); p = float(chi2dist.sf(lr, KPC))
        print(f"   {ph:11s}: incremental LRT chi2={lr:.1f} (df={KPC}) p={p:.2e} -> within-AMR PCs "
              f"{'DO' if p < 0.05 else 'do NOT'} add phenotype-relevant structure beyond global gPCs", flush=True)
        rows.append(dict(section="C", item=ph, metric="incremental_LRT_chi2", value=lr, extra=f"p={p:.3e}"))

    print("\nD. BOUND — the most of the g~y association the within-AMR PC subspace could mediate", flush=True)
    print("   |corr(P_C g, P_C y)| <= sqrt(R2_g) * sqrt(R2_y)   (Cauchy-Schwarz; LINEAR/LPM bound)", flush=True)
    dos = H.load_sv_dosages(svs, iidx)
    Zw = np.column_stack([np.ones(na), wz])
    sv_r2 = {}
    for sv in svs:
        g = dos.get(sv)
        if g is None:
            continue
        ga = np.asarray(g, float)[amr][ok]
        sv_r2[sv] = r2_lin(ga, Zw)
    if sv_r2 and pheno_r2:
        gmax = max(sv_r2.values()); ymax = max(pheno_r2.values())
        smax = max(sv_r2, key=sv_r2.get); pmax = max(pheno_r2, key=pheno_r2.get)
        bound = float(np.sqrt(max(gmax, 0)) * np.sqrt(max(ymax, 0)))
        print(f"   max R2(SV dosage ~ within-AMR PCs)  = {gmax:.4f}   ({smax.split('_best')[0]})", flush=True)
        print(f"   max R2(phenotype ~ within-AMR PCs)  = {ymax:.4f}   ({pmax})", flush=True)
        print(f"   => the PC subspace can mediate at most |r| = {bound:.4f} of any SV-phenotype "
              f"correlation in AMR.", flush=True)
        rows.append(dict(section="D", item="max_R2_SV_on_PCs", metric="R2", value=gmax, extra=smax))
        rows.append(dict(section="D", item="max_R2_pheno_on_PCs", metric="R2", value=ymax, extra=pmax))
        rows.append(dict(section="D", item="max_mediated_correlation", metric="cauchy_schwarz_bound",
                         value=bound))
        obs = []
        for r in csv.DictReader(open(RESCUE), delimiter="\t"):
            sv, ph = r["sv_locus"], r["phecode"]
            g = dos.get(sv); yf = YV.get(ph)
            if g is None or yf is None:
                continue
            ga = np.asarray(g, float)[amr][ok]; ya = yf[amr][ok]
            pr = ya >= 0
            if pr.sum() < 100:
                continue
            c = np.corrcoef(ga[pr], ya[pr].astype(float))[0, 1]
            if np.isfinite(c):
                obs.append((abs(c), sv, ph))
        if obs:
            obs.sort(reverse=True)
            med = float(np.median([o[0] for o in obs]))
            print(f"   observed |r(SV, phenotype)| in AMR over {len(obs)} rescue cells: "
                  f"median {med:.4f}, max {obs[0][0]:.4f} ({obs[0][1].split('_best')[0]} x {obs[0][2]})",
                  flush=True)
            n_below = sum(1 for o in obs if o[0] <= bound)
            print(f"   cells whose ENTIRE observed correlation could in principle be PC-mediated: "
                  f"{n_below}/{len(obs)}", flush=True)
            rows.append(dict(section="D", item="observed_abs_r_median", metric="corr", value=med))
            rows.append(dict(section="D", item="observed_abs_r_max", metric="corr", value=obs[0][0],
                             extra=f"{obs[0][1]} x {obs[0][2]}"))
            rows.append(dict(section="D", item="cells_within_bound", metric="count", value=n_below,
                             extra=f"of {len(obs)}"))
            print("", flush=True)
            if n_below == 0:
                print("   READ: no rescue cell's observed AMR correlation is small enough to be fully", flush=True)
                print("   explained by the within-AMR PC subspace. Combined with 197/197 SURVIVING the", flush=True)
                print("   adjustment, within-AMR GENETIC stratification is refuted as the sole explanation.", flush=True)
            else:
                print(f"   READ: {n_below} cell(s) have an observed correlation at or below the bound —", flush=True)
                print("   for those, PC-mediated stratification cannot be excluded on magnitude alone.", flush=True)
    print("\n   ⚠️ SCOPE: this bounds only what the within-AMR PC subspace can mediate. Environmental", flush=True)
    print("   confounding and an AMR-enriched tag haplotype remain open (same caveat as step61).", flush=True)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    cols = ["section", "item", "metric", "value", "extra", "n", "mcfadden_r2", "linear_r2"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: (f"{r[c]:.6g}" if isinstance(r.get(c), float) and np.isfinite(r[c])
                            else r.get(c, "")) for c in cols})
    print(f"\nwrote -> {OUT}", flush=True)
    print("PC_VERIFY_EXT2_DONE", flush=True)


if __name__ == "__main__":
    main()
