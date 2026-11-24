#!/usr/bin/env python3
"""Step 86 (ext2 arm) — is the chrX signal an artifact of unmodelled ploidy?"""
import os, sys, csv, argparse
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S
import step51b_exact_firth_refit_ext2 as H

RV = S.RV
MCC2 = f"{RV}/step30_phewide_ext2"
MCC1 = f"{RV}/step30_phewide_mcc1_ext2"
OUT = f"{MCC2}/per_ancestry/triage/chrx_ploidy_check_ext2.tsv"
MIN_CASES = 20


def chrx_sig(root, tag):
    p = f"{root}/disease_matrices_firth/disease_hits_final.tsv"
    if not os.path.exists(p):
        print(f"  ({tag}: no hit table at {p} — skipped)", flush=True)
        return []
    out = []
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r.get("tier_label") == "significant" and r["sv_locus"].startswith("chrX"):
            out.append(dict(sv=r["sv_locus"], ph=r["phecode"], label=r.get("label", ""),
                            arm=tag, q=r.get("q_bh_tier", ""), p=r.get("p_final", "")))
    print(f"  {tag}: {len(out)} significant chrX cells "
          f"({len({c['sv'] for c in out})} SVs x {len({c['ph'] for c in out})} phecodes)", flush=True)
    return out


def r2_on(vec, Z):
    beta, *_ = np.linalg.lstsq(Z, vec, rcond=None)
    ss = np.sum((vec - vec.mean()) ** 2)
    return float(1 - np.sum((vec - Z @ beta) ** 2) / ss) if ss > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="both", choices=["both", "mcc2", "mcc1"])
    a = ap.parse_args()

    print("=" * 96, flush=True)
    print("STEP 86 — chrX ploidy diagnostic (targets derived from the significant sets)", flush=True)
    print("=" * 96, flush=True)
    cells = []
    if a.arm in ("both", "mcc2"):
        cells += chrx_sig(MCC2, "MCC2")
    if a.arm in ("both", "mcc1"):
        cells += chrx_sig(MCC1, "MCC1")
    if not cells:
        print("no significant chrX cells — nothing to test (that is itself the answer).", flush=True)
        with open(OUT, "w") as f:
            f.write("# no significant chrX cells in the selected arm(s)\n")
        return

    byc = {}
    for c in cells:
        k = (c["sv"], c["ph"])
        byc.setdefault(k, dict(c, arms=set()))["arms"].add(c["arm"])
    cells = list(byc.values())
    print(f"\n{len(cells)} distinct chrX cells to test\n", flush=True)

    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    if "sex" not in names:
        raise SystemExit(f"ABORT: no 'sex' covariate among {names} — cannot stratify.")
    sx = X[:, names.index("sex")].astype(float)
    lv = np.unique(sx[np.isfinite(sx)])
    if len(lv) != 2:
        raise SystemExit(f"ABORT: 'sex' is not binary (levels={lv[:6]}) — refusing to guess the coding.")
    g0, g1 = (sx == lv[0]), (sx == lv[1])
    print(f"sex covariate levels {lv.tolist()}: "
          f"n(level={lv[0]})={int(g0.sum())}, n(level={lv[1]})={int(g1.sum())}", flush=True)
    keep = [j for j, c in enumerate(names) if c not in ("intercept", "sex")]
    Xs = X[:, keep]
    Xall = X[:, [j for j, c in enumerate(names) if c != "intercept"]]

    svs = sorted({c["sv"] for c in cells})
    dos = H.load_sv_dosages(svs, iidx)

    PHENO_INP = {"MCC2": f"{S.CARRY}/inputs", "MCC1": f"{RV}/inputs_mcc1_ext2"}
    Yarm, _saved = {}, S.INP
    for arm in sorted({x for c in cells for x in c["arms"]}):
        need = sorted({c["ph"] for c in cells if arm in c["arms"]})
        S.INP = PHENO_INP[arm]
        if not os.path.exists(f"{S.INP}/disease.pheno"):
            raise SystemExit(f"ABORT: {arm} phenotypes missing at {S.INP}/disease.pheno")
        Yarm[arm] = H.load_pheno_cols(need, iidx)
        print(f"  loaded {len(Yarm[arm])} phenotype cols for {arm} from {S.INP}", flush=True)
    S.INP = _saved
    print(f"loaded {len(dos)} chrX dosages\n", flush=True)

    Zsex = np.column_stack([np.ones(N), sx])
    rows = []
    print(f"{'SV':<34}{'phecode':<11}{'R2(dos~sex)':>12}{'sex0 OR(p)':>20}{'sex1 OR(p)':>20}  verdict", flush=True)
    for c in sorted(cells, key=lambda c: (c["sv"], c["ph"])):
        sv, ph = c["sv"], c["ph"]
        arm = "MCC2" if "MCC2" in c["arms"] else sorted(c["arms"])[0]
        g = dos.get(sv); y = Yarm.get(arm, {}).get(ph)
        if g is None or y is None:
            print(f"  WARN missing {sv} / {ph}", flush=True); continue
        gg = np.asarray(g, float)
        r2sex = r2_on(gg, Zsex)

        def fit(mask):
            pr = (y >= 0) & mask
            nc = int((y[pr] == 1).sum()); nn = int((y[pr] == 0).sum())
            if nc < MIN_CASES or nn < MIN_CASES:
                return (np.nan, np.nan, False, nc)
            Xu = Xall if np.all(mask) else Xs
            b, p, conv = H.firth(Xu, gg, y, pr)
            return (float(np.exp(b)) if np.isfinite(b) else np.nan, p, bool(conv), nc)

        o0, p0, c0, n0 = fit(g0)
        o1, p1, c1, n1 = fit(g1)
        oP, pP, cP, nP = fit(np.ones(N, bool))
        s0 = np.isfinite(p0) and p0 < 0.05 and c0
        s1 = np.isfinite(p1) and p1 < 0.05 and c1
        same = (np.isfinite(o0) and np.isfinite(o1) and np.sign(np.log(o0)) == np.sign(np.log(o1)))
        if not (c0 and c1):
            v = ("NON-CONVERGED (no verdict; Firth hit max iterations in "
                 + ("both sexes)" if not (c0 or c1) else ("sex=%g)" % (lv[0] if not c0 else lv[1]))))
        elif s0 and s1 and same:
            v = "CONCORDANT (ploidy not sole driver)"
        elif s0 and s1:
            v = "BOTH-SEXES-OPPOSITE-SIGN"
        elif s0 or s1:
            v = "ONE-SEX-ONLY (ploidy/sex artifact consistent)"
        else:
            v = "NEITHER (pooled signal does not survive stratification)"
        f0 = f"{o0:.2f}({p0:.1e})" if np.isfinite(o0) else "na"
        f1 = f"{o1:.2f}({p1:.1e})" if np.isfinite(o1) else "na"
        print(f"{sv.split('_best')[0][:34]:<34}{ph:<11}{r2sex:>12.4f}{f0:>20}{f1:>20}  {v}", flush=True)
        rows.append(dict(sv_locus=sv, phecode=ph, label=c["label"], arms=",".join(sorted(c["arms"])),
                         pheno_arm=arm,
                         r2_dosage_on_sex=r2sex,
                         sex0_level=lv[0], sex0_cases=n0, sex0_OR=o0, sex0_p=p0, sex0_conv=int(c0),
                         sex1_level=lv[1], sex1_cases=n1, sex1_OR=o1, sex1_p=p1, sex1_conv=int(c1),
                         pooled_cases=nP, pooled_OR=oP, pooled_p=pP, verdict=v))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    cols = list(rows[0].keys())
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) and np.isfinite(v)
                            else ("" if isinstance(v, float) else v)) for k, v in r.items()})

    from collections import Counter
    tal = Counter(r["verdict"] for r in rows)
    print(f"\nverdicts: {dict(tal)}", flush=True)
    hi = [r for r in rows if np.isfinite(r["r2_dosage_on_sex"]) and r["r2_dosage_on_sex"] > 0.10]
    print(f"cells whose dosage tracks sex at R^2>0.10: {len(hi)}/{len(rows)}"
          + (f"  -> {sorted({r['sv_locus'].split('_best')[0] for r in hi})}" if hi else ""), flush=True)
    print(f"\nwrote -> {OUT}", flush=True)
    print("⚠️ This is a diagnostic, not a fix: correct chrX dosage needs a ploidy-aware LUT in step50.",
          flush=True)
    print("CHRX_PLOIDY_CHECK_DONE", flush=True)


if __name__ == "__main__":
    main()
