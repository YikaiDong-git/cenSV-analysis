#!/usr/bin/env python3
"""Step 81 (ext2 arm) — Tier-1/2 robustness triage of the CONFIRMED cross-ancestry disease signals."""
import os, sys, csv, argparse
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S
import step51b_exact_firth_refit_ext2 as H
import step56_ancestry_disease_ext2 as P

REL = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/"
                         "aux/relatedness/samples_relatedness_flagged_samples.tsv")
PANC = f"{S.OUTA}/per_ancestry"
TRIAGE = f"{PANC}/triage"
TARGETS = f"{TRIAGE}/q0510_confirmed_ext2.tsv"
ENRICH = f"{S.OUTA}/phewide_disease_hits_enriched_flagged.tsv"
OUT = f"{TRIAGE}/robustness_triage_ext2.tsv"
WELL = ("EUR", "AFR", "AMR")


def short_name(sv, ph):
    a = sv.split("_")
    return f"{a[0]}_{a[1]}({ph})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default=TARGETS)
    a = ap.parse_args()

    if not os.path.exists(a.targets):
        raise SystemExit(f"ABORT: no target list at {a.targets} — run step80 first "
                         f"(it writes the CONFIRMED set this step triages).")
    tgt = list(csv.DictReader(open(a.targets), delimiter="\t"))
    if not tgt:
        print(f"No CONFIRMED cross-ancestry signals in {a.targets} — nothing to triage.\n"
              f"(That is itself a result: report it, do not silently skip the section.)", flush=True)
        with open(OUT, "w") as f:
            f.write("# no CONFIRMED cross-ancestry band signals under ext2 weights — triage not applicable\n")
        print("TRIAGE_EXT2_DONE (empty)", flush=True)
        return

    SIG = [(short_name(t["sv_locus"], t["phecode"]), t["sv_locus"], t["phecode"],
            t.get("label", "")) for t in tgt]
    TIER = {(t["sv_locus"], t["phecode"]): t.get("tier", "") for t in tgt}
    from collections import Counter
    print(f"triaging {len(SIG)} CONFIRMED cells over {len({s[1] for s in SIG})} distinct SVs "
          f"(targets: {a.targets})", flush=True)
    print(f"  tier composition: {dict(Counter(TIER.values()))}", flush=True)
    nx = sum(1 for s in SIG if s[1].startswith("chrX"))
    if nx:
        print(f"  ⚠️ {nx} cell(s) on chrX — step28's 2-bit LUT has no ploidy handling, so these are "
              f"hypothesis-only regardless of how they triage.", flush=True)
    print("", flush=True)

    cvr2 = {}
    if os.path.exists(ENRICH):
        for r in csv.DictReader(open(ENRICH), delimiter="\t"):
            try:
                cvr2[r["sv_locus"]] = float(r["cv_r2"])
            except (ValueError, KeyError):
                pass

    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx); Xcov = X[:, 1:]
    lab = P.ancestry_labels(X, names, N).astype(str)
    uniq_sv = sorted({s[1] for s in SIG})
    dos = H.load_sv_dosages(uniq_sv, iidx)
    missing = [s for s in uniq_sv if s not in dos]
    if missing:
        raise SystemExit(f"ABORT: dosages missing for {missing} — cannot triage a signal whose "
                         f"predictor cannot be loaded.")
    Mx = np.vstack([dos[s] for s in uniq_sv]).astype(np.float64)

    print("=== 1. INDEPENDENCE (pairwise dosage r across the confirmed SVs) ===", flush=True)
    C = np.corrcoef(Mx) if len(uniq_sv) > 1 else np.ones((1, 1))
    lbl = [s.split("_")[0] + "_" + s.split("_")[1] for s in uniq_sv]
    print("        " + "".join(f"{k[-9:]:>11}" for k in lbl), flush=True)
    for i, k in enumerate(lbl):
        print(f"  {k[-9:]:<7}" + "".join(f"{C[i, j]:>11.2f}" for j in range(len(uniq_sv))), flush=True)
    seen, blocks = set(), []
    for i in range(len(uniq_sv)):
        if i in seen:
            continue
        blk = [i]; seen.add(i)
        for j in range(len(uniq_sv)):
            if j not in seen and abs(C[i, j]) >= 0.5:
                blk.append(j); seen.add(j)
        blocks.append(blk)
    print(f"\n  -> {len(SIG)} confirmed CELLS  /  {len(uniq_sv)} distinct SVs  /  "
          f"{len(blocks)} independent loci at |r|>=0.5", flush=True)
    for b in blocks:
        if len(b) > 1:
            print(f"     correlated block: {[lbl[i] for i in b]}", flush=True)

    print("\n=== 2. PER-ANCESTRY DOSAGE INFORMATIVENESS (mean ~ 2*MAF ; sd = information) ===", flush=True)
    print("  a predictor with real spread in one ancestry and ~constant in another produces an", flush=True)
    print("  'ancestry-specific effect' that is differential imputation, not biology.", flush=True)
    print(f"  {'SV':<24}{'cv_r2':>7}" + "".join(f"{a2 + ' m/sd':>17}" for a2 in WELL), flush=True)
    info = {}
    for i, sv in enumerate(uniq_sv):
        g = Mx[i]; row = f"  {lbl[i]:<24}{cvr2.get(sv, float('nan')):>7.2f}"
        info[sv] = {}
        for anc in WELL:
            m = lab == anc
            mu, sd = (float(g[m].mean()), float(g[m].std())) if m.any() else (np.nan, np.nan)
            info[sv][anc] = (mu, sd)
            row += f"{mu:>8.2f}/{sd:<8.2f}"
        print(row, flush=True)
    print("", flush=True)
    for sv in uniq_sv:
        sds = [info[sv][a2][1] for a2 in WELL if np.isfinite(info[sv][a2][1])]
        if sds and max(sds) > 0:
            ratio = min(sds) / max(sds)
            if ratio < 0.5:
                print(f"  ⚠️ {sv.split('_best')[0]}: sd varies {ratio:.2f}x across EUR/AFR/AMR "
                      f"— differential-imputation risk; an ancestry-specific effect here is not "
                      f"interpretable as biology.", flush=True)

    print("\n=== 3. RELATEDNESS (AoU v8 drop-set) + unrelated-only exact-Firth refit ===", flush=True)
    with open(REL) as f:
        first = f.readline().strip()
        flag = {ln.strip() for ln in f if ln.strip()}
        if first and first != "sample_id":
            flag.add(first)
    rel = np.array([iids[i] in flag for i in range(N)])
    print(f"  flagged samples in file: {len(flag)};  present in cohort: {int(rel.sum())} "
          f"({100 * rel.mean():.1f}% of {N})", flush=True)

    Y = H.load_pheno_cols(sorted({s[2] for s in SIG}), iidx)
    rows = []
    print(f"\n  {'signal':<26}{'FULL OR(p)':>22}{'UNRELATED OR(p)':>24}{'dOR%':>8}", flush=True)
    for short, full, ph, label in SIG:
        g = dos[full]; y = Y.get(ph)
        if y is None:
            print(f"  WARN missing phenotype {ph}", flush=True); continue

        def fit(mask):
            pr = (y >= 0) & mask
            b, p, c = H.firth(Xcov, g, y, pr)
            return (float(np.exp(b)) if np.isfinite(b) else np.nan), p, bool(c), int((y[pr] == 1).sum())

        of, pf, cf, nf = fit(np.ones(N, bool))
        ou, pu, cu, nu = fit(~rel)
        d = 100.0 * (ou - of) / of if (np.isfinite(of) and of != 0 and np.isfinite(ou)) else np.nan
        print(f"  {short:<26}{of:>10.2f}({pf:>9.1e}){ou:>12.2f}({pu:>9.1e}){d:>8.1f}", flush=True)
        rows.append(dict(sv_locus=full, tier=TIER.get((full, ph), ""), chrX=int(full.startswith("chrX")),
                         phecode=ph, label=label, short=short,
                         cv_r2=cvr2.get(full, ""),
                         EUR_mean=info[full]["EUR"][0], EUR_sd=info[full]["EUR"][1],
                         AFR_mean=info[full]["AFR"][0], AFR_sd=info[full]["AFR"][1],
                         AMR_mean=info[full]["AMR"][0], AMR_sd=info[full]["AMR"][1],
                         full_cases=nf, full_OR=of, full_p=pf, full_conv=int(cf),
                         unrel_cases=nu, unrel_OR=ou, unrel_p=pu, unrel_conv=int(cu),
                         dOR_pct=d,
                         survives_unrelated=int(np.isfinite(pu) and pu < 0.05 and cu
                                                and np.isfinite(of) and np.isfinite(ou)
                                                and np.sign(np.log(ou)) == np.sign(np.log(of)))))
    surv = sum(r["survives_unrelated"] for r in rows)
    print(f"\n  -> {surv}/{len(rows)} confirmed cells survive the unrelated-only refit "
          f"(nominal p<0.05, converged, same direction)", flush=True)

    with open(OUT, "w", newline="") as f:
        cols = list(rows[0].keys())
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) and np.isfinite(v)
                            else ("" if isinstance(v, float) else v)) for k, v in r.items()})
    print(f"\nwrote -> {OUT}", flush=True)
    print("TRIAGE_EXT2_DONE", flush=True)


if __name__ == "__main__":
    main()
