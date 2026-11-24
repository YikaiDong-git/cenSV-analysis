#!/usr/bin/env python3
"""step30t — per-ancestry EXACT FIRTH of the q in [0.05,0.1) suggestive-band REPLICATES signals."""
import os, sys, csv
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S
import step30h_firth_refit_hits as H
import step30r_ancestry_quant as R

WELL = ["EUR", "AFR", "AMR"]
SCOPES = ["POOLED", "EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
TARGETS = sys.argv[1] if len(sys.argv) > 1 else "/tmp/q0510_firth_targets.tsv"
OUT = f"{S.OUTA}/q0510_firth_perancestry.tsv"
MIN_CASES = 20


def main():
    tgt = list(csv.DictReader(open(TARGETS), delimiter="\t"))
    svs = sorted({t["sv_locus"] for t in tgt}); phes = sorted({t["phecode"] for t in tgt})
    print(f"{len(tgt)} target cells: {len(svs)} SVs x {len(phes)} phecodes", flush=True)

    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    Xcov = X[:, 1:]
    lab = R.ancestry_labels(X, names, N).astype(str)
    dos = H.load_sv_dosages(svs, iidx)
    Y = H.load_pheno_cols(phes, iidx)
    print(f"loaded {len(dos)} dosages, {len(Y)} pheno cols\n", flush=True)

    rows = []
    for t in tgt:
        sv, ph, label = t["sv_locus"], t["phecode"], t["label"]
        g = dos.get(sv); y = Y.get(ph)
        if g is None or y is None:
            print(f"WARN missing {sv} / {ph}"); continue
        for scope in SCOPES:
            present = (y >= 0) if scope == "POOLED" else ((y >= 0) & (lab == scope))
            nc = int((y[present] == 1).sum()); nn = int((y[present] == 0).sum())
            if nc < MIN_CASES or nn < MIN_CASES:
                rows.append((sv, ph, label, scope, nc, nn, np.nan, np.nan, False)); continue
            bf, pf, conv = H.firth(Xcov, g, y, present)
            orf = float(np.exp(bf)) if np.isfinite(bf) else np.nan
            rows.append((sv, ph, label, scope, nc, nn, orf, pf, conv))

    with open(OUT, "w") as f:
        f.write("sv_locus\tphecode\tlabel\tscope\tcases\tcontrols\tOR_firth\tp_firth\tconverged\n")
        for sv, ph, label, sc, nc, nn, orf, pf, conv in rows:
            f.write(f"{sv}\t{ph}\t{label}\t{sc}\t{nc}\t{nn}\t{'' if not np.isfinite(orf) else f'{orf:.4g}'}\t"
                    f"{'' if not np.isfinite(pf) else f'{pf:.3g}'}\t{int(bool(conv))}\n")
    print(f"wrote {OUT}\n", flush=True)

    cur = None
    for sv, ph, label, sc, nc, nn, orf, pf, conv in rows:
        if (sv, ph) != cur:
            cur = (sv, ph)
            print(f"\n{sv}  {ph} ({label[:38]})", flush=True)
            print(f"  {'scope':<8}{'cases':>7}{'OR':>9}{'p_firth':>11}{'conv':>6}")
        o = f"{orf:.3f}" if np.isfinite(orf) else "na"
        p = f"{pf:.2e}" if np.isfinite(pf) else "na"
        print(f"  {sc:<8}{nc:>7}{o:>9}{p:>11}{int(bool(conv)):>6}", flush=True)

    print("\n=== per-(SV,phecode) WELL verdict (exact Firth) ===", flush=True)
    bycell = {}
    for sv, ph, label, sc, nc, nn, orf, pf, conv in rows:
        bycell.setdefault((sv, ph, label), {})[sc] = (orf, pf, nc)
    for (sv, ph, label), d in bycell.items():
        sig = [a for a in WELL if a in d and np.isfinite(d[a][1]) and d[a][1] < 0.05
               and np.isfinite(d[a][0]) and d[a][0] > 0]
        signs = {np.sign(np.log(d[a][0])) for a in sig}
        pooled = d.get("POOLED", (np.nan,))[0]
        conc = len(sig) >= 2 and len(signs) == 1
        tag = "CONFIRMED(>=2 well, same sign)" if conc else ("1-well-only" if len(sig) == 1 else "not-confirmed")
        print(f"  {sv.split('_best')[0]:<24}{ph:<10}{label[:26]:<27} sig={sig} pooledOR={pooled:.2f} -> {tag}",
              flush=True)
    print("Q0510_FIRTH_DONE", flush=True)


if __name__ == "__main__":
    main()
