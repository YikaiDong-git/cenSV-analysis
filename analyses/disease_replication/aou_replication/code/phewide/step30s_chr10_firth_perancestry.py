#!/usr/bin/env python3
import os, sys, csv
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S
import step30h_firth_refit_hits as H
import step30r_ancestry_quant as R

SVS = ["chr10_sv697_INS_mon_best8_fksize500k", "chr10_sv698_INS_ori_best8_fksize500k"]
PHE = "DE_664.42"; LABEL = "Psoriatic arthropathy"
ANCS = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
WELL = ["EUR", "AFR", "AMR"]
OUT = f"{S.OUTA}/chr10_firth_perancestry.tsv"


def main():
    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    Xcov = X[:, 1:]
    lab = R.ancestry_labels(X, names, N).astype(str)
    dos = H.load_sv_dosages(SVS, iidx)
    Y = H.load_pheno_cols([PHE], iidx); y = Y.get(PHE)
    if y is None:
        raise SystemExit(f"ABORT: {PHE} not a column in disease.pheno")
    print(f"phecode {PHE} ({LABEL}): pooled cases={int((y==1).sum())} controls={int((y==0).sum())}\n", flush=True)

    rows = []
    for sv in SVS:
        g = dos.get(sv)
        if g is None:
            print(f"WARN: {sv} dosage not found"); continue
        for scope in ["POOLED"] + ANCS:
            if scope == "POOLED":
                present = (y >= 0)
            else:
                present = (y >= 0) & (lab == scope)
            nc = int((y[present] == 1).sum()); nn = int((y[present] == 0).sum())
            if nc < 10 or nn < 10:
                rows.append((sv, scope, nc, nn, np.nan, np.nan, False)); continue
            bf, pf, conv = H.firth(Xcov, g, y, present)
            orf = float(np.exp(bf)) if np.isfinite(bf) else np.nan
            rows.append((sv, scope, nc, nn, orf, pf, conv))

    with open(OUT, "w") as f:
        f.write("sv_locus\tscope\tcases\tcontrols\tOR_firth\tp_firth\tconverged\n")
        for sv, sc, nc, nn, orf, pf, conv in rows:
            f.write(f"{sv}\t{sc}\t{nc}\t{nn}\t{'' if not np.isfinite(orf) else f'{orf:.4g}'}\t"
                    f"{'' if not np.isfinite(pf) else f'{pf:.3g}'}\t{int(bool(conv))}\n")
    print(f"wrote {OUT}\n", flush=True)

    print(f"  {'SV':<40}{'scope':<8}{'cases':>7}{'OR':>9}{'p_firth':>11}{'conv':>6}")
    for sv, sc, nc, nn, orf, pf, conv in rows:
        o = f"{orf:.3f}" if np.isfinite(orf) else "na"
        p = f"{pf:.2e}" if np.isfinite(pf) else "na"
        print(f"  {sv[:39]:<40}{sc:<8}{nc:>7}{o:>9}{p:>11}{int(bool(conv)):>6}", flush=True)

    print("\n=== per-SV WELL-stratum verdict (exact Firth) ===", flush=True)
    for sv in SVS:
        wr = {sc: (orf, pf, nc) for (s, sc, nc, nn, orf, pf, conv) in rows if s == sv and sc in WELL}
        sig = [a for a in WELL if a in wr and np.isfinite(wr[a][1]) and wr[a][1] < 0.05]
        signs = {np.sign(np.log(wr[a][0])) for a in sig if np.isfinite(wr[a][0]) and wr[a][0] > 0}
        pooled_or = next((orf for (s, sc, nc, nn, orf, pf, conv) in rows if s == sv and sc == "POOLED"), np.nan)
        msg = (f"  {sv[:34]:<36} sig strata={sig} sign-concordant={len(signs)==1 and len(sig)>=1} "
               f"pooled_OR={pooled_or:.2f}")
        print(msg, flush=True)
    print("CHR10_FIRTH_DONE", flush=True)


if __name__ == "__main__":
    main()
