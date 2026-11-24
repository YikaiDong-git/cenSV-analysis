#!/usr/bin/env python3
"""Step 67 (ext2 arm) — per-ancestry EXACT Firth for the chr10 → psoriatic arthropathy cells (← step30s)."""
import os, sys, csv
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S
import step51b_exact_firth_refit_ext2 as H
import step56_ancestry_disease_ext2 as P

PHE = "DE_664.42"; LABEL = "Psoriatic arthropathy"
PUBLISHED = ["chr10_sv697_INS_mon_best8_fksize500k", "chr10_sv698_INS_ori_best8_fksize500k"]
ANCS = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
WELL = ["EUR", "AFR", "AMR"]
FLAG = f"{S.OUTA}/phewide_disease_hits_enriched_flagged.tsv"
OUT = f"{S.OUTA}/chr10_firth_perancestry_ext2.tsv"


def targets():
    found = []
    if os.path.exists(FLAG):
        found = sorted({r["sv_locus"] for r in csv.DictReader(open(FLAG), delimiter="\t")
                        if r["sv_locus"].startswith("chr10") and r["phecode"] == PHE})
    svo = {l.split("\t")[0] for l in open(f"{S.MATDIR}/disease_svorder.txt") if l.strip()}
    out = list(found)
    for p in PUBLISHED:
        if p not in out:
            out.append(p)
    missing = [s for s in out if s not in svo]
    out = [s for s in out if s in svo]
    if missing:
        print(f"⚠️ not scored in ext2 (skipped): {missing}", flush=True)
    return out, found


def main():
    svs, found = targets()
    print(f"targets: {len(svs)} chr10 × {PHE} loci "
          f"({len(found)} from the ext2 hit set, published pair force-included)", flush=True)
    for s in svs:
        print(f"  {s}{'   [published]' if s in PUBLISHED else '   [ext2-new]'}")

    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    Xcov = X[:, 1:]
    lab = P.ancestry_labels(X, names, N).astype(str)
    dos = H.load_sv_dosages(svs, iidx)
    Y = H.load_pheno_cols([PHE], iidx); y = Y.get(PHE)
    if y is None:
        raise SystemExit(f"ABORT: {PHE} not a column in disease.pheno")
    print(f"\nphecode {PHE} ({LABEL}): pooled cases={int((y==1).sum())} controls={int((y==0).sum())}\n",
          flush=True)

    rows = []
    for sv in svs:
        g = dos.get(sv)
        if g is None:
            print(f"WARN: {sv} dosage not found"); continue
        for scope in ["POOLED"] + ANCS:
            present = (y >= 0) if scope == "POOLED" else ((y >= 0) & (lab == scope))
            nc = int((y[present] == 1).sum()); nn = int((y[present] == 0).sum())
            if nc < 10 or nn < 10:
                rows.append((sv, scope, nc, nn, np.nan, np.nan, False)); continue
            bf, pf, conv = H.firth(Xcov, g, y, present)
            orf = float(np.exp(bf)) if np.isfinite(bf) else np.nan
            rows.append((sv, scope, nc, nn, orf, pf, conv))

    with open(OUT, "w") as f:
        f.write("sv_locus\tscope\tcases\tcontrols\tOR_firth\tp_firth\tconverged\tpublished_locus\n")
        for sv, sc, nc, nn, orf, pf, conv in rows:
            f.write(f"{sv}\t{sc}\t{nc}\t{nn}\t{'' if not np.isfinite(orf) else f'{orf:.4g}'}\t"
                    f"{'' if not np.isfinite(pf) else f'{pf:.3g}'}\t{int(bool(conv))}\t"
                    f"{int(sv in PUBLISHED)}\n")
    print(f"wrote {OUT}\n", flush=True)

    print(f"  {'SV':<38}{'scope':<8}{'cases':>7}{'OR':>10}{'p_firth':>11}{'conv':>6}")
    for sv, sc, nc, nn, orf, pf, conv in rows:
        o = f"{orf:.3f}" if np.isfinite(orf) else "na"
        p = f"{pf:.2e}" if np.isfinite(pf) else "na"
        print(f"  {sv[:37]:<38}{sc:<8}{nc:>7}{o:>10}{p:>11}{int(bool(conv)):>6}", flush=True)

    print("\n=== per-SV WELL-stratum verdict (exact Firth) ===", flush=True)
    for sv in svs:
        wr = {sc: (orf, pf, nc) for (s, sc, nc, nn, orf, pf, conv) in rows if s == sv and sc in WELL}
        sig = [a for a in WELL if a in wr and np.isfinite(wr[a][1]) and wr[a][1] < 0.05]
        signs = {np.sign(np.log(wr[a][0])) for a in sig if np.isfinite(wr[a][0]) and wr[a][0] > 0}
        pooled = next((orf for (s, sc, nc, nn, orf, pf, conv) in rows if s == sv and sc == "POOLED"), np.nan)
        conc = (len(signs) == 1 and len(sig) >= 1)
        sigs = ",".join(sig) if sig else "-"
        print(f"  {sv[:34]:<36} sig_WELL={sigs:<12} n_sig={len(sig)} sign-concordant={conc}  "
              f"pooled_OR={pooled:.2f}", flush=True)
    print("EXT2_CHR10_FIRTH_DONE", flush=True)


if __name__ == "__main__":
    main()
