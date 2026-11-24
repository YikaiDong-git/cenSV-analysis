#!/usr/bin/env python3
"""step30v — Tier-1/2 robustness triage of the confirmed cross-ancestry disease signals."""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step30g_numpy_disease_scoretest as S
import step30h_firth_refit_hits as H
import step30r_ancestry_quant as R

REL = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/"
                         "aux/relatedness/samples_relatedness_flagged_samples.tsv")
SIG = [
 ("chr10_sv697(PsA)",   "chr10_sv697_INS_mon_best8_fksize500k",    "DE_664.42"),
 ("chr19_sv1121(hep)",  "chr19_sv1121_INS_ori_best8_midpoint1000", "GI_540.1"),
 ("chr19_sv438(hep)",   "chr19_sv438_INS_ori_best8_midpoint1000",  "GI_540"),
 ("chr19_sv880(hep)",   "chr19_sv880_DEL_ori_best8_midpoint500",   "GI_540"),
 ("chr19_sv229(GIobs)", "chr19_sv229_INS_non_best8_edge500",       "GI_514"),
 ("chr12_sv823(metsyn)","chr12_sv823_INS_ori_best8_fksize100k",    "EM_203"),
 ("chr3_sv49(VTE)",     "chr3_sv49_INS_mon_best6_alt_fksize100k",  "CV_440"),
 ("chr11_sv121(HTN)",   "chr11_sv121_DEL_non_best8_fksize100k",    "CV_401"),
]


def main():
    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx); Xcov = X[:, 1:]
    lab = R.ancestry_labels(X, names, N).astype(str)
    dos = H.load_sv_dosages([s[1] for s in SIG], iidx)
    keys = [s[0] for s in SIG]
    M = np.vstack([dos[s[1]] for s in SIG]).astype(np.float64)

    print("=== 1. INDEPENDENCE (pairwise dosage r) ===")
    C = np.corrcoef(M)
    print("        " + "".join(f"{k.split('(')[0][-6:]:>8}" for k in keys))
    for i, k in enumerate(keys):
        print(f"  {k.split('(')[0][-7:]:<7}" + "".join(f"{C[i,j]:>8.2f}" for j in range(len(keys))))

    print("\n=== 2. PER-ANCESTRY DOSAGE INFORMATIVENESS (mean ~2*MAF ; sd = information) ===")
    print(f"  {'signal':<22}{'EUR m/sd':>15}{'AFR m/sd':>15}{'AMR m/sd':>15}")
    for i, k in enumerate(keys):
        g = M[i]; row = f"  {k:<22}"
        for a in ("EUR", "AFR", "AMR"):
            m = lab == a; row += f"{g[m].mean():>7.2f}/{g[m].std():<7.2f}"
        print(row)

    print("\n=== 3. RELATEDNESS (AoU v8 drop-set) + unrelated-only refit ===")
    flag = set(l.strip() for l in open(REL))
    rel = np.array([iids[i] in flag for i in range(N)])
    print(f"  related in cohort: {int(rel.sum())} ({100*rel.mean():.1f}%)")
    Y = H.load_pheno_cols([s[2] for s in SIG], iidx)
    print(f"  {'signal':<22}{'FULL OR(p)':>20}{'UNRELATED OR(p)':>22}")
    for short, full, ph in SIG:
        g = dos[full]; y = Y[ph]

        def fit(mask):
            pr = (y >= 0) & mask
            b, p, c = H.firth(Xcov, g, y, pr)
            return (np.exp(b) if np.isfinite(b) else np.nan), p
        of, pf = fit(np.ones(N, bool)); ou, pu = fit(~rel)
        print(f"  {short:<22}{of:>8.2f}({pf:>8.1e}){ou:>12.2f}({pu:>8.1e})")
    print("TRIAGE_DONE")


if __name__ == "__main__":
    main()
