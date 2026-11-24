#!/usr/bin/env python3
import os, csv, numpy as np
from step10_snp_conditioning import read_pgen_iids, load_mode, ols_dose
from step14_leads_conditional import read_leads_chrom, TRAIT_LOINC, TRAIT_NAME, LEADS, OUT


def main():
    groups = {}
    for l in csv.DictReader(open(LEADS), delimiter="\t"):
        groups.setdefault((l["sv_locus"], l["trait_code"]), []).append(l)
    pgen_iids = read_pgen_iids(); Npgen = len(pgen_iids)
    iids, idx, C, N, dose, lirnt = load_mode("eur")
    pgen_pos = {i: k for k, i in enumerate(pgen_iids)}
    sel = np.array([pgen_pos[i] for i in iids], dtype=np.int64)
    hdr = ("sv\ttrait\tpos\tr2_1kgp\tis_predictor\tr2_aou_dosage\tlead_beta_trait\tlead_p_trait\tinterpretation\n")
    out = []
    for (sv, tc), ls in groups.items():
        ls = sorted(ls, key=lambda l: -float(l["r2_with_predictor"]))
        loinc = TRAIT_LOINC[tc]; y = lirnt[loinc]; dsv = dose[sv].astype(np.float64)
        cols = read_leads_chrom(ls[0]["chrom"], ls, sel, Npgen)
        my = np.isfinite(y) & np.isfinite(dsv)
        for l in ls[:3]:
            pos = int(l["pos_hg38"])
            if pos not in cols:
                continue
            g = cols[pos]
            r2_aou = float(np.corrcoef(dsv[my], g[my])[0, 1] ** 2)
            X = np.column_stack([np.ones(my.sum()), g[my], C[my]])
            keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9
            lb, lse, lp = ols_dose(X[:, keep], y[my])
            r2_1kgp = float(l["r2_with_predictor"])
            weak = (r2_aou < 0.05) or (lp > 0.05)
            interp = "WEAK-INSTRUMENT" if weak else "informative"
            out.append([sv[:20], TRAIT_NAME[tc], pos, f"{r2_1kgp:.3f}", l["is_predictor_weight_snp"],
                        f"{r2_aou:.3f}", f"{lb:.4g}", f"{lp:.3g}", interp])
            print(f"  {TRAIT_NAME[tc]:7} {pos} r2_1kgp={r2_1kgp:.2f} pred={l['is_predictor_weight_snp']:5} "
                  f"| r2_AoU(dosage)={r2_aou:.3f} lead~trait b={lb:+.4f} p={lp:.1e} -> {interp}", flush=True)
    with open(f"{OUT}/lead_diagnostic_eur.tsv", "w") as fo:
        fo.write(hdr)
        for r in out:
            fo.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote {OUT}/lead_diagnostic_eur.tsv")


if __name__ == "__main__":
    main()
