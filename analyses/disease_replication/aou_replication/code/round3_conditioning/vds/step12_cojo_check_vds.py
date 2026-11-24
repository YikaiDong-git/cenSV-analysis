#!/usr/bin/env python3
"""AoU individual-level COJO analog — verify discovery's locus-specificity claim (commit 83f73be)."""
import os, csv, numpy as np, pgenlib
from step10_snp_conditioning import read_pgen_iids, load_mode, ols_dose

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
H = f"{BASE}/phewide_222_handoff"
PGEN = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen")
OUT = f"{BASE}/aou_replication/results_vds/step10_cond"
SV = "chr3_sv20_DEL_non_best6_alt_edge1000"
CHROM = "chr3"
LOINC = "29463-7"
KS = [1, 3, 5, 10, 30]


def top_predictors():
    rows = [r for r in csv.DictReader(open(f"{H}/sv_weights_long_222.tsv"), delimiter="\t")
            if r["sv_locus"] == SV]
    rows.sort(key=lambda r: -abs(float(r["beta_std"])))
    return [(int(r["pos"]), r["ref"], r["alt"]) for r in rows]


def read_snp_cols(preds, sel, Npgen):
    want = {p for p, _, _ in preds}
    offsets = [0]; pos_map = {}; vidx = 0
    with open(f"{PGEN}/acaf_threshold.{CHROM}.pvar") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t", 5)
            alts = c[4].split(",")
            offsets.append(offsets[-1] + 1 + len(alts))
            if c[1] in want or int(c[1]) in want:
                pos_map.setdefault(int(c[1]), []).append((vidx, c[3], alts))
            vidx += 1
    offs = np.array(offsets, dtype=np.uintp)
    reader = pgenlib.PgenReader(f"{PGEN}/acaf_threshold.{CHROM}.pgen".encode(),
                                raw_sample_ct=Npgen, variant_ct=vidx, allele_idx_offsets=offs)
    buf = np.empty(Npgen, dtype=np.int32)
    cols, used = [], []
    for pos, ref, alt in preds:
        hit = None
        for (vi, pref, palts) in pos_map.get(pos, []):
            if pref == ref and alt in palts:
                hit = (vi, 1 + palts.index(alt)); break
            if pref == alt and ref in palts:
                hit = (vi, 0); break
        if hit is None:
            continue
        vi, ai = hit
        reader.read(vi, buf, allele_idx=ai)
        g = buf[sel].astype(np.float64); g[g < 0] = np.nan
        m = np.nanmean(g); s = np.nanstd(g)
        if not s > 0:
            continue
        g = np.where(np.isnan(g), m, g)
        cols.append(((g - m) / s)); used.append(pos)
    reader.close()
    return cols, used


def main():
    preds = top_predictors()
    pgen_iids = read_pgen_iids(); Npgen = len(pgen_iids)
    iids, idx, C, N, dose, lirnt = load_mode("eur")
    pgen_pos = {i: k for k, i in enumerate(pgen_iids)}
    sel = np.array([pgen_pos[i] for i in iids], dtype=np.int64)
    dsv = dose[SV].astype(np.float64); y = lirnt[LOINC]
    snpcols, used = read_snp_cols(preds, sel, Npgen)
    print(f"read {len(snpcols)}/{len(preds)} predictor SNPs as conditioning covariates (top-|weight| first)\n")

    mask = np.isfinite(y) & np.isfinite(dsv)
    n = int(mask.sum()); yy = y[mask]
    def fit(extra_cols):
        cols = [np.ones(n), dsv[mask], C[mask]] + [c[mask] for c in extra_cols]
        X = np.column_stack(cols)
        keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9
        return ols_dose(X[:, keep], yy)
    b0, se0, p0 = fit([])
    print(f"M0  IRNT(weight) ~ SV_dosage + covars:   beta={b0:+.5f}  se={se0:.5f}  p={p0:.2e}  (n={n})")
    rows = [("M0", 0, b0, se0, p0, 0.0)]
    for K in KS:
        k = min(K, len(snpcols))
        bK, seK, pK = fit(snpcols[:k])
        att = (abs(b0) - abs(bK)) / abs(b0) * 100
        print(f"Mc  + top-{k:<2} predictor SNPs:           beta={bK:+.5f}  se={seK:.5f}  p={pK:.2e}  atten={att:+5.1f}%")
        rows.append((f"cond_top{k}", k, bK, seK, pK, att))
    k3 = min(3, len(snpcols))
    X3 = np.column_stack([np.ones(n), C[mask]] + [c[mask] for c in snpcols[:k3]])
    beta3 = np.linalg.lstsq(X3, yy, rcond=None)[0]
    print(f"\n(the top-{k3} predictor SNPs are the known-locus weight SNPs per discovery COJO)")
    with open(f"{OUT}/cojo_check_chr3_sv20_weight.tsv", "w") as fo:
        fo.write("model\tn_cond_snps\tbeta_sv\tse_sv\tp_sv\tpct_atten\n")
        for m, k, b, se, p, a in rows:
            fo.write(f"{m}\t{k}\t{b:.5g}\t{se:.5g}\t{p:.3g}\t{a:.1f}\n")
    print(f"\nwrote {OUT}/cojo_check_chr3_sv20_weight.tsv")


if __name__ == "__main__":
    main()
