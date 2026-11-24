#!/usr/bin/env python3
import os, csv, numpy as np, pgenlib
from step10_snp_conditioning import read_pgen_iids, load_mode, ols_dose

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
LEADS = f"{BASE}/aou_handoff/leads_for_aou_conditional.tsv"
PGEN = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen")
OUT = f"{BASE}/aou_replication/results/step10_cond"
TRAIT_LOINC = {"21002_irnt": "29463-7", "50_irnt": "8302-2", "BMI_IRN": "39156-5", "30110_irnt": "776-5"}
TRAIT_NAME = {"21002_irnt": "weight", "50_irnt": "height", "BMI_IRN": "BMI", "30110_irnt": "PDW"}


def read_leads_chrom(chrom, leads, sel, Npgen):
    want = {int(l["pos_hg38"]) for l in leads}
    offs = [0]; pmap = {}; vidx = 0
    with open(f"{PGEN}/acaf_threshold.{chrom}.pvar") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t", 5)
            alts = c[4].split(",")
            offs.append(offs[-1] + 1 + len(alts))
            if int(c[1]) in want:
                pmap.setdefault(int(c[1]), []).append((vidx, c[3], alts))
            vidx += 1
    reader = pgenlib.PgenReader(f"{PGEN}/acaf_threshold.{chrom}.pgen".encode(), raw_sample_ct=Npgen,
                                variant_ct=vidx, allele_idx_offsets=np.array(offs, dtype=np.uintp))
    buf = np.empty(Npgen, dtype=np.int32); out = {}
    for l in leads:
        pos, ref, alt = int(l["pos_hg38"]), l["ref"], l["alt_effect"]
        hit = None
        for (vi, pref, palts) in pmap.get(pos, []):
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
        if s > 0:
            out[pos] = np.where(np.isnan(g), m, (g - m) / s)
    reader.close()
    return out


def main():
    rows_leads = list(csv.DictReader(open(LEADS), delimiter="\t"))
    groups = {}
    for l in rows_leads:
        groups.setdefault((l["sv_locus"], l["trait_code"]), []).append(l)
    pgen_iids = read_pgen_iids(); Npgen = len(pgen_iids)
    HDR = ("sv\ttrait\tmode\tn\tn_leads_indep\tn_leads_all\tbeta_m0\tp_m0\tbeta_clean\tp_clean\tatten_clean\t"
           "beta_full\tp_full\tatten_full\ttop_lead_r2\ttop_lead_is_predictor\n")
    for mode in ["eur", "all"]:
        print(f"=== mode={mode} ===", flush=True)
        iids, idx, C, N, dose, lirnt = load_mode(mode)
        pgen_pos = {i: k for k, i in enumerate(pgen_iids)}
        sel = np.array([pgen_pos[i] for i in iids], dtype=np.int64)
        out_rows = []
        for (sv, tc), ls in groups.items():
            loinc = TRAIT_LOINC[tc]; y = lirnt[loinc]; dsv = dose[sv].astype(np.float64)
            cols = read_leads_chrom(ls[0]["chrom"], ls, sel, Npgen)
            indep = [cols[int(l["pos_hg38"])] for l in ls
                     if l["is_predictor_weight_snp"] != "True" and int(l["pos_hg38"]) in cols]
            allc = [cols[int(l["pos_hg38"])] for l in ls if int(l["pos_hg38"]) in cols]
            mask = np.isfinite(y) & np.isfinite(dsv); n = int(mask.sum()); yy = y[mask]

            def fit(extra):
                X = np.column_stack([np.ones(n), dsv[mask], C[mask]] + [c[mask] for c in extra])
                keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9
                return ols_dose(X[:, keep], yy)
            b0, s0, p0 = fit([]); bc, sc, pc = fit(indep); bf, sf, pf = fit(allc)
            ac = (abs(b0) - abs(bc)) / abs(b0) * 100 if b0 else float("nan")
            af = (abs(b0) - abs(bf)) / abs(b0) * 100 if b0 else float("nan")
            top = max(ls, key=lambda l: float(l["r2_with_predictor"]))
            out_rows.append([sv[:22], TRAIT_NAME[tc], mode, n, len(indep), len(allc),
                             f"{b0:.4g}", f"{p0:.3g}", f"{bc:.4g}", f"{pc:.3g}", f"{ac:.1f}",
                             f"{bf:.4g}", f"{pf:.3g}", f"{af:.1f}", top["r2_with_predictor"], top["is_predictor_weight_snp"]])
            print(f"  {sv[:22]:24}{TRAIT_NAME[tc]:8} n={n} | M0 b={b0:+.4f} p={p0:.1e} | clean(+{len(indep)} indep) "
                  f"b={bc:+.4f} p={pc:.1e} at={ac:+.0f}% | full(+{len(allc)}) b={bf:+.4f} p={pf:.1e} at={af:+.0f}%", flush=True)
        with open(f"{OUT}/leads_conditional_{mode}.tsv", "w") as fo:
            fo.write(HDR)
            for r in out_rows:
                fo.write("\t".join(str(x) for x in r) + "\n")
        print(f"  wrote leads_conditional_{mode}.tsv", flush=True)


if __name__ == "__main__":
    main()
