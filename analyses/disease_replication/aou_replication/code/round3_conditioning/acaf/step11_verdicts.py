#!/usr/bin/env python3
"""Aggregate the SNP-conditioning cells (step10) into one verdict per (SV, trait)."""
import csv, os, statistics as st

R = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results/step10_cond")
ORDER = [("chr3_sv20", "weight"), ("chr4_sv349", "height"), ("chr11_sv263", "BMI"),
         ("chr4_sv398", "MPV_30100"), ("chr4_sv398", "PDW_30110")]
BH = {"weight", "height", "BMI"}


def load(mode):
    f = f"{R}/conditioning_{mode}.tsv"
    return list(csv.DictReader(open(f), delimiter="\t")) if os.path.exists(f) else []


def num(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def main():
    eur, allm = load("eur"), load("all")
    out = []
    for svshort, trait in ORDER:
        prim = [r for r in eur if r["sv"].startswith(svshort) and r["trait"] == trait
                and r["flank_kb"] == "0" and int(r["k_pcs"]) <= 40]
        if not prim:
            continue
        p0 = num(prim[0]["p_m0"]); b0 = num(prim[0]["beta_m0"])
        pm1 = [num(r["p_m1"]) for r in prim]
        att = [num(r["pct_atten"]) for r in prim]
        frac_sig = sum(p < 0.05 for p in pm1) / len(pm1)
        att_med = st.median(att)
        if frac_sig >= 0.60 and att_med < 35:
            verdict = "SURVIVES"
        elif frac_sig <= 0.34 and att_med > 40:
            verdict = "COLLAPSES"
        else:
            verdict = "PARTIAL"
        ac = [r for r in allm if r["sv"].startswith(svshort) and r["trait"] == trait
              and r["flank_kb"] == "0" and r["k_pcs"] == "20"]
        ap = num(ac[0]["p_m1"]) if ac else None
        aa = num(ac[0]["pct_atten"]) if ac else None
        fc = [r for r in eur if r["sv"].startswith(svshort) and r["trait"] == trait
              and r["flank_kb"] == "10" and r["k_pcs"] == "20"]
        fp = num(fc[0]["p_m1"]) if fc else None
        fa = num(fc[0]["pct_atten"]) if fc else None
        out.append(dict(sv=svshort, trait=trait, bh=("BH" if trait in BH else "nominal"),
                        p_m0=p0, beta_m0=b0, n_prim=len(prim), frac_sig=frac_sig,
                        atten_med=att_med, p_m1_min=min(pm1), p_m1_max=max(pm1),
                        all_p_m1=ap, all_atten=aa, flank_p_m1=fp, flank_atten=fa, verdict=verdict))
    cols = ["sv", "trait", "bh", "verdict", "beta_m0", "p_m0", "n_prim", "frac_sig",
            "atten_med", "p_m1_min", "p_m1_max", "all_p_m1", "all_atten", "flank_p_m1", "flank_atten"]
    with open(f"{R}/conditioning_verdicts.tsv", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in out:
            fo.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    print(f"{'SV':12}{'trait':11}{'tier':8}{'verdict':10}{'M0 p':>9}  {'sig cells':>9}  {'med atten':>9}  "
          f"{'cond p range':>18}  {'all p':>8}  {'flank p':>8}")
    for r in out:
        print(f"{r['sv']:12}{r['trait']:11}{r['bh']:8}{r['verdict']:10}{r['p_m0']:9.1e}  "
              f"{r['frac_sig']*r['n_prim']:.0f}/{r['n_prim']:<7}  {r['atten_med']:+8.0f}%  "
              f"{r['p_m1_min']:.1e}-{r['p_m1_max']:.1e}  {r['all_p_m1']:8.1e}  {r['flank_p_m1']:8.1e}")
    print(f"\nwrote {R}/conditioning_verdicts.tsv")


if __name__ == "__main__":
    main()
