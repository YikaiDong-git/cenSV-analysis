#!/usr/bin/env python3
"""Step 55 — compare the full ext2 rescore against the published 6,513-predictor run."""
import os, sys, csv, json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phewide_ext"))
import step42_gate2_headline_rescore as G

RV = G.RV
PUB = f"{RV}/step30_phewide/disease_matrices_firth/disease_hits_final.tsv"
EXT = f"{RV}/step30_phewide_ext2/disease_matrices_firth/disease_hits_final.tsv"
GATE2 = f"{RV}/step42_gate2/results/gate2_bh_recount_persd.tsv"
QC = f"{RV}/step28_phewide_ext2/dosage_qc_ext2.tsv"
OUT = f"{RV}/step30_phewide_ext2"
CMP = f"{OUT}/comparison/allancestry"


def fnum(x, d=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def load(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def main():
    os.makedirs(CMP, exist_ok=True)
    pub = load(PUB); ext = load(EXT)
    pub_sig = {(r["sv_locus"], r["phecode"]): r for r in pub if fnum(r["q_bh_tier"], 1) < 0.05}
    ext_by = {(r["sv_locus"], r["phecode"]): r for r in ext}
    ext_sig = {k: r for k, r in ext_by.items() if fnum(r["q_bh_tier"], 1) < 0.05}
    print(f"published hits rows {len(pub):,} ({len(pub_sig)} q<0.05)   "
          f"ext2 hits rows {len(ext):,} ({len(ext_sig)} q<0.05)", flush=True)

    scored = {r["sv_locus"] for r in csv.DictReader(open(QC), delimiter="\t") if r["scored"] == "1"}
    g2 = {(r["sv_locus"], r["phecode"]): r for r in load(GATE2)} if os.path.exists(GATE2) else {}

    rows = []
    for k, pr in sorted(pub_sig.items(), key=lambda kv: fnum(kv[1]["p_final"])):
        sv, ph = k
        er = ext_by.get(k)
        in_pgen = sv in scored
        q_e = fnum(er["q_bh_tier"], float("nan")) if er else float("nan")
        p_e = fnum(er["p_final"], float("nan")) if er else float("nan")
        surv = int(np.isfinite(q_e) and q_e < 0.05)
        gr = g2.get(k)
        rows.append(dict(
            sv_locus=sv, phecode=ph, label=pr["label"], cases=pr["cases"],
            p_published=fnum(pr["p_final"]), q_published=fnum(pr["q_bh_tier"]),
            p_ext2=p_e, q_ext2=q_e,
            in_ext2_hits=int(er is not None), sv_scored_in_ext2=int(in_pgen),
            survives_full_rescore=surv,
            gate2_predicted=int(gr["survives_exact_bh"]) if gr else -1,
            agrees_with_gate2=(int(surv == int(gr["survives_exact_bh"])) if gr else -1)))
    n_surv = sum(r["survives_full_rescore"] for r in rows)
    n_g2 = sum(1 for r in rows if r["gate2_predicted"] == 1)
    n_agree = sum(1 for r in rows if r["agrees_with_gate2"] == 1)
    print(f"\n1. FATE OF THE PUBLISHED {len(rows)} q<0.05 CELLS", flush=True)
    print(f"   survive the FULL rescore: {n_surv}/{len(rows)}   "
          f"(Gate 2's held-fixed estimate was {n_g2}/{len(rows)}; agreement {n_agree}/{len(rows)})",
          flush=True)
    for r in rows:
        print(f"   {'KEEP' if r['survives_full_rescore'] else 'LOST'}  "
              f"p {r['p_published']:.3g} -> {r['p_ext2']:.3g}   q {r['q_published']:.3g} -> "
              f"{r['q_ext2']:.3g}   {r['sv_locus'][:34]:34s} {r['phecode']}", flush=True)

    new = {k: r for k, r in ext_sig.items() if k not in pub_sig}
    from collections import Counter
    by_tier = Counter(r["tier"] for r in new.values())
    old_pred = {k: r for k, r in new.items() if r["tier"] in ("signal", "null_ext")}
    new_pred = {k: r for k, r in new.items() if r["tier"] not in ("signal", "null_ext")}
    print(f"\n2. NEW q<0.05 CELLS: {len(new)}", flush=True)
    print(f"   on already-published predictors (corrected weights): {len(old_pred)}", flush=True)
    print(f"   on never-before-tested predictors (new coverage)   : {len(new_pred)}", flush=True)
    print(f"   by tier: {dict(by_tier)}", flush=True)
    print("   strongest 15 new cells:", flush=True)
    for k, r in sorted(new.items(), key=lambda kv: fnum(kv[1]["p_final"]))[:15]:
        print(f"     p={fnum(r['p_final']):.3g} q={fnum(r['q_bh_tier']):.3g} [{r['tier']:9s}] "
              f"{r['sv_locus'][:36]:36s} {r['phecode']:9s} {r['label'][:34]}", flush=True)

    print("\n3. PER-TIER", flush=True)
    tiers = sorted({r["tier"] for r in ext})
    tier_rows = []
    for t in tiers:
        sig_t = [r for r in ext_sig.values() if r["tier"] == t]
        cut = max((fnum(r["p_final"]) for r in sig_t), default=float("nan"))
        tier_rows.append(dict(tier=t, n_hit_rows=sum(1 for r in ext if r["tier"] == t),
                              n_q05=len(sig_t), bh_cutoff_p=cut))
        print(f"   {t:11s} hit-rows {tier_rows[-1]['n_hit_rows']:>6,}  q<0.05 {len(sig_t):>5}  "
              f"largest passing p {cut:.4g}", flush=True)

    HEAD = sorted({k[0] for k in pub_sig})
    print("\n4. THE 8 GATE-2 HEADLINE SVs (published q<0.05 cells kept / total)", flush=True)
    head_rows = []
    for sv in HEAD:
        tot = sum(1 for k in pub_sig if k[0] == sv)
        kept = sum(1 for r in rows if r["sv_locus"] == sv and r["survives_full_rescore"])
        allnew = sum(1 for k, r in ext_sig.items() if k[0] == sv)
        head_rows.append(dict(sv_locus=sv, published_cells=tot, kept=kept, ext2_total_cells=allnew))
        print(f"   {kept}/{tot} kept  ({allnew} q<0.05 cells on this SV in the rescore)  {sv}", flush=True)

    def w(path, data):
        ks = list(data[0].keys())
        with open(G.guard(path), "w") as f:
            f.write("\t".join(ks) + "\n")
            for r in data:
                f.write("\t".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k])
                                  for k in ks) + "\n")
    w(f"{CMP}/ext2_vs_published_published28.tsv", rows)
    w(f"{CMP}/ext2_vs_published_pertier.tsv", tier_rows)
    w(f"{CMP}/ext2_vs_published_headline_svs.tsv", head_rows)
    if new:
        w(f"{CMP}/ext2_new_hits.tsv",
          [dict(sv_locus=k[0], phecode=k[1], tier=r["tier"], label=r["label"], cases=r["cases"],
                p_final=fnum(r["p_final"]), q_bh_tier=fnum(r["q_bh_tier"]),
                beta_final_logOR=fnum(r["beta_final_logOR"]), p_source=r["p_source"],
                separation=r["separation"],
                on_new_predictor=int(r["tier"] not in ("signal", "null_ext")))
           for k, r in sorted(new.items(), key=lambda kv: fnum(kv[1]["p_final"]))])
    with open(G.guard(f"{CMP}/step55_summary.json"), "w") as f:
        json.dump(dict(published_q05=len(pub_sig), ext2_q05=len(ext_sig),
                       published28_surviving=n_surv, gate2_predicted=n_g2, agreement=n_agree,
                       new_cells=len(new), new_on_published_predictors=len(old_pred),
                       new_on_new_predictors=len(new_pred), by_tier=dict(by_tier),
                       per_tier=tier_rows), f, indent=1)
    print(f"\nSTEP55_DONE -> {OUT}/ext2_vs_published_*.tsv", flush=True)


if __name__ == "__main__":
    main()
