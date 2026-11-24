#!/usr/bin/env python3
"""Step 59 (ext2 arm) — per-ancestry: ext2 vs published."""
import os, sys, csv, json
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S

RV = S.RV
PUB_PA = f"{RV}/step30_phewide/disease_matrices_firth/per_ancestry"
EXT_PA = f"{S.OUTA}/disease_matrices_firth/per_ancestry"
PUB_ALL = f"{RV}/step30_phewide/disease_matrices_firth/disease_hits_final.tsv"
EXT_ALL = f"{S.OUTA}/disease_matrices_firth/disease_hits_final.tsv"
NEWCELLS = f"{S.OUTA}/comparison/allancestry/ext2_new_hits.tsv"
PUB28 = f"{S.OUTA}/comparison/allancestry/ext2_vs_published_published28.tsv"
CMP = f"{S.OUTA}/comparison/perancestry"
ANCS = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
WELL = ["EUR", "AFR", "AMR"]


def load_pa(root):
    sig = {}; allrows = {}
    for anc in ANCS:
        p = f"{root}/disease_hits_final_{anc}.tsv"
        if not os.path.exists(p):
            sig[anc] = {}; allrows[anc] = {}; continue
        rows = list(csv.DictReader(open(p), delimiter="\t"))
        allrows[anc] = {(r["sv_locus"], r["phecode"]): r for r in rows}
        sig[anc] = {(r["sv_locus"], r["phecode"]): r for r in rows if r["tier_label"] == "significant"}
    return sig, allrows


def classify(cell, sig):
    hit = [a for a in WELL if cell in sig.get(a, {})]
    if not hit:
        return "none", hit
    if hit == ["AMR"]:
        return "AMR_only", hit
    if hit == ["EUR"]:
        return "EUR_only", hit
    if hit == ["AFR"]:
        return "AFR_only", hit
    return "cross_ancestry", hit


def power_adjudication_unbiased(ref="AMR", others=("EUR", "AFR")):
    import numpy as np
    from scipy.stats import binomtest
    MAT = f"{S.OUTA}/disease_matrices/per_ancestry"
    hitf = f"{EXT_PA}/disease_hits_final_{ref}.tsv"
    if not (os.path.exists(f"{MAT}/disease_U_{ref}.npy") and os.path.exists(hitf)):
        print("(unbiased power adjudication skipped — missing score U/V or per-ancestry Firth hits)")
        return None
    svs = [l.split("\t")[0] for l in open(f"{MAT}/disease_svorder_{ref}.txt") if l.strip()]
    meta = list(csv.DictReader(open(f"{MAT}/disease_meta_{ref}.tsv"), delimiter="\t"))
    si = {s: i for i, s in enumerate(svs)}; pi = {m["phecode"]: j for j, m in enumerate(meta)}
    U = {}; V = {}
    for a in (ref,) + tuple(others):
        U[a] = np.load(f"{MAT}/disease_U_{a}.npy"); V[a] = np.load(f"{MAT}/disease_V_{a}.npy")
    sig = [r for r in csv.DictReader(open(hitf), delimiter="\t")
           if r["tier_label"] == "significant" and r.get("tier") != "ext2_chrX"]
    best = {}
    for r in sig:
        k = r["sv_locus"]
        if k not in best or float(r["p_final"]) < float(best[k]["p_final"]):
            best[k] = r
    print(f"\n=== UNBIASED EXCLUSIVITY TEST: cells selected ONLY on the {ref} per-ancestry Firth ===")
    print(f"  {len(sig)} autosomal {ref}-significant cells over {len(best)} distinct SVs "
          f"(selection never touches {'/'.join(others)} data)")
    out = {"ref": ref, "n_cells": len(sig), "n_svs": len(best)}
    for a in others:
        res = {}
        for tag, rows in (("cell", sig), ("per_sv", list(best.values()))):
            tot = c = hi = hic = 0; obs = []
            for r in rows:
                i = si.get(r["sv_locus"]); j = pi.get(r["phecode"])
                if i is None or j is None:
                    continue
                va, ua = float(V[ref][i, j]), float(U[ref][i, j])
                v, u = float(V[a][i, j]), float(U[a][i, j])
                if not (np.isfinite(va) and va > 0 and np.isfinite(v) and v > 0):
                    continue
                ba = ua / va; same = int(np.sign(u / v) == np.sign(ba)); exp = ba * ba * v
                tot += 1; c += same; obs.append(u * u / v)
                if exp >= 8.0:
                    hi += 1; hic += same
            res[tag] = dict(n=tot, conc=c, pct=100 * c / max(tot, 1),
                            p=float(binomtest(c, tot, 0.5).pvalue) if tot else None,
                            powered=hi, powered_conc=hic,
                            powered_pct=100 * hic / max(hi, 1),
                            powered_p=float(binomtest(hic, hi, 0.5).pvalue) if hi else None,
                            median_obs_chi2=float(np.median(obs)) if obs else None)
        cl, sv = res["cell"], res["per_sv"]
        print(f"  {a}: cell-level  {cl['conc']}/{cl['n']} ({cl['pct']:.0f}%) p={cl['p']:.3f} | "
              f"powered {cl['powered_conc']}/{cl['powered']} ({cl['powered_pct']:.0f}%) "
              f"p={cl['powered_p']:.3f} | median obs chi2 {cl['median_obs_chi2']:.3f} (null 0.455)")
        print(f"  {a}: PER-SV      {sv['conc']}/{sv['n']} ({sv['pct']:.0f}%) p={sv['p']:.3f} | "
              f"powered {sv['powered_conc']}/{sv['powered']} ({sv['powered_pct']:.0f}%) "
              f"p={sv['powered_p']:.3f}   <-- quote THIS one")
        out[a] = res
    print("  Per-SV sign concordance at chance + observed chi2 at the null => the effect is ABSENT in those")
    print("  strata, not attenuated (attenuation preserves sign). Stratum-PRIVATE: AMR-private real biology")
    print("  OR an AMR-specific confound. step60/step61 (within-AMR PCs) separate those two.")
    return out


def power_adjudication(sv_prefix="chr16_", ph_prefix="MB_28", ref="AMR", others=("EUR", "AFR")):
    import numpy as np
    from scipy.stats import chi2 as chi2dist
    MAT = f"{S.OUTA}/disease_matrices/per_ancestry"
    if not os.path.exists(f"{MAT}/disease_U_{ref}.npy"):
        print("(power adjudication skipped — no score U/V matrices)"); return None
    svs = [l.split("\t")[0] for l in open(f"{MAT}/disease_svorder_{ref}.txt") if l.strip()]
    meta = list(csv.DictReader(open(f"{MAT}/disease_meta_{ref}.tsv"), delimiter="\t"))
    si = {s: i for i, s in enumerate(svs)}; pi = {m["phecode"]: j for j, m in enumerate(meta)}
    U = {}; V = {}
    for a in (ref,) + tuple(others):
        U[a] = np.load(f"{MAT}/disease_U_{a}.npy"); V[a] = np.load(f"{MAT}/disease_V_{a}.npy")
    hitf = f"{S.OUTA}/per_ancestry/phewide_disease_full_{ref}.tsv"
    hits = [r for r in csv.DictReader(open(hitf), delimiter="\t")
            if r["q_bh_tier"] and float(r["q_bh_tier"]) < 0.05
            and r["sv_locus"].startswith(sv_prefix) and r["phecode"].startswith(ph_prefix)]
    recs = []
    for h in hits:
        i = si.get(h["sv_locus"]); j = pi.get(h["phecode"])
        if i is None or j is None:
            continue
        ur, vr = float(U[ref][i, j]), float(V[ref][i, j])
        if not (np.isfinite(vr) and vr > 0):
            continue
        b_ref = ur / vr
        rec = {"b_ref": b_ref}
        for a in others:
            u, v = float(U[a][i, j]), float(V[a][i, j])
            if np.isfinite(v) and v > 0:
                rec[a] = (b_ref * b_ref * v, u * u / v, u / v)
            else:
                rec[a] = (np.nan, np.nan, np.nan)
        recs.append(rec)
    print(f"\n=== POWER ADJUDICATION: is {ref}-exclusivity real, or just power? "
          f"({len(recs)} {ref}-significant {sv_prefix}* x {ph_prefix}* cells) ===")
    out = {"ref": ref, "n_cells": len(recs)}
    for a in others:
        e = np.array([r[a][0] for r in recs]); o = np.array([r[a][1] for r in recs])
        b = np.array([r[a][2] for r in recs]); br = np.array([r["b_ref"] for r in recs])
        fin = np.isfinite(e) & np.isfinite(o)
        pw8 = fin & (e >= 8.0)
        nullish = int((pw8 & (o < 3.841459)).sum())
        conc = int((np.sign(b[pw8]) == np.sign(br[pw8])).sum())
        print(f"  {a}: median exp chi2 {np.nanmedian(e[fin]):7.2f} | strongly powered (>=8) "
              f"{int(pw8.sum())}/{int(fin.sum())} | of those observed-null {nullish}/{int(pw8.sum())} | "
              f"median observed chi2 {np.nanmedian(o[pw8]):.3f} (null 0.455) | "
              f"sign-concordant {conc}/{int(pw8.sum())} ({100*conc/max(1,int(pw8.sum())):.0f}%)")
        out[a] = dict(median_exp_chi2=float(np.nanmedian(e[fin])), strongly_powered=int(pw8.sum()),
                      testable=int(fin.sum()), observed_null=nullish,
                      median_obs_chi2=float(np.nanmedian(o[pw8])), sign_concordant=conc)
    print("  Sign concordance ~50% is incompatible with differential-imputation attenuation (which preserves")
    print("  sign) -> the effect is ABSENT in those strata, not merely weakened. That is a stratum-PRIVATE")
    print("  signature: AMR-private real biology OR an AMR-specific confound; step61 separates those.")
    return out


def main():
    os.makedirs(CMP, exist_ok=True)
    pub_sig, _pub_all = load_pa(PUB_PA)
    ext_sig, ext_all = load_pa(EXT_PA)
    have = [a for a in ANCS if ext_sig.get(a)]
    print(f"ext2 per-ancestry files present for: {have}\n", flush=True)

    rows = []
    print("=== per-ancestry significant cells: PUBLISHED (old wts, 6,513, autosomal) vs EXT2 (11,129) ===")
    print(f"{'anc':<6}{'pub_sig':>9}{'ext2_auto':>11}{'ext2_chrX':>10}{'ext2_all':>9}"
          f"{'pub_c16xMB28':>14}{'ext2_c16xMB28':>15}")
    for anc in ANCS:
        ps, es = pub_sig.get(anc, {}), ext_sig.get(anc, {})
        ex_x = {k: v for k, v in es.items() if v.get("tier") == "ext2_chrX"}
        ex_a = {k: v for k, v in es.items() if v.get("tier") != "ext2_chrX"}
        pc16 = [k for k in ps if k[0].startswith("chr16_")]
        ec16 = [k for k in ex_a if k[0].startswith("chr16_")]
        pmb = [k for k in pc16 if k[1].startswith("MB_28")]
        emb = [k for k in ec16 if k[1].startswith("MB_28")]
        rows.append(dict(ancestry=anc, published_sig=len(ps), ext2_sig_autosomal=len(ex_a),
                         ext2_sig_chrX=len(ex_x), ext2_sig_total=len(es),
                         published_chr16=len(pc16), ext2_chr16_autosomal=len(ec16),
                         published_chr16_MB28=len(pmb), ext2_chr16_MB28=len(emb)))
        print(f"{anc:<6}{len(ps):>9}{len(ex_a):>11}{len(ex_x):>10}{len(es):>9}{len(pmb):>14}{len(emb):>15}")
    print("  (ext2_chrX is a NEW tier with no published counterpart and no ploidy handling — "
          "read the autosomal column for the like-for-like comparison)")
    with open(f"{CMP}/perancestry_counts_ext2.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

    ext_all_sig = [r for r in csv.DictReader(open(EXT_ALL), delimiter="\t")
                   if r["tier_label"] == "significant"]
    newset = set()
    onnew = {}
    if os.path.exists(NEWCELLS):
        for r in csv.DictReader(open(NEWCELLS), delimiter="\t"):
            newset.add((r["sv_locus"], r["phecode"])); onnew[(r["sv_locus"], r["phecode"])] = r.get("on_new_predictor", "")
    out = []
    tally = defaultdict(Counter)
    for r in ext_all_sig:
        cell = (r["sv_locus"], r["phecode"])
        lab, hit = classify(cell, ext_sig)
        grp = "new" if cell in newset else "published_predictor"
        tally[grp][lab] += 1
        rec = dict(sv_locus=cell[0], phecode=cell[1], label=r.get("label", ""), tier=r.get("tier", ""),
                   pooled_p_final=r.get("p_final", ""), pooled_q=r.get("q_bh_tier", ""),
                   is_new_cell=int(cell in newset), on_new_predictor=onnew.get(cell, ""),
                   attribution=lab, sig_well_strata=",".join(hit))
        for a in ANCS:
            rr = ext_sig.get(a, {}).get(cell) or ext_all.get(a, {}).get(cell)
            rec[f"{a}_p"] = rr["p_final"] if rr else ""
            rec[f"{a}_beta"] = rr["beta_final_logOR"] if rr else ""
            rec[f"{a}_sig"] = int(cell in ext_sig.get(a, {}))
        out.append(rec)
    with open(f"{CMP}/perancestry_new_cell_attribution_ext2.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)

    print(f"\n=== ancestry attribution of the {len(ext_all_sig)} ext2 all-ancestry SIGNIFICANT cells ===")
    labs = ["cross_ancestry", "AMR_only", "EUR_only", "AFR_only", "none"]
    print(f"  {'group':<22}{'n':>5}" + "".join(f"{l:>16}" for l in labs))
    for grp in ("published_predictor", "new"):
        tot = sum(tally[grp].values())
        print(f"  {grp:<22}{tot:>5}" + "".join(f"{tally[grp].get(l,0):>16}" for l in labs))

    fate = []
    if os.path.exists(PUB28):
        for r in csv.DictReader(open(PUB28), delimiter="\t"):
            cell = (r["sv_locus"], r["phecode"])
            plab, phit = classify(cell, pub_sig)
            elab, ehit = classify(cell, ext_sig)
            fate.append(dict(sv_locus=cell[0], phecode=cell[1], label=r.get("label", ""),
                             survives_full_rescore=r.get("survives_full_rescore", ""),
                             published_attribution=plab, published_sig_strata=",".join(phit),
                             ext2_attribution=elab, ext2_sig_strata=",".join(ehit),
                             attribution_changed=int(plab != elab)))
        with open(f"{CMP}/perancestry_published28_fate_ext2.tsv", "w") as f:
            w = csv.DictWriter(f, fieldnames=list(fate[0]), delimiter="\t"); w.writeheader(); w.writerows(fate)
        ch = sum(r["attribution_changed"] for r in fate)
        print(f"\n=== published 28 cells: per-ancestry attribution ===")
        print(f"  attribution UNCHANGED vs published: {len(fate)-ch}/{len(fate)};  changed: {ch}")
        for r in fate:
            if r["attribution_changed"]:
                print(f"    {r['sv_locus'][:26]:<28}{r['phecode']:<10} "
                      f"{r['published_attribution']} -> {r['ext2_attribution']}")

    c16 = [r for r in out if r["sv_locus"].startswith("chr16_") and r["phecode"].startswith("MB_28")]
    c16_amr = [r for r in c16 if r["attribution"] == "AMR_only"]
    c16_cross = [r for r in c16 if r["attribution"] == "cross_ancestry"]
    eur_sig_any = [r for r in out if r[f"EUR_sig"] == 1]
    afr_sig_any = [r for r in out if r[f"AFR_sig"] == 1]
    print(f"\n=== Q1: is the chr16 -> substance-use cluster still AMR-RESTRICTED under ext2 weights? ===")
    print(f"  chr16 x MB_28x cells significant all-ancestry : {len(c16)}")
    print(f"    AMR-only                                    : {len(c16_amr)}")
    print(f"    cross-ancestry (>=2 WELL strata)            : {len(c16_cross)}")
    print(f"    EUR-significant anywhere (any cell)         : {len(eur_sig_any)}")
    print(f"    AFR-significant anywhere (any cell)         : {len(afr_sig_any)}")
    print(f"  published comparison: EUR 1 / AFR 0 / AMR 150 sig; 122 chr16xMB_28 AMR cells, "
          f"122/122 survived within-AMR PCs")

    newrows = [r for r in out if r["is_new_cell"]]
    ntal = Counter(r["attribution"] for r in newrows)
    print(f"\n=== Q2: ancestry of the {len(newrows)} NEW ext2 cells ===")
    for l in labs:
        if ntal.get(l):
            print(f"    {l:<18}{ntal[l]:>5}")

    power = power_adjudication()
    power_ub = power_adjudication_unbiased()

    summary = dict(
        power_adjudication_unbiased=power_ub,
        power_adjudication_scorelevel=power,
        per_ancestry_counts=rows,
        ext2_all_sig=len(ext_all_sig), new_cells=len(newrows),
        attribution_published_predictor=dict(tally["published_predictor"]),
        attribution_new=dict(tally["new"]),
        chr16_MB28_sig=len(c16), chr16_MB28_AMR_only=len(c16_amr), chr16_MB28_cross=len(c16_cross),
        eur_sig_cells=len(eur_sig_any), afr_sig_cells=len(afr_sig_any),
        published28_attribution_changed=sum(r["attribution_changed"] for r in fate) if fate else None,
    )
    json.dump(summary, open(f"{CMP}/step59_summary.json", "w"), indent=2)
    print(f"\nwrote {CMP}/perancestry_*_ext2.tsv + step59_summary.json")
    print("EXT2_PERANC_COMPARE_DONE", flush=True)


if __name__ == "__main__":
    main()
