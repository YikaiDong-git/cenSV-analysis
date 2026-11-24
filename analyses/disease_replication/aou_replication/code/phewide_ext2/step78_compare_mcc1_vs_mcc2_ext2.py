#!/usr/bin/env python3
"""Step 78 (ext2 arm) — MCC=1 vs MCC=2 disease-PheWAS sensitivity comparison, both on ext2 weights."""
import os, sys, csv, json, hashlib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S

RV = S.RV
CARRY = S.CARRY
MCC2 = f"{RV}/step30_phewide_ext2"
MCC1 = f"{RV}/step30_phewide_mcc1_ext2"
MCC1_INP = f"{RV}/inputs_mcc1_ext2"
OUT = f"{MCC2}/comparison/mcc1"

C16 = "chr16_"
SUBST = "MB_28"

_log_lines = []


def log(msg=""):
    print(msg, flush=True)
    _log_lines.append(str(msg))


def md5_file(path, buf=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(buf), b""):
            h.update(b)
    return h.hexdigest()


def load_arm(root):
    mat = f"{root}/disease_matrices"
    fir = f"{root}/disease_matrices_firth"
    svloc, tiers = [], []
    for ln in open(f"{mat}/disease_svorder.txt"):
        if ln.strip():
            a = ln.rstrip("\n").split("\t")
            svloc.append(a[0]); tiers.append(a[1] if len(a) > 1 else "")
    meta = list(csv.DictReader(open(f"{mat}/disease_meta.tsv"), delimiter="\t"))
    ph_row = {m["phecode"]: m for m in meta}
    hits = {}
    for r in csv.DictReader(open(f"{fir}/disease_hits_final.tsv"), delimiter="\t"):
        hits[(r["sv_locus"], r["phecode"])] = r
    return dict(
        root=root, svloc=svloc, tiers=tiers,
        svidx={s: i for i, s in enumerate(svloc)},
        phidx={m["phecode"]: j for j, m in enumerate(meta)},
        ph_row=ph_row, meta=meta, hits=hits,
        Pv=np.load(f"{fir}/disease_Pv_final.npy", mmap_mode="r"),
        Be=np.load(f"{fir}/disease_Be_final.npy", mmap_mode="r"),
    )


def cell_status(arm, sv, ph):
    row = arm["ph_row"].get(ph)
    if row is None or row.get("tested") != "1":
        return ("not_tested", np.nan, np.nan, np.nan, 0)
    h = arm["hits"].get((sv, ph))
    i = arm["svidx"].get(sv); j = arm["phidx"].get(ph)
    p = float(arm["Pv"][i, j]) if (i is not None and j is not None) else np.nan
    b = float(arm["Be"][i, j]) if (i is not None and j is not None) else np.nan
    if h is None:
        return ("sub-threshold", np.nan, p, b, 0)
    q = float(h["q_bh_tier"]) if h["q_bh_tier"] else np.nan
    if not np.isfinite(p) or p == 0.0:
        try:
            p = float(h["p_final"])
        except ValueError:
            pass
    if not np.isfinite(b):
        try:
            b = float(h["beta_final_logOR"])
        except ValueError:
            pass
    return (h["tier_label"], q, p, b, int(h["separation"] or 0))


def is_sig(lbl):
    return lbl == "significant"


def is_ind_row(r):
    try:
        q = float(r["q_bh_tier"])
    except (ValueError, KeyError, TypeError):
        return False
    return 0.05 <= q < 0.10


def gates(a2, a1):
    log("=" * 96)
    log("GATES (fail-closed)")
    log("=" * 96)
    fail = []

    m2 = md5_file(f"{CARRY}/inputs/covar.txt")
    m1 = md5_file(f"{MCC1_INP}/covar.txt")
    ok = (m1 == m2)
    log(f"G1 covar.txt md5   MCC2={m2}  MCC1={m1}  -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        fail.append("G1 covariate/sample set differs between arms")

    ok = (a2["svloc"] == a1["svloc"]) and (a2["tiers"] == a1["tiers"])
    log(f"G2 SV order+tiers  MCC2 n={len(a2['svloc'])}  MCC1 n={len(a1['svloc'])}  "
        f"identical={ok} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        fail.append("G2 predictor set/order/tier differs between arms")

    p2 = {m["phecode"]: int(m["cases"]) for m in a2["meta"] if m.get("tested") == "1"}
    p1 = {m["phecode"]: int(m["cases"]) for m in a1["meta"] if m.get("tested") == "1"}
    missing = sorted(set(p2) - set(p1))
    shared = sorted(set(p2) & set(p1))
    viol = [ph for ph in shared if p1[ph] < p2[ph]]
    ratio = np.array([p1[ph] / p2[ph] for ph in shared], float)
    ok = (not missing) and (not viol)
    log(f"G3 phecodes        MCC2 tested={len(p2)}  MCC1 tested={len(p1)}  "
        f"MCC1-only=+{len(set(p1) - set(p2))}  MCC2-not-in-MCC1={len(missing)}")
    log(f"   monotonicity    violations={len(viol)}  inflation median={np.median(ratio):.3f}x  "
        f"range {ratio.min():.3f}-{ratio.max():.3f}  -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        fail.append(f"G3 non-monotone relaxation ({len(missing)} lost, {len(viol)} shrank)")

    for nm, a in (("MCC2", a2), ("MCC1", a1)):
        blank = sum(1 for t in a["tiers"] if not t)
        cnt = {}
        for t in a["tiers"]:
            cnt[t] = cnt.get(t, 0) + 1
        log(f"G4 {nm} tiers      {cnt}  blank={blank} -> {'PASS' if blank == 0 else 'FAIL'}")
        if blank:
            fail.append(f"G4 {nm} has {blank} SVs with no tier (fail-open: excluded from FDR)")

    for nm, a in (("MCC2", a2), ("MCC1", a1)):
        exp = (len(a["svloc"]), len(a["meta"]))
        ok = (tuple(a["Pv"].shape) == exp and tuple(a["Be"].shape) == exp)
        log(f"G5 {nm} matrices   Pv{tuple(a['Pv'].shape)} Be{tuple(a['Be'].shape)} vs meta{exp} "
            f"-> {'PASS' if ok else 'FAIL'}")
        if not ok:
            fail.append(f"G5 {nm} matrix shape != metadata")

    if fail:
        log("")
        for f in fail:
            log(f"  ABORT: {f}")
        raise SystemExit("GATES FAILED — comparison not written.")
    log("all gates PASS\n")
    return dict(covar_md5=m2, n_sv=len(a2["svloc"]),
                n_ph_mcc2=len(p2), n_ph_mcc1=len(p1),
                n_ph_mcc1_only=len(set(p1) - set(p2)),
                infl_median=float(np.median(ratio)),
                infl_min=float(ratio.min()), infl_max=float(ratio.max()))


def main():
    os.makedirs(OUT, exist_ok=True)
    for d, nm in ((MCC2, "MCC2"), (MCC1, "MCC1")):
        for rel in ("disease_matrices_firth/disease_hits_final.tsv",
                    "disease_matrices_firth/disease_Pv_final.npy",
                    "disease_matrices_firth/disease_Be_final.npy",
                    "disease_matrices/disease_meta.tsv",
                    "disease_matrices/disease_svorder.txt"):
            f = f"{d}/{rel}"
            if not os.path.exists(f):
                raise SystemExit(f"ABORT: {nm} arm incomplete — missing {f}")
    a2, a1 = load_arm(MCC2), load_arm(MCC1)
    g = gates(a2, a1)

    sig2 = {k for k, r in a2["hits"].items() if is_sig(r["tier_label"])}
    sig1 = {k for k, r in a1["hits"].items() if is_sig(r["tier_label"])}
    ind2 = {k for k, r in a2["hits"].items() if is_ind_row(r)}
    ind1 = {k for k, r in a1["hits"].items() if is_ind_row(r)}

    union = sorted(sig2 | sig1)
    rows = []
    for sv, ph in union:
        l2, q2, p2, b2, s2 = cell_status(a2, sv, ph)
        l1, q1, p1, b1, s1 = cell_status(a1, sv, ph)
        in2, in1 = is_sig(l2), is_sig(l1)
        fate = "HELD" if (in2 and in1) else ("DROPPED" if in2 else "APPEARED")
        c2 = a2["ph_row"].get(ph, {}).get("cases", "")
        c1 = a1["ph_row"].get(ph, {}).get("cases", "")
        rows.append(dict(
            sv_locus=sv, tier=a2["tiers"][a2["svidx"][sv]], phecode=ph,
            label=(a2["ph_row"].get(ph) or a1["ph_row"].get(ph, {})).get("label", ""),
            fate=fate,
            mcc2_status=l2, mcc2_q=q2, mcc2_p=p2, mcc2_beta=b2, mcc2_cases=c2, mcc2_sep=s2,
            mcc1_status=l1, mcc1_q=q1, mcc1_p=p1, mcc1_beta=b1, mcc1_cases=c1, mcc1_sep=s1,
            case_ratio=(int(c1) / int(c2)) if (c1 and c2 and int(c2)) else np.nan,
            logp_shift=(np.log10(p2) - np.log10(p1)) if (np.isfinite(p2) and np.isfinite(p1)
                                                         and p2 > 0 and p1 > 0) else np.nan,
            beta_sign_same=int(np.sign(b2) == np.sign(b1)) if (np.isfinite(b2) and np.isfinite(b1)) else "",
            chr16_subst=int(sv.startswith(C16) and ph.startswith(SUBST)),
        ))
    cols = list(rows[0].keys()) if rows else []
    with open(f"{OUT}/mcc1_vs_mcc2_cells.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["fate"], r["mcc2_p"] if np.isfinite(r["mcc2_p"]) else 1.0)):
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) and np.isfinite(v)
                            else ("" if isinstance(v, float) else v)) for k, v in r.items()})

    held = [r for r in rows if r["fate"] == "HELD"]
    dropped = [r for r in rows if r["fate"] == "DROPPED"]
    appeared = [r for r in rows if r["fate"] == "APPEARED"]

    log("=" * 96)
    log("HEADLINE")
    log("=" * 96)
    log(f"{'metric':<58}{'MCC=2 (primary)':>18}{'MCC=1 (this)':>18}")
    log(f"{'significant cells (q<0.05, Firth-recal)':<58}{len(sig2):>18}{len(sig1):>18}")
    log(f"{'indicative cells (q 0.05-0.1)':<58}{len(ind2):>18}{len(ind1):>18}")
    log(f"{'distinct SVs with >=1 significant cell':<58}"
        f"{len({k[0] for k in sig2}):>18}{len({k[0] for k in sig1}):>18}")
    log(f"{'distinct phecodes with >=1 significant cell':<58}"
        f"{len({k[1] for k in sig2}):>18}{len({k[1] for k in sig1}):>18}")
    log(f"{'chr16 x MB_28* significant cells':<58}"
        f"{len([k for k in sig2 if k[0].startswith(C16) and k[1].startswith(SUBST)]):>18}"
        f"{len([k for k in sig1 if k[0].startswith(C16) and k[1].startswith(SUBST)]):>18}")
    from collections import Counter
    t2 = Counter(a2["tiers"][a2["svidx"][k[0]]] for k in sig2)
    t1 = Counter(a1["tiers"][a1["svidx"][k[0]]] for k in sig1)
    log("")
    log(f"{'tier composition of significant cells':<58}{str(dict(t2)):>18}")
    log(f"{'  (MCC=1)':<58}{str(dict(t1)):>18}")
    log(f"{'NEGATIVE-CONTROL (null_ext) significant cells':<58}"
        f"{t2.get('null_ext', 0):>18}{t1.get('null_ext', 0):>18}   <- FPR readout")
    log("")
    log(f"MCC=2-significant cells that STAY significant at MCC=1 : {len(held)} / {len(sig2)}")
    log(f"MCC=2-significant cells that DROP at MCC=1             : {len(dropped)}")
    log(f"cells significant ONLY at MCC=1 (appeared)             : {len(appeared)}")

    log("")
    log("-" * 96)
    log("THE INFORMATIVE DIRECTION — cells that DROP despite MCC=1 giving them MORE cases")
    log("(holding is the expected outcome once power rises, so it is the drops that carry")
    log(" information: they mark associations specific to the stricter phenotype.)")
    log("-" * 96)
    if not dropped:
        log("  none — every MCC=2-significant cell retained significance.")
    for r in sorted(dropped, key=lambda r: -(r["case_ratio"] if np.isfinite(r["case_ratio"]) else 0)):
        log(f"  {r['sv_locus'][:44]:<44} {r['phecode']:<12} cases {r['mcc2_cases']:>7} -> "
            f"{r['mcc1_cases']:>7} ({r['case_ratio']:.2f}x)  q {r['mcc2_q']:.3g} -> "
            f"{r['mcc1_q'] if np.isfinite(r['mcc1_q']) else float('nan'):.3g}  "
            f"[MCC1 {r['mcc1_status']}]")

    sh = np.array([r["logp_shift"] for r in held if np.isfinite(r["logp_shift"])])
    if sh.size:
        log("")
        log(f"held cells: p STRENGTHENED in {int((sh > 0).sum())}/{sh.size}, weakened in "
            f"{int((sh < 0).sum())}/{sh.size} (median -log10p gain {np.median(sh):+.2f})")
    bs = [r["beta_sign_same"] for r in held if r["beta_sign_same"] != ""]
    if bs:
        log(f"held cells: effect-direction agreement {sum(bs)}/{len(bs)} "
            f"({100.0 * sum(bs) / len(bs):.1f}%)  [a sign flip would invalidate 'held']")

    log("")
    log("-" * 96)
    log("THE ASCERTAINMENT WATCH — cells significant ONLY at MCC=1")
    log("-" * 96)
    cat = {}
    for r in appeared:
        sv, ph = r["sv_locus"], r["phecode"]
        mcc2_absent = (r["mcc2_status"] == "not_tested")
        if sv.startswith(C16) and ph.startswith(SUBST):
            c = "chr16_substance_expansion"
        elif sv.startswith(C16):
            c = "chr16_other_phenotype"
        else:
            c = "new_axis"
        r["appear_class"] = c
        r["mcc2_phecode_absent"] = int(mcc2_absent)
        r["was_indicative_at_mcc2"] = int(np.isfinite(r["mcc2_q"]) and 0.05 <= r["mcc2_q"] < 0.10)
        cat.setdefault(c, []).append(r)
    for c in sorted(cat, key=lambda c: -len(cat[c])):
        rs = cat[c]
        nabs = sum(r["mcc2_phecode_absent"] for r in rs)
        nind = sum(r["was_indicative_at_mcc2"] for r in rs)
        log(f"  {c:<28} {len(rs):>4} cells  ({len({r['sv_locus'] for r in rs})} SVs, "
            f"{len({r['phecode'] for r in rs})} phecodes)  "
            f"[{nabs} on MCC=1-only phecodes (no MCC=2 counterpart), {nind} already indicative at MCC=2]")
    if appeared:
        with open(f"{OUT}/mcc1_appeared_coherence.tsv", "w", newline="") as f:
            ac = ["sv_locus", "tier", "phecode", "label", "appear_class", "mcc2_phecode_absent",
                  "was_indicative_at_mcc2", "mcc2_status", "mcc2_q", "mcc1_q", "mcc1_p", "mcc1_beta",
                  "mcc2_cases", "mcc1_cases", "case_ratio"]
            w = csv.DictWriter(f, fieldnames=ac, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in sorted(appeared, key=lambda r: (r["appear_class"], r["mcc1_p"])):
                w.writerow({k: (f"{r[k]:.4g}" if isinstance(r[k], float) and np.isfinite(r[k])
                                else ("" if isinstance(r[k], float) else r[k])) for k in ac})

    na = sorted([r for r in appeared if r["appear_class"] == "new_axis"], key=lambda r: r["mcc1_p"])
    if na:
        log("")
        log("  top 20 NEW-AXIS appeared cells (highest ascertainment suspicion):")
        for r in na[:20]:
            flag = " [MCC=1-only phecode]" if r["mcc2_phecode_absent"] else ""
            log(f"    {r['sv_locus'][:40]:<40} {r['phecode']:<12} {r['label'][:34]:<34} "
                f"cases {r['mcc1_cases']:>7}  q={r['mcc1_q']:.3g}{flag}")

    svs = sorted({k[0] for k in sig2 | sig1})
    with open(f"{OUT}/mcc1_persv.tsv", "w") as f:
        f.write("sv_locus\ttier\tn_sig_mcc2\tn_sig_mcc1\tn_held\tn_dropped\tn_appeared\n")
        for sv in svs:
            f.write(f"{sv}\t{a2['tiers'][a2['svidx'][sv]]}\t"
                    f"{len([k for k in sig2 if k[0] == sv])}\t{len([k for k in sig1 if k[0] == sv])}\t"
                    f"{len([r for r in held if r['sv_locus'] == sv])}\t"
                    f"{len([r for r in dropped if r['sv_locus'] == sv])}\t"
                    f"{len([r for r in appeared if r['sv_locus'] == sv])}\n")

    summ = dict(
        n_sv=g["n_sv"], covar_md5=g["covar_md5"],
        n_phecodes_mcc2=g["n_ph_mcc2"], n_phecodes_mcc1=g["n_ph_mcc1"],
        n_phecodes_mcc1_only=g["n_ph_mcc1_only"],
        case_inflation_median=g["infl_median"],
        case_inflation_min=g["infl_min"], case_inflation_max=g["infl_max"],
        sig_mcc2=len(sig2), sig_mcc1=len(sig1), ind_mcc2=len(ind2), ind_mcc1=len(ind1),
        svs_mcc2=len({k[0] for k in sig2}), svs_mcc1=len({k[0] for k in sig1}),
        held=len(held), dropped=len(dropped), appeared=len(appeared),
        chr16_subst_mcc2=len([k for k in sig2 if k[0].startswith(C16) and k[1].startswith(SUBST)]),
        chr16_subst_mcc1=len([k for k in sig1 if k[0].startswith(C16) and k[1].startswith(SUBST)]),
        appeared_by_class={c: len(v) for c, v in cat.items()},
        tier_sig_mcc2=dict(t2), tier_sig_mcc1=dict(t1),
        null_ext_sig_mcc2=t2.get("null_ext", 0), null_ext_sig_mcc1=t1.get("null_ext", 0),
        appeared_on_mcc1_only_phecodes=sum(r["mcc2_phecode_absent"] for r in appeared),
        held_p_strengthened=int((sh > 0).sum()) if sh.size else 0,
        held_p_total=int(sh.size),
        held_sign_agree=int(sum(bs)) if bs else 0, held_sign_total=len(bs),
    )
    with open(f"{OUT}/mcc1_vs_mcc2_summary.tsv", "w") as f:
        f.write("metric\tvalue\n")
        for k, v in summ.items():
            f.write(f"{k}\t{json.dumps(v) if isinstance(v, dict) else v}\n")

    with open(f"{OUT}/mcc1_vs_mcc2_gates.txt", "w") as f:
        f.write("\n".join(_log_lines) + "\n")
    log("")
    log(f"wrote -> {OUT}/{{mcc1_vs_mcc2_cells.tsv, mcc1_vs_mcc2_summary.tsv, "
        f"mcc1_appeared_coherence.tsv, mcc1_persv.tsv, mcc1_vs_mcc2_gates.txt}}")
    with open(f"{OUT}/mcc1_vs_mcc2_gates.txt", "w") as f:
        f.write("\n".join(_log_lines) + "\n")


if __name__ == "__main__":
    main()
