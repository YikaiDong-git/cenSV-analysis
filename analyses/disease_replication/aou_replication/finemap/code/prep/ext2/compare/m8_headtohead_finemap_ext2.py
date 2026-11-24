#!/usr/bin/env python3
"""M8-ext2 — the fine-mapping head-to-head, executed to the pre-agreed protocol."""
import os
import re
import csv
import sys
import json
import argparse
import statistics
import collections

FM = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
      "disease_replication/aou_replication/finemap")
PUB = f"{FM}/results"
EXT = f"{FM}/results_ext2"
PUB_TGT = f"{FM}/targets/target_manifest.tsv"
EXT_TGT = f"{FM}/targets_ext2/target_manifest.tsv"
OUT = f"{EXT}/headtohead"
TAB = f"{OUT}/tables"
VER = f"{OUT}/verdicts"
DOC = f"{OUT}/docs"

BASE_RE = re.compile(r"(chr[0-9XY]+_sv[0-9]+)")

EXPECT = {
    "n_fit_modes": 338, "n_fits": 169,
    "tag_SNP_only": 273, "SV_in_CS": 55, "unresolved": 10,
    "n_base_svs": 290, "n_raw_svs": 304, "competitors": 6446,
    "median_maxr": 0.9111190366, "n_below_085": 65,
    "median_maxr_nsv1": 0.784, "median_maxr_nsvgt1": 0.982,
    "n_nsv1": 96, "n_nsvgt1": 73,
    "chain": [55, 30, 7, 4, 4, 3, 2],
    "purity_any_fits": 10, "purity_any_svs": 6,
    "purity_2of3_fits": 9, "purity_2of3_svs": 5,
    "gap_rank5": 0.5599, "gap_rank6": 0.3580, "gap_chr16_in_top5": 5,
}

PUB_FLAGSHIP = "chr5.SO_394.21.EUR"
EXT_SIBLING = "chr5.SO_394.2.EUR"
CHR8_FIT = "chr8.SO_397.all-anc"
CHR8_NEW = "chr8.MS_712.AMR"
SV396 = "chr16_sv396"

LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def load(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def base(sv):
    m = BASE_RE.match(sv)
    return m.group(1) if m else sv


def r_fit_modes(M):
    return len(M)


def r_fits(M):
    return len({r["fit_id"] for r in M})


def r_verdicts(M):
    return collections.Counter(r["verdict"] for r in M)


def r_candidate_svs(manifest):
    raw = set()
    for r in manifest:
        raw.update(s for s in r["sv_ids"].split(",") if s)
    return len(raw), len({base(s) for s in raw})


def r_competitors(M):
    per = collections.defaultdict(set)
    for r in M:
        if r["mode"] == "all_window":
            per[r["chrom"]].add(int(num(r["n_sv"])))
    bad = {c: v for c, v in per.items() if len(v) > 1}
    return sum(max(v) for v in per.values()), len(per), bad


def r_distinctness(M):
    H = [r for r in M if r["mode"] == "hit_only" and r["max_abs_sv_snp_r"] not in ("", "NA")]
    out = {}
    for lab, sel in (("all", lambda r: True),
                     ("n_sv==1", lambda r: int(num(r["n_sv"])) == 1),
                     ("n_sv>1", lambda r: int(num(r["n_sv"])) > 1)):
        g = [num(r["max_abs_sv_snp_r"]) for r in H if sel(r)]
        out[lab] = dict(n=len(g),
                        median=statistics.median(g) if g else float("nan"),
                        below085=sum(1 for v in g if v < 0.85))
    return out


def r_filter_chain(M):
    steps = []
    s = [r for r in M if r["verdict"] == "SV_in_CS"]
    steps.append(("S0 verdict==SV_in_CS", len(s)))
    s = [r for r in s if r["mode"] == "hit_only"]
    steps.append(("S1 mode==hit_only", len(s)))
    s = [r for r in s if num(r["max_abs_sv_snp_r"], 1.0) < 0.85]
    steps.append(("S2 max|r| < 0.85", len(s)))
    s = [r for r in s if (r.get("best_sv_cs_variant") or "") != ""]
    steps.append(("S3 SV leads a CS", len(s)))
    s = [r for r in s if all(r.get(f"converged_{e}") == "True"
                             for e in ("indiv", "rss_lpm", "rss_logit"))]
    steps.append(("S4 3/3 engines converged (NON-BINDING)", len(s)))
    s = [r for r in s if num(r["ncase"]) >= 500]
    steps.append(("S5 ncase >= 500", len(s)))
    surv = {(base(r["hit_sv_ids"].split("|")[0]), r["stratum"]) for r in s}
    steps.append(("S6 collapse (base_sv, stratum)", len(surv)))
    return steps, sorted(surv)


def r_purity(M, root):
    any_f, any_sv, two_f, two_sv = set(), set(), set(), set()
    for r in M:
        if r["mode"] != "all_window":
            continue
        p = f"{root}/fits/{r['stratum']}/{r['fit_id']}/all_window/susie_credible_sets.tsv"
        if not os.path.exists(p):
            continue
        eng = collections.defaultdict(set)
        for cs in load(p):
            if cs["top_type"] == "sv":
                eng[base(cs["top_variant"])].add(cs["engine"])
        for sv, es in eng.items():
            any_f.add(r["fit_id"])
            any_sv.add(sv)
            if len(es) >= 2:
                two_f.add(r["fit_id"])
                two_sv.add(sv)
    return dict(any_fits=len(any_f), any_svs=len(any_sv), any_sv_list=sorted(any_sv),
                two_fits=len(two_f), two_svs=len(two_sv), two_sv_list=sorted(two_sv))


def r_matched(PM, EM):
    pk = {(r["fit_id"], r["mode"]): r for r in PM}
    ek = {(r["fit_id"], r["mode"]): r for r in EM}
    keys = []
    for k in set(pk) & set(ek):
        a, b = pk[k], ek[k]
        if (a["n_sv"] == b["n_sv"] and a["n_variants"] == b["n_variants"]
                and a["mcc"] == b["mcc"] and a["ncase"] == b["ncase"]
                and a["stratum"] == b["stratum"]):
            keys.append(k)
    return sorted(keys), pk, ek


def r_gap(M):
    by = collections.defaultdict(dict)
    for r in M:
        by[r["fit_id"]][r["mode"]] = r
    gaps = []
    for fid, d in by.items():
        if "all_window" in d and "hit_only" in d:
            a, h = num(d["all_window"]["best_top_sv_pip"], None), num(d["hit_only"]["best_top_sv_pip"], None)
            if a is None or h is None:
                continue
            gaps.append((a - h, fid))
    gaps.sort(reverse=True)
    return gaps


def gate_published():
    say("=" * 100)
    say("GATE — reproduce the published numbers with the recipes this comparison will use")
    say("=" * 100)
    M = load(f"{PUB}/finemap_master.tsv")
    T = load(PUB_TGT)
    rows = []
    ok = True

    def chk(name, got, want, tol=0.0):
        nonlocal ok
        good = (abs(got - want) <= tol) if isinstance(want, float) else (got == want)
        ok &= good
        gs = f"{got:.10f}" if isinstance(got, float) else str(got)
        ws = f"{want:.10f}" if isinstance(want, float) else str(want)
        say(f"  {name:<40} published {ws:>14}   recomputed {gs:>14}   "
            f"{'PASS' if good else '*** FAIL ***'}")
        rows.append(dict(statistic=name, published=ws, recomputed=gs, verdict="PASS" if good else "FAIL"))

    chk("fit-modes", r_fit_modes(M), EXPECT["n_fit_modes"])
    chk("fits", r_fits(M), EXPECT["n_fits"])
    V = r_verdicts(M)
    for k in ("tag_SNP_only", "SV_in_CS", "unresolved"):
        chk(f"verdict {k}", V[k], EXPECT[k])
    raw, bs = r_candidate_svs(T)
    chk("candidate SVs (raw sv_locus)", raw, EXPECT["n_raw_svs"])
    chk("candidate SVs (base-collapsed)", bs, EXPECT["n_base_svs"])
    comp, nchrom, bad = r_competitors(M)
    chk("competitor pool (per-chrom sum)", comp, EXPECT["competitors"])
    if bad:
        say(f"  ⚠️ chroms with >1 distinct all_window n_sv: {bad}")
    D = r_distinctness(M)
    chk("median max|r| (hit_only, all)", D["all"]["median"], EXPECT["median_maxr"], 1e-9)
    chk("count max|r| < 0.85", D["all"]["below085"], EXPECT["n_below_085"])
    chk("n fits with n_sv==1", D["n_sv==1"]["n"], EXPECT["n_nsv1"])
    chk("n fits with n_sv>1", D["n_sv>1"]["n"], EXPECT["n_nsvgt1"])
    chk("median max|r| n_sv==1", D["n_sv==1"]["median"], EXPECT["median_maxr_nsv1"], 5e-3)
    chk("median max|r| n_sv>1", D["n_sv>1"]["median"], EXPECT["median_maxr_nsvgt1"], 5e-3)
    steps, surv = r_filter_chain(M)
    got = [n for _, n in steps]
    chk("filter chain 55->2", got, EXPECT["chain"])
    say(f"  survivors: {surv}")
    P = r_purity(M, PUB)
    chk("purity any-engine (fits)", P["any_fits"], EXPECT["purity_any_fits"])
    chk("purity any-engine (SVs)", P["any_svs"], EXPECT["purity_any_svs"])
    chk("purity >=2-of-3 (fits)", P["two_fits"], EXPECT["purity_2of3_fits"])
    chk("purity >=2-of-3 (SVs)", P["two_svs"], EXPECT["purity_2of3_svs"])
    say(f"  dropped at >=2-of-3: {sorted(set(P['any_sv_list']) - set(P['two_sv_list']))}")
    G = r_gap(M)
    chk("gap rank5", G[4][0], EXPECT["gap_rank5"], 5e-4)
    chk("gap rank6", G[5][0], EXPECT["gap_rank6"], 5e-4)
    chk("chr16 in gap top5", sum(1 for _, f in G[:5] if f.startswith("chr16.")),
        EXPECT["gap_chr16_in_top5"])

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(VER, exist_ok=True)
    with open(f"{VER}/gate_published.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=["statistic", "published", "recomputed", "verdict"],
                           delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    say("")
    say(f"GATE {'PASS — every recipe reproduces its published value' if ok else 'FAILED'}")
    return ok, M, T


def wtsv(name, rows, fields):
    d = VER if name in ("attribution_ledger.tsv",) else TAB
    os.makedirs(d, exist_ok=True)
    with open(f"{d}/{name}", "w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def pip_of(root, M, fit_id, mode="hit_only"):
    row = next((r for r in M if r["fit_id"] == fit_id and r["mode"] == mode), None)
    if row is None:
        return None
    p = f"{root}/fits/{row['stratum']}/{fit_id}/{mode}/susie_fit_summary.tsv"
    s = load(p)[0] if os.path.exists(p) else {}
    nsv = int(num(row["n_sv"]))
    return dict(fit_id=fit_id, mode=mode, stratum=row["stratum"], mcc=row.get("mcc", ""),
                n_sv=nsv, nsv1_ok=(nsv == 1), ncase=int(num(row["ncase"])),
                top_sv_pip_indiv=num(s.get("top_sv_pip_indiv"), float("nan")),
                sum_sv_pip_indiv=num(row.get("sum_sv_pip_indiv"), float("nan")),
                max_abs_sv_snp_r=num(row.get("max_abs_sv_snp_r"), float("nan")),
                verdict=row["verdict"])


def ext2_side(PM, PT):
    EM, ET = load(f"{EXT}/finemap_master.tsv"), load(EXT_TGT)
    shared = {r["fit_id"] for r in PM} & {r["fit_id"] for r in EM}
    AUTO = [r for r in EM if r["chrom"] != "chrX"]
    SH_E = [r for r in EM if r["fit_id"] in shared]
    SH_P = [r for r in PM if r["fit_id"] in shared]
    led = []

    def L(stat, pub, ext, cls, why):
        led.append(dict(statistic=stat, published=pub, ext2=ext, attribution=cls, reason=why))

    say("\n" + "=" * 100)
    say("STEP 2 — DENOMINATORS (reported before any statistic)")
    say("=" * 100)
    den = []
    for lab, arm in (("published", PM), ("ext2", EM)):
        st = collections.Counter(r["stratum"] for r in arm)
        den.append(dict(arm=lab, fit_modes=len(arm), fits=len({r["fit_id"] for r in arm}),
                        chroms=len({r["chrom"] for r in arm}),
                        chrX_fits=len({r["fit_id"] for r in arm if r["chrom"] == "chrX"}),
                        AFR_fits=len({r["fit_id"] for r in arm if r["stratum"] == "AFR"}),
                        strata=";".join(f"{k}={v}" for k, v in sorted(st.items()))))
    for d in den:
        say(f"  {d['arm']:<10} fit-modes={d['fit_modes']:<5} fits={d['fits']:<5} chroms={d['chroms']:<4} "
            f"chrX_fits={d['chrX_fits']:<4} AFR_fits={d['AFR_fits']:<4} [{d['strata']}]")
    say(f"  shared fit_ids (like-for-like universe) = {len(shared)}")
    say("    NOTE the protocol's 128 is the MANIFEST overlap; 127 is the MASTER overlap. The one")
    say("    difference, chr7.MS_728.2.AMR, is _SKIP in BOTH arms (below case floor) — identical")
    say("    treatment on both sides, so it is a concordance, not a loss.")
    wtsv("denominators.tsv", den, ["arm", "fit_modes", "fits", "chroms", "chrX_fits", "AFR_fits", "strata"])
    L("fit-modes", len(PM), len(EM), "scope", "479 vs 170 manifest targets; +chrX +AFR")
    L("chrX fits", 0, len({r["fit_id"] for r in EM if r["chrom"] == "chrX"}), "scope",
      "old predictor store held ZERO chrX SVs; ext2 adds them")
    L("AFR fits", 0, len({r["fit_id"] for r in EM if r["stratum"] == "AFR"}), "scope",
      "published builder hardcoded STRATA without AFR")

    say("\n" + "=" * 100)
    say("STEP 3 — VERDICT PARTITION (three ways; only (c) is like-for-like)")
    say("=" * 100)
    vr = []
    for lab, arm in (("(a) ext2 all", EM), ("(b) ext2 autosomes only", AUTO),
                     ("(c) ext2 shared fit_ids", SH_E), ("(c') published shared fit_ids", SH_P),
                     ("published all", PM)):
        V = r_verdicts(arm)
        tot = sum(V.values())
        assert tot == len(arm), f"partition does not sum for {lab}"
        vr.append(dict(subset=lab, n=tot, tag_SNP_only=V["tag_SNP_only"],
                       SV_in_CS=V["SV_in_CS"], unresolved=V["unresolved"],
                       pct_SV_in_CS=round(100 * V["SV_in_CS"] / tot, 2) if tot else 0))
        say(f"  {lab:<30} n={tot:<5} tag_SNP_only={V['tag_SNP_only']:<5} SV_in_CS={V['SV_in_CS']:<4} "
            f"unresolved={V['unresolved']:<4} ({100 * V['SV_in_CS'] / tot:.1f}% SV_in_CS)")
    wtsv("verdicts.tsv", vr, ["subset", "n", "tag_SNP_only", "SV_in_CS", "unresolved", "pct_SV_in_CS"])
    sp, se = r_verdicts(SH_P), r_verdicts(SH_E)
    flips = sum(1 for f in shared for m in ("hit_only", "all_window")
                if (next((r["verdict"] for r in PM if r["fit_id"] == f and r["mode"] == m), None)
                    != next((r["verdict"] for r in EM if r["fit_id"] == f and r["mode"] == m), None)))
    say(f"  VERDICT FLIPS on shared fit-modes: {flips}")

    MK, pk, ek = r_matched(PM, EM)
    say("")
    say(f"  (d) MATCHED-PANEL subset = {len(MK)} fit-modes "
        f"({dict(collections.Counter(k[1] for k in MK))})")
    say("      same fit_id + stratum + mcc + ncase + n_sv + n_variants. 'Shared fit_id' is NOT enough:")
    say(f"      {len(shared) * 2} fit-modes share a fit_id but only {len(MK)} share a PANEL, so a delta")
    say("      read off the shared set silently absorbs B3 panel growth.")
    mflip = [k for k in MK if pk[k]["verdict"] != ek[k]["verdict"]]
    mv_p = collections.Counter(pk[k]["verdict"] for k in MK)
    mv_e = collections.Counter(ek[k]["verdict"] for k in MK)
    say(f"      published {dict(mv_p)}")
    say(f"      ext2      {dict(mv_e)}")
    say(f"      VERDICT FLIPS within matched panels: {len(mflip)} / {len(MK)}")
    for k in sorted(mflip):
        say(f"        {k[0]:<26} [{k[1]}] {pk[k]['verdict']} -> {ek[k]['verdict']}  "
            f"PIP {num(pk[k]['best_top_sv_pip']):.4f} -> {num(ek[k]['best_top_sv_pip']):.4f}  "
            f"max|r| {num(pk[k]['max_abs_sv_snp_r']):.4f} -> {num(ek[k]['max_abs_sv_snp_r']):.4f}")
    dl = sorted(((abs(num(ek[k]["best_top_sv_pip"]) - num(pk[k]["best_top_sv_pip"])), k) for k in MK),
                reverse=True)
    med = statistics.median([d for d, _ in dl]) if dl else float("nan")
    say(f"      median |ΔPIP| on matched panels = {med:.4f}   (identical geometry ⟹ attributable to"
        " the corrected DOSAGE alone)")
    say("      largest matched-panel PIP moves:")
    for d, k in dl[:6]:
        say(f"        {k[0]:<26} {num(pk[k]['best_top_sv_pip']):.4f} -> "
            f"{num(ek[k]['best_top_sv_pip']):.4f}  |Δ|={d:.4f}")
    wtsv("matched_panel.tsv",
         [dict(fit_id=k[0], mode=k[1], stratum=pk[k]["stratum"], mcc=pk[k]["mcc"],
               n_sv=pk[k]["n_sv"], n_variants=pk[k]["n_variants"], ncase=pk[k]["ncase"],
               pub_verdict=pk[k]["verdict"], ext_verdict=ek[k]["verdict"],
               pub_pip=round(num(pk[k]["best_top_sv_pip"]), 6),
               ext_pip=round(num(ek[k]["best_top_sv_pip"]), 6),
               pub_maxr=round(num(pk[k]["max_abs_sv_snp_r"]), 6),
               ext_maxr=round(num(ek[k]["max_abs_sv_snp_r"]), 6),
               flipped=pk[k]["verdict"] != ek[k]["verdict"]) for k in MK],
         ["fit_id", "mode", "stratum", "mcc", "n_sv", "n_variants", "ncase", "pub_verdict",
          "ext_verdict", "pub_pip", "ext_pip", "pub_maxr", "ext_maxr", "flipped"])
    L("SV_in_CS (shared fit_id)", sp["SV_in_CS"], se["SV_in_CS"], "composition",
      f"shared fit_id is NOT a matched panel — only {len(MK)} of {len(shared) * 2} share geometry")
    L("verdict flips (MATCHED panels)", f"0 of {len(MK)}", f"{len(mflip)} of {len(MK)}", "genuine",
      "identical geometry/mcc/ncase; the corrected dosage is the only difference")

    say("\n" + "=" * 100)
    say("STEP 4 — CANDIDATE / COMPETITOR SCOPE (units labelled — §3.2, §3.3)")
    say("=" * 100)
    praw, pbase = r_candidate_svs(PT)
    eraw, ebase = r_candidate_svs(ET)
    pc, pn, _ = r_competitors(PM)
    ec, en, ebad = r_competitors(EM)
    say(f"  candidate SVs   RAW  sv_locus : published {praw:<6} ext2 {eraw}")
    say(f"  candidate SVs   BASE collapsed: published {pbase:<6} ext2 {ebase}")
    say(f"  phecodes                      : published {len({r['phecode'] for r in PT}):<6} "
        f"ext2 {len({r['phecode'] for r in ET})}   <-- §3.2: ext2's 290 is PHECODES, "
        f"published's 290 is BASE SVs. Never 290 vs 290.")
    say(f"  competitor pool (per-chrom sum): published {pc:<6} ({pn} chroms)  ext2 {ec} ({en} chroms)")
    if ebad:
        say(f"  *** chroms with >1 distinct all_window n_sv: {ebad} — per-chrom sum ill-defined")
    wtsv("scope.tsv", [dict(unit="candidate SVs (raw sv_locus)", published=praw, ext2=eraw),
                       dict(unit="candidate SVs (base-collapsed)", published=pbase, ext2=ebase),
                       dict(unit="phecodes", published=len({r["phecode"] for r in PT}),
                            ext2=len({r["phecode"] for r in ET})),
                       dict(unit="competitor pool (per-chrom sum)", published=pc, ext2=ec),
                       dict(unit="chroms", published=pn, ext2=en)],
         ["unit", "published", "ext2"])
    L("candidate SVs (base)", pbase, ebase, "scope", "different target set by construction")
    L("competitor pool", pc, ec, "scope", "+3 chroms incl. chrX; larger all_window panels")

    say("\n" + "=" * 100)
    say("STEP 5 — DISTINCTNESS max|r| (four ways — B3: overall alone is a MIXTURE statistic)")
    say("=" * 100)
    DP, DE, DS = r_distinctness(PM), r_distinctness(EM), r_distinctness(SH_E)
    DSP = r_distinctness(SH_P)
    dr = []
    for lab, dp, de in (("overall", DP["all"], DE["all"]),
                        ("n_sv==1", DP["n_sv==1"], DE["n_sv==1"]),
                        ("n_sv>1", DP["n_sv>1"], DE["n_sv>1"])):
        say(f"  {lab:<10} published n={dp['n']:<5} median={dp['median']:.4f} <0.85={dp['below085']:<5}"
            f"   ext2 n={de['n']:<5} median={de['median']:.4f} <0.85={de['below085']}")
        dr.append(dict(subset=lab, pub_n=dp["n"], pub_median=round(dp["median"], 6),
                       pub_below085=dp["below085"], ext_n=de["n"],
                       ext_median=round(de["median"], 6), ext_below085=de["below085"]))
    say(f"  {'shared':<10} published n={DSP['all']['n']:<5} median={DSP['all']['median']:.4f} "
        f"<0.85={DSP['all']['below085']:<5}   ext2 n={DS['all']['n']:<5} "
        f"median={DS['all']['median']:.4f} <0.85={DS['all']['below085']}")
    dr.append(dict(subset="shared fit_ids", pub_n=DSP["all"]["n"],
                   pub_median=round(DSP["all"]["median"], 6), pub_below085=DSP["all"]["below085"],
                   ext_n=DS["all"]["n"], ext_median=round(DS["all"]["median"], 6),
                   ext_below085=DS["all"]["below085"]))
    say(f"  n_sv==1 : n_sv>1 ratio  published {DP['n_sv==1']['n']}:{DP['n_sv>1']['n']}   "
        f"ext2 {DE['n_sv==1']['n']}:{DE['n_sv>1']['n']}")
    MH = [k for k in MK if k[1] == "hit_only"]
    mmp = statistics.median([num(pk[k]["max_abs_sv_snp_r"]) for k in MH]) if MH else float("nan")
    mme = statistics.median([num(ek[k]["max_abs_sv_snp_r"]) for k in MH]) if MH else float("nan")
    say(f"  {'matched':<10} published n={len(MH):<5} median={mmp:.4f}"
        f"{'':<14}   ext2 n={len(MH):<5} median={mme:.4f}")
    dr.append(dict(subset="MATCHED panels", pub_n=len(MH), pub_median=round(mmp, 6), pub_below085="",
                   ext_n=len(MH), ext_median=round(mme, 6), ext_below085=""))
    say(f"  n_sv==1 : n_sv>1 ratio  published {DP['n_sv==1']['n']}:{DP['n_sv>1']['n']}   "
        f"ext2 {DE['n_sv==1']['n']}:{DE['n_sv>1']['n']}")
    d_overall = DE["all"]["median"] - DP["all"]["median"]
    d_shared = DS["all"]["median"] - DSP["all"]["median"]
    d_matched = mme - mmp
    say(f"  => Δmedian overall={d_overall:+.4f}   shared fit_id={d_shared:+.4f}   "
        f"MATCHED panel={d_matched:+.4f}")
    comp = abs(d_matched) < 0.01 and abs(d_overall) > 0.02
    say("  => " + ("the overall drop is COMPOSITIONAL: on matched panels tagging is unchanged "
                   f"({d_matched:+.4f}). The overall move tracks the {DP['n_sv==1']['n']}:"
                   f"{DP['n_sv>1']['n']} -> {DE['n_sv==1']['n']}:{DE['n_sv>1']['n']} mix, not tagging."
                   if comp else "matched-panel tagging ALSO moved — not purely compositional."))
    wtsv("distinctness.tsv", dr, ["subset", "pub_n", "pub_median", "pub_below085",
                                  "ext_n", "ext_median", "ext_below085"])
    L("median max|r| (overall)", round(DP["all"]["median"], 4), round(DE["all"]["median"], 4),
      "composition" if comp else "genuine",
      f"B3 mixture statistic; matched-panel Δ={d_matched:+.4f}, shared-fit Δ={d_shared:+.4f}")

    say("\n" + "=" * 100)
    say("STEP 6 — FILTER CHAIN, step by step (S4 is NON-BINDING by construction, §2.1)")
    say("=" * 100)
    ps, psurv = r_filter_chain(PM)
    es, esurv = r_filter_chain(EM)
    cr = []
    for (n1, c1), (n2, c2) in zip(ps, es):
        say(f"  {n1:<40} published {c1:<6} ext2 {c2}")
        cr.append(dict(step=n1, published=c1, ext2=c2))
    say(f"  published survivors: {psurv}")
    say(f"  ext2 survivors     : {esurv}")
    say(f"  survivors retained : {sorted(set(psurv) & set(esurv))}")
    say(f"  survivors lost     : {sorted(set(psurv) - set(esurv))}")
    say(f"  survivors new      : {sorted(set(esurv) - set(psurv))}")
    wtsv("filter_chain.tsv", cr, ["step", "published", "ext2"])
    L("chain endpoint (survivors)", len(psurv), len(esurv), "mixed",
      "scope (bigger target set) + composition (S2 is max|r|, a panel maximum)")

    say("\n" + "=" * 100)
    say("STEP 7 — NAMED SIGNALS (B1/B2 caveats are part of the result, not footnotes)")
    say("=" * 100)
    ns = []

    def named(tag, fit, arm_M, root, note=""):
        d = pip_of(root, arm_M, fit)
        if d is None:
            say(f"  {tag:<26} {fit:<26} ABSENT from this arm's master   {note}")
            return None
        aw = next((r for r in arm_M if r["fit_id"] == fit and r["mode"] == "all_window"), {})
        d.update(tag=tag, all_window_n_sv=int(num(aw.get("n_sv"))), note=note)
        flag = "" if d["nsv1_ok"] else "  *** n_sv>1: sum_sv_pip_indiv is NO LONGER the published quantity"
        say(f"  {tag:<26} {fit:<26} n_sv={d['n_sv']:<3} top_sv_pip_indiv={d['top_sv_pip_indiv']:.5f} "
            f"max|r|={d['max_abs_sv_snp_r']:.4f} ncase={d['ncase']:<7} mcc={d['mcc']:<5} "
            f"aw_n_sv={d['all_window_n_sv']:<5}{flag}")
        ns.append(d)
        return d

    say("  -- chr8_sv5 --")
    named("published", CHR8_FIT, PM, PUB)
    named("ext2", CHR8_FIT, EM, EXT)
    named("ext2 NEW fit", CHR8_NEW, EM, EXT, note="ext2-only; report as NEW, not as a change")
    say("  -- chr5_sv324 --")
    named("published FLAGSHIP", PUB_FLAGSHIP, PM, PUB)
    say(f"  B1: {PUB_FLAGSHIP} is NOT an ext2 target — "
        f"{'confirmed absent' if PUB_FLAGSHIP not in {r['fit_id'] for r in ET} else '*** PRESENT ***'}"
        " from targets_ext2. Report as 'did not reach ext2 candidacy', NEVER as a changed PIP.")
    named("published sibling", EXT_SIBLING, PM, PUB)
    named("ext2 sibling", EXT_SIBLING, EM, EXT)
    say("  -- chr16_sv396 (B2) --")
    cand = sum(1 for r in ET if SV396 in r["sv_ids"])
    hit_files = sum(1 for r in EM if r["mode"] == "hit_only" and SV396 in r["hit_sv_ids"])
    say(f"  B2: {SV396} reaches candidacy in {cand} ext2 manifest fits (published walk-back preserved, "
        f"but for a DIFFERENT reason: lost discovery candidacy, not the within-fit hit filter).")
    say(f"      {SV396} in ext2 hit_only hit_sv_ids: {hit_files} fit-modes")
    for arm_lab, arm_M, root in (("published", PM, PUB), ("ext2", EM, EXT)):
        for fit in ("chr16.MB_280.1.AMR", "chr16.MB_280.1.all-anc"):
            p = f"{root}/fits/{fit.split('.')[-1]}/{fit}/all_window/variant_meta.tsv"
            if not os.path.exists(p):
                continue
            row = next((v for v in load(p) if v["variant_id"].startswith(SV396 + "_")), None)
            if row:
                say(f"      {arm_lab:<10} {fit:<24} {SV396} marginal z={row['z']} p={row['p']}")
    for arm_lab, arm_M in (("published", PM), ("ext2", EM)):
        h = [r for r in arm_M if r["mode"] == "hit_only" and r["stratum"] == "AMR" and r["chrom"] == "chr16"]
        if h:
            say(f"      {arm_lab:<10} chr16 AMR hit_only: n fits={len(h)} "
                f"n_sv range={min(int(num(r['n_sv'])) for r in h)}-{max(int(num(r['n_sv'])) for r in h)} "
                f"median max|r|={statistics.median(num(r['max_abs_sv_snp_r']) for r in h):.7f}")
    wtsv("named_signals.tsv", ns, ["tag", "fit_id", "mode", "stratum", "mcc", "n_sv", "nsv1_ok",
                                   "ncase", "top_sv_pip_indiv", "sum_sv_pip_indiv",
                                   "max_abs_sv_snp_r", "all_window_n_sv", "verdict", "note"])
    L(f"{PUB_FLAGSHIP} PIP", "0.92714", "N/A — not an ext2 target", "scope",
      "B1: left the target set; a PIP delta here would be meaningless")
    L(f"{SV396} candidacy", "0 hit_only of 169", f"{cand} manifest fits", "scope",
      "B2: same conclusion, different mechanism")

    say("\n" + "=" * 100)
    say("STEP 8 — GAP RANKING (Δ best_top_sv_pip); the rank5/rank6 MARGIN is the claim")
    say("=" * 100)
    gr = []
    for lab, G in (("published", r_gap(PM)), ("ext2", r_gap(EM))):
        top = G[:5]
        r5, r6 = G[4][0], G[5][0]
        say(f"  {lab:<10} top5: " + ", ".join(f"{f}={g:.4f}" for g, f in top))
        say(f"  {lab:<10} rank5={r5:.4f} rank6={r6:.4f} ({G[5][1]}) margin={r5 / r6:.2f}x  "
            f"chr16 in top5={sum(1 for _, f in top if f.startswith('chr16.'))}/5")
        gr.append(dict(arm=lab, rank5=round(r5, 4), rank6=round(r6, 4), rank6_fit=G[5][1],
                       margin=round(r5 / r6, 2),
                       chr16_in_top5=sum(1 for _, f in top if f.startswith("chr16.")),
                       top5=";".join(f for _, f in top)))
    wtsv("gap_ranking.tsv", gr, ["arm", "rank5", "rank6", "rank6_fit", "margin", "chr16_in_top5", "top5"])
    ptf = {r["fit_id"] for r in PT}
    EG = r_gap(EM)
    newcomers = [f for _, f in EG[:6] if f not in ptf]
    say(f"  of ext2's top6, {len(newcomers)} are NEW ext2-only targets (never eligible to rank in the"
        f" published arm): {newcomers}")
    L("gap margin rank5/rank6", gr[0]["margin"], gr[1]["margin"],
      "scope" if newcomers else "genuine",
      f"{len(newcomers)} of ext2's top6 are new targets that could not have ranked before")

    say("\n" + "=" * 100)
    say("STEP 9 — PURITY (SV leads a CS), any-engine AND >=2-of-3 (§2.3)")
    say("=" * 100)
    PP, PE = r_purity(PM, PUB), r_purity(EM, EXT)
    say(f"  published  any-engine {PP['any_fits']:>3} fits / {PP['any_svs']:>3} SVs   "
        f">=2-of-3 {PP['two_fits']:>3} fits / {PP['two_svs']:>3} SVs")
    say(f"  ext2       any-engine {PE['any_fits']:>3} fits / {PE['any_svs']:>3} SVs   "
        f">=2-of-3 {PE['two_fits']:>3} fits / {PE['two_svs']:>3} SVs")
    say(f"  published >=2-of-3 SVs: {PP['two_sv_list']}")
    say(f"  ext2      >=2-of-3 SVs: {PE['two_sv_list']}")
    say(f"  retained: {sorted(set(PP['two_sv_list']) & set(PE['two_sv_list']))}")
    say(f"  lost    : {sorted(set(PP['two_sv_list']) - set(PE['two_sv_list']))}")
    say("  NOTE §2.3: 2-of-3 engines form ZERO purity CSs on chr8 partly because ~242 mutually collinear"
        " SVs push every 95% set below min_abs_corr=0.5. A panel-size change flips this for reasons"
        " unrelated to the signal.")
    wtsv("purity.tsv", [dict(arm="published", **{k: v for k, v in PP.items() if not k.endswith("_list")}),
                        dict(arm="ext2", **{k: v for k, v in PE.items() if not k.endswith("_list")})],
         ["arm", "any_fits", "any_svs", "two_fits", "two_svs"])
    L("purity >=2-of-3 (SVs)", PP["two_svs"], PE["two_svs"], "mixed",
      "scope (more all_window fits) + composition (panel size drives min_abs_corr)")

    say("\n" + "=" * 100)
    say("STEP 10 — ATTRIBUTION LEDGER (only 'genuine' is reportable as a finding)")
    say("=" * 100)
    for r in led:
        say(f"  [{r['attribution']:<12}] {r['statistic']:<34} {str(r['published']):>14} -> "
            f"{str(r['ext2']):<22} {r['reason']}")
    wtsv("attribution_ledger.tsv", led, ["statistic", "published", "ext2", "attribution", "reason"])
    tally = collections.Counter(r["attribution"] for r in led)
    say(f"\n  ledger tally: {dict(tally)}")
    say("  NOT REPORTABLE AS CHANGED (prior agreement, §2.1/§2.2): engine concordance (filtered 0 of 338"
        " published rows), the curated '10 signals' intermediate, any all_window-only PIP.")
    say("  LIMITATION (carried into every findings doc): no external replication cohort is available.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate-only", action="store_true",
                    help="published-side reproduction only; works before the ext2 fits exist")
    a = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(DOC, exist_ok=True)
    os.makedirs(VER, exist_ok=True)
    ok, PM, PT = gate_published()
    if not ok:
        say("\nRefusing to proceed: a recipe that cannot reproduce the published value cannot be used to "
            "claim an ext2 difference.")
        os.makedirs(DOC, exist_ok=True)
        with open(f"{DOC}/SUMMARY.md", "w") as f:
            f.write("\n".join(LOG) + "\n")
        sys.exit(1)

    if a.gate_only:
        say("\n--gate-only: stopping before the ext2 side.")
        os.makedirs(DOC, exist_ok=True)
        with open(f"{DOC}/SUMMARY.md", "w") as f:
            f.write("\n".join(LOG) + "\n")
        say("M8_GATE_ONLY_DONE")
        return

    em = f"{EXT}/finemap_master.tsv"
    if not os.path.exists(em):
        say(f"\next2 master table not found at {em} — run code_ext2/collate_finemap.py first.")
        sys.exit(2)
    ext2_side(PM, PT)
    os.makedirs(DOC, exist_ok=True)
    with open(f"{DOC}/SUMMARY.md", "w") as f:
        f.write("\n".join(LOG) + "\n")
    say("M8_DONE")


if __name__ == "__main__":
    main()
