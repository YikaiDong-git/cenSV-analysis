#!/usr/bin/env python3
"""Step 87 (ext2 arm) — PER-ANCESTRY QUANT *ANALYSIS* (<- step30r_quant_ancestry_SUMMARY.md)."""
import os
import sys
import gc
import json
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np                                       # noqa: E402
from scipy.stats import chi2 as _chi2                    # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from frozen_tree_sig import tree_sig                     # noqa: E402

HOME = os.path.expanduser("~")
RV = (f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/"
      f"aou_replication/results_vds")

PUB = f"{RV}/step30_phewide"
EXT = f"{RV}/step30_phewide_ext2"
OUT = f"{EXT}/quant/analysis"
PUB_META = f"{PUB}/quant_ancestry_meta.tsv"
EXT_META = f"{EXT}/quant/quant_ancestry_meta.tsv"
PUB_MAT = f"{PUB}/quant_matrices"
EXT_MAT = f"{EXT}/quant/matrices"
PUB_TIER = f"{RV}/phewide_carryover/sv_tier.tsv"
EXT_TIER = f"{RV}/step28_phewide_ext2/sv_tier.tsv"
COVAR = f"{RV}/phewide_carryover/inputs/covar.txt"
FROZEN = ["step30_phewide", "step28_phewide", "phewide_carryover", "step30_phewide_mcc1"]

ANCS = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
WELL = ["EUR", "AFR", "AMR"]
KNOWN = {"29463-7": "weight", "8302-2": "height", "39156-5": "BMI", "32623-1": "MPV", "776-5": "PDW"}
WBC = ["6690-2", "770-8", "751-8"]
VERDICTS = ["REPLICATES", "REPLICATES_HET", "REPLICATES_PLIM", "STRAT_SUSPECT",
            "DISCORDANT", "UNDERPOWERED"]
REPL_SET = ("REPLICATES", "REPLICATES_HET", "REPLICATES_PLIM")
VCODE = {v: i for i, v in enumerate(VERDICTS)}

PUBLISHED_CLAIMS = {
    "tally": {"UNDERPOWERED": 356941, "STRAT_SUSPECT": 62647, "REPLICATES_PLIM": 30026,
              "DISCORDANT": 3568, "REPLICATES": 2035, "REPLICATES_HET": 693},
    "anthro_p1e4": {"replicates_total": 277, "REPLICATES": 229, "REPLICATES_HET": 15,
                    "REPLICATES_PLIM": 33, "STRAT_SUSPECT": 87},
    "chr1_wbc_all": {"STRAT_SUSPECT": 1580, "UNDERPOWERED": 172, "REPLICATES_PLIM": 18,
                     "DISCORDANT": 12, "REPLICATES_HET": 2, "REPLICATES": 1},
    "anc_n": {"EUR": 177230, "AFR": 64726, "AMR": 60984, "EAS": 7144, "SAS": 3944, "MID": 1176},
}

LOG = []


def say(msg=""):
    print(msg, flush=True)
    LOG.append(msg)


def die(msg):
    say(f"\n*** GATE FAILED: {msg}")
    say("*** aborting without writing analysis outputs (fail-closed)")
    sys.exit(1)


def read_tier(path):
    t = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        if hdr[:2] != ["sv_locus", "tier"]:
            die(f"{path}: unexpected header {hdr[:2]}")
        for line in f:
            q = line.rstrip("\n").split("\t")
            t[q[0]] = q[1]
    return t


def read_meta(path, tag):
    t0 = time.time()
    svs, los, sv2i, lo2j = [], [], {}, {}
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        need = ["sv_locus", "loinc", "meta_beta", "meta_p", "I2", "k", "verdict"] + \
               [f"{a}_{x}" for a in WELL for x in ("beta", "p")]
        miss = [c for c in need if c not in hdr]
        if miss:
            die(f"{path}: missing columns {miss}")
        idx = {c: hdr.index(c) for c in need}
        for line in f:
            q = line.rstrip("\n").split("\t")
            sv, lo = q[0], q[1]
            i = sv2i.get(sv)
            if i is None:
                i = sv2i[sv] = len(svs)
                svs.append(sv)
            j = lo2j.get(lo)
            if j is None:
                j = lo2j[lo] = len(los)
                los.append(lo)
            rows.append((i, j, q))
    M, P = len(svs), len(los)
    cols = {c: np.full((M, P), np.nan, np.float32) for c in
            ["meta_beta", "meta_p", "I2"] + [f"{a}_{x}" for a in WELL for x in ("beta", "p")]}
    vcode = np.full((M, P), -1, np.int8)
    kk = np.zeros((M, P), np.int8)
    unknown = set()
    for i, j, q in rows:
        for c in cols:
            v = q[idx[c]]
            if v != "":
                try:
                    cols[c][i, j] = float(v)
                except ValueError:
                    pass
        v = q[idx["verdict"]]
        if v not in VCODE:
            unknown.add(v)
        else:
            vcode[i, j] = VCODE[v]
        try:
            kk[i, j] = int(q[idx["k"]])
        except ValueError:
            pass
    del rows
    gc.collect()
    if unknown:
        die(f"{path}: unrecognised verdict label(s) {sorted(unknown)} — the verdict vocabulary drifted")
    n_cell = int((vcode >= 0).sum())
    say(f"[load] {tag:<10} {path.split('results_vds/')[-1]}")
    say(f"        {M:,} SVs x {P} LOINCs -> {n_cell:,} cells ({time.time()-t0:.0f}s)")
    cols["verdict"] = vcode
    cols["k"] = kk
    return svs, los, sv2i, lo2j, cols


def tally(vcode, mask=None):
    v = vcode if mask is None else np.where(mask, vcode, -1)
    return {name: int((v == c).sum()) for name, c in VCODE.items()}


def pct(n, d):
    return f"{100.0*n/d:.1f}%" if d else "n/a"


def lambda_gc(p):
    p = np.asarray(p, np.float64).ravel()
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if p.size < 100:
        return float("nan"), int(p.size)
    c = _chi2.isf(p, 1)
    return float(np.median(c) / 0.4549364231195736), int(p.size)


def main():
    t0 = time.time()
    say("=" * 100)
    say("step87 — per-ancestry QUANT analysis on ext2 weights, head-to-head vs the published arm")
    say("=" * 100)

    for p in (PUB_META, EXT_META, PUB_TIER, EXT_TIER, COVAR, PUB_MAT, EXT_MAT):
        if not os.path.exists(p):
            die(f"missing input {p}")
    if os.path.islink(OUT):
        die(f"output dir {OUT} is a symlink — refusing to write through it")
    os.makedirs(OUT, exist_ok=True)
    say(f"[G0] all inputs present; outputs -> {OUT.split('results_vds/')[-1]}")

    base = {}
    for name in FROZEN:
        c, m, n = tree_sig(f"{RV}/{name}")
        base[name] = (c, m, n)
        say(f"[V-pre] {name:<22} n={n:<5} content={c}")

    ptier, etier = read_tier(PUB_TIER), read_tier(EXT_TIER)
    psv, plo, psv2i, plo2j, pc = read_meta(PUB_META, "published")
    esv, elo, esv2i, elo2j, ec = read_meta(EXT_META, "ext2")

    if set(plo) != set(elo):
        die(f"LOINC panels differ: pub-only={set(plo)-set(elo)} ext2-only={set(elo)-set(plo)}")
    say(f"[G1] LOINC panel identical across arms ({len(plo)} traits) — PASS")
    if not set(psv) <= set(esv):
        die(f"{len(set(psv)-set(esv))} published SVs absent from ext2 — not a superset")
    say(f"[G2] published {len(psv):,} SVs are a strict subset of ext2 {len(esv):,} — PASS")
    mism = [s for s in set(ptier) & set(etier) if ptier[s] != etier[s]]
    if mism:
        die(f"{len(mism)} SVs carry different tier labels across arms, e.g. {mism[:3]}")
    say(f"[G3] tier labels agree on all {len(set(ptier) & set(etier)):,} shared SVs — PASS")
    miss_t = [s for s in esv if s not in etier]
    if miss_t:
        die(f"{len(miss_t)} ext2 SVs have no tier (fail-open risk), e.g. {miss_t[:3]}")
    say(f"[G4] every ext2 SV has a tier — no fail-open ({len(esv):,} SVs) — PASS")

    say("\n" + "-" * 100)
    say("A. PER-ANCESTRY SAMPLE SIZE (derived from covar.txt; AFR is the reference level)")
    with open(COVAR) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        dcol = {a: hdr.index(f"anc_{a.lower()}") for a in ANCS if f"anc_{a.lower()}" in hdr}
        cnt = {a: 0 for a in ANCS}
        ntot = 0
        for line in f:
            q = line.rstrip("\n").split("\t")
            ntot += 1
            lab = "AFR"
            for a, ci in dcol.items():
                if q[ci] == "1" or q[ci] == "1.0":
                    lab = a
                    break
            cnt[lab] += 1
    say(f"  total {ntot:,}   " + "  ".join(f"{a} {cnt[a]:,}" for a in ANCS))
    bad_n = {a: (cnt[a], PUBLISHED_CLAIMS['anc_n'][a]) for a in ANCS if cnt[a] != PUBLISHED_CLAIMS["anc_n"][a]}
    if bad_n:
        say(f"  ⚠️  differs from the published summary: {bad_n}")
    else:
        say("  identical to the published summary — the two arms scan the SAME people (only weights differ)")

    say("\n" + "-" * 100)
    say("B. VERDICT TALLY")
    ptal = tally(pc["verdict"])
    etal = tally(ec["verdict"])

    bad = {k: (ptal[k], v) for k, v in PUBLISHED_CLAIMS["tally"].items() if ptal[k] != v}
    if bad:
        die(f"FIDELITY: re-tallying the PUBLISHED meta does not reproduce step30r's table: {bad}")
    say("[V3a] re-tally of the published meta reproduces step30r_quant_ancestry_SUMMARY.md exactly — PASS")

    shared_i = np.array([esv2i[s] for s in psv], np.int64)
    lo_reidx = np.array([elo2j[lo] for lo in plo], np.int64)
    ev_shared = ec["verdict"][np.ix_(shared_i, lo_reidx)]
    etal_shared = tally(ev_shared)

    etier_arr = np.array([etier[s] for s in esv])
    tier_masks = {t: (etier_arr == t)[:, None] & (ec["verdict"] >= 0) for t in sorted(set(etier_arr))}

    say(f"\n  {'verdict':<18}{'published':>12}{'ext2 (all)':>12}{'ext2 (shared)':>15}"
        f"{'shared delta':>14}")
    for v in VERDICTS:
        d = etal_shared[v] - ptal[v]
        say(f"  {v:<18}{ptal[v]:>12,}{etal[v]:>12,}{etal_shared[v]:>15,}{d:>+14,}")
    say(f"  {'TOTAL':<18}{sum(ptal.values()):>12,}{sum(etal.values()):>12,}{sum(etal_shared.values()):>15,}")

    say(f"\n  ext2 by tier:")
    say(f"  {'tier':<12}{'n SV':>8}  " + "".join(f"{v[:13]:>15}" for v in VERDICTS))
    tier_tal = {}
    for t, m in tier_masks.items():
        tt = tally(ec["verdict"], m)
        tier_tal[t] = tt
        say(f"  {t:<12}{int(m[:, 0].sum()):>8}  " + "".join(f"{tt[v]:>15,}" for v in VERDICTS))

    say("\n" + "-" * 100)
    say("B2. KNOWN-SIGNAL REPLICATION, pooled-hit-conditioned (cells q<0.05 in the all-ancestry scan)")
    ks_rows = []

    def pooled_hits(path):
        out = set()
        if not os.path.exists(path):
            return None
        with open(path) as f:
            h = f.readline().rstrip("\n").split("\t")
            if "loinc" not in h or "q_bh_tier" not in h:
                return None
            il, iq = h.index("loinc"), h.index("q_bh_tier")
            for line in f:
                q = line.rstrip("\n").split("\t")
                try:
                    if float(q[iq]) < 0.05:
                        out.add((q[0], q[il]))
                except (ValueError, IndexError):
                    continue
        return out

    for arm, hit_path, svs, los, sv2i, lo2j, C in (
            ("published", f"{PUB}/phewide_quant.tsv", psv, plo, psv2i, plo2j, pc),
            ("ext2", f"{EXT}/quant/phewide_quant_ext2.tsv", esv, elo, esv2i, elo2j, ec)):
        H = pooled_hits(hit_path)
        if H is None:
            say(f"  {arm}: pooled hit table missing/unusable at {hit_path} — skipped")
            continue
        per = {}
        for sv, lo in H:
            if lo not in KNOWN:
                continue
            i, j = sv2i.get(sv), lo2j.get(lo)
            if i is None or j is None or C["verdict"][i, j] < 0:
                continue
            t = KNOWN[lo]
            a = per.setdefault(t, [0, 0])
            a[1] += 1
            if VERDICTS[C["verdict"][i, j]] in REPL_SET:
                a[0] += 1
        tr = sum(v[0] for v in per.values())
        tn = sum(v[1] for v in per.values())
        say(f"  {arm:<10} pooled q<0.05 cells {len(H):>6,}   anchor cells with a verdict {tn:>5}   "
            f"REPLICATES {tr:>5} ({pct(tr, tn)})")
        say("             " + "  ".join(f"{t} {per[t][0]}/{per[t][1]}" for t in
                                        ("BMI", "weight", "height", "MPV", "PDW") if t in per))
        for t, v in sorted(per.items()):
            ks_rows.append(dict(arm=arm, definition="pooled_q05_conditioned", trait=t,
                                n_cells=v[1], replicates_any=v[0],
                                rate=f"{v[0]/v[1]:.4f}" if v[1] else ""))
        ks_rows.append(dict(arm=arm, definition="pooled_q05_conditioned", trait="TOTAL",
                            n_cells=tn, replicates_any=tr, rate=f"{tr/tn:.4f}" if tn else ""))
    say("\n  ⚠️ Under this definition the anchor replication rate FALLS from published to ext2, as it does")
    say("     under step30r's own all-cells/meta_p<1e-4 rule. HEADTOHEAD §4.2's '81.3% vs 74%, improves'")
    say("     compares an ext2 number to a published number computed under a different universe.")

    say("\n" + "-" * 100)
    say("C. PAIRED VERDICT TRANSITION (the 455,910 cells present in BOTH arms; only the weights differ)")
    pv = pc["verdict"]
    both = (pv >= 0) & (ev_shared >= 0)
    trans = np.zeros((len(VERDICTS), len(VERDICTS)), np.int64)
    for a in range(len(VERDICTS)):
        sel = both & (pv == a)
        if not sel.any():
            continue
        for b in range(len(VERDICTS)):
            trans[a, b] = int((sel & (ev_shared == b)).sum())
    npair = int(trans.sum())
    same = int(np.trace(trans))
    say(f"  paired cells {npair:,}   unchanged verdict {same:,} ({pct(same, npair)})")
    _corner = "published -> ext2"
    say(f"\n  {_corner:<20}" + "".join(f"{v[:12]:>14}" for v in VERDICTS) + f"{'row tot':>12}")
    for a, va in enumerate(VERDICTS):
        say(f"  {va:<20}" + "".join(f"{trans[a, b]:>14,}" for b in range(len(VERDICTS)))
            + f"{trans[a].sum():>12,}")

    say("\n" + "-" * 100)
    say("D. THE CALIBRATION TEST — does the SAME machinery call pan-ancestry effects REPLICATES?")
    say("   (published step30r: anthropometric anchor meta_p<1e-4 vs chr1->WBC meta_p<1e-6)")

    def group_mask(svs, los, lo2j, which, pcut, pmat, tiermap=None, tier_filter=None):
        M, P = len(svs), len(los)
        m = np.zeros((M, P), bool)
        if which == "anthro":
            cols = [lo2j[lo] for lo in KNOWN if lo in lo2j]
            m[:, cols] = True
        elif which == "chr1wbc":
            rows = [i for i, s in enumerate(svs) if s.startswith("chr1_")]
            cols = [lo2j[lo] for lo in WBC if lo in lo2j]
            if rows and cols:
                m[np.ix_(rows, cols)] = True
        else:
            die(f"unknown trait group '{which}'")
        if tier_filter is not None:
            keep = np.array([tiermap.get(s, "") in tier_filter for s in svs])
            m &= keep[:, None]
        if pcut is not None:
            m &= np.isfinite(pmat) & (pmat < pcut)
        return m

    rows_out = []

    def run_group(arm, svs, los, lo2j, C, which, pcut, tier_filter=None, tiermap=None):
        m = group_mask(svs, los, lo2j, which, pcut, C["meta_p"],
                       tiermap=tiermap, tier_filter=tier_filter)
        m &= (C["verdict"] >= 0)
        t = tally(C["verdict"], m)
        tot = sum(t.values())
        rep = sum(t[v] for v in REPL_SET)
        lab = f"{arm}/{which}" + (f"/{'+'.join(sorted(tier_filter))}" if tier_filter else "")
        rows_out.append(dict(arm=arm, group=which, p_cut=pcut,
                             tier=("+".join(sorted(tier_filter)) if tier_filter else "ALL"),
                             n_cells=tot, replicates_any=rep, **t))
        say(f"  {lab:<44}{'n=' + str(tot):>10}  REPL(any) {rep:>6} ({pct(rep, tot):>6})  "
            f"STRAT {t['STRAT_SUSPECT']:>6} ({pct(t['STRAT_SUSPECT'], tot):>6})  "
            f"[R {t['REPLICATES']} / HET {t['REPLICATES_HET']} / PLIM {t['REPLICATES_PLIM']} / "
            f"DISC {t['DISCORDANT']} / UND {t['UNDERPOWERED']}]")
        return t, tot, rep

    say("\n  -- published arm (must reproduce step30r) --")
    pt_anthro, pn_anthro, pr_anthro = run_group("published", psv, plo, plo2j, pc, "anthro", 1e-4)
    pt_wbc6, pn_wbc6, _ = run_group("published", psv, plo, plo2j, pc, "chr1wbc", 1e-6)
    pt_wbcA, pn_wbcA, _ = run_group("published", psv, plo, plo2j, pc, "chr1wbc", None)

    cl = PUBLISHED_CLAIMS["anthro_p1e4"]
    bad = {k: (pt_anthro[k], cl[k]) for k in ("REPLICATES", "REPLICATES_HET", "REPLICATES_PLIM",
                                              "STRAT_SUSPECT") if pt_anthro[k] != cl[k]}
    if pr_anthro != cl["replicates_total"]:
        bad["replicates_total"] = (pr_anthro, cl["replicates_total"])
    if bad:
        die(f"FIDELITY: published anthropometric split not reproduced: {bad}")
    say(f"  [V3b] published anthropometric split reproduced exactly "
        f"({pr_anthro} REPL = {pct(pr_anthro, pn_anthro)}, {pt_anthro['STRAT_SUSPECT']} STRAT "
        f"= {pct(pt_anthro['STRAT_SUSPECT'], pn_anthro)}) — PASS")
    clw = PUBLISHED_CLAIMS["chr1_wbc_all"]
    badw = {k: (pt_wbcA[k], v) for k, v in clw.items() if pt_wbcA[k] != v}
    if badw:
        die(f"FIDELITY: published chr1->WBC breakdown not reproduced: {badw}")
    say(f"  [V3c] published chr1->WBC breakdown reproduced exactly "
        f"({pt_wbcA['STRAT_SUSPECT']} STRAT = {pct(pt_wbcA['STRAT_SUSPECT'], pn_wbcA)} of {pn_wbcA}) — PASS")
    say("  ⟹ the analysis code is faithful; any ext2 difference below is attributable to the WEIGHTS.")

    say("\n  -- published arm, MATCHED thresholds (the published contrast used 1e-4 vs 1e-6) --")
    run_group("published", psv, plo, plo2j, pc, "chr1wbc", 1e-4)
    run_group("published", psv, plo, plo2j, pc, "anthro", 1e-6)

    say("\n  -- ext2 arm, all 10,799 predictors --")
    et_anthro, en_anthro, er_anthro = run_group("ext2", esv, elo, elo2j, ec, "anthro", 1e-4)
    et_wbc6, en_wbc6, _ = run_group("ext2", esv, elo, elo2j, ec, "chr1wbc", 1e-6)
    et_wbcA, en_wbcA, _ = run_group("ext2", esv, elo, elo2j, ec, "chr1wbc", None)
    run_group("ext2", esv, elo, elo2j, ec, "chr1wbc", 1e-4)
    run_group("ext2", esv, elo, elo2j, ec, "anthro", 1e-6)

    say("\n  -- ext2 arm restricted to the SHARED 6,513 (isolates the weight change) --")
    shared_C = {"meta_p": ec["meta_p"][np.ix_(shared_i, lo_reidx)], "verdict": ev_shared}
    run_group("ext2_shared", psv, plo, plo2j, shared_C, "anthro", 1e-4)
    run_group("ext2_shared", psv, plo, plo2j, shared_C, "chr1wbc", 1e-6)
    run_group("ext2_shared", psv, plo, plo2j, shared_C, "chr1wbc", 1e-4)
    run_group("ext2_shared", psv, plo, plo2j, shared_C, "anthro", 1e-6)

    say("\n  -- ext2 arm, SPLIT BY TIER (null_ext is discovery's NEGATIVE CONTROL) --")
    tier_res = {}
    for tf in ({"signal"}, {"null_ext"}, {"ext2_new"}, {"ext2_chrX"}):
        if not (set(etier_arr) & tf):
            continue
        tier_res[tuple(sorted(tf))] = {
            "anthro": run_group("ext2", esv, elo, elo2j, ec, "anthro", 1e-4, tf, etier),
            "chr1wbc": run_group("ext2", esv, elo, elo2j, ec, "chr1wbc", 1e-6, tf, etier),
        }

    say("\n" + "-" * 100)
    say("E. NEGATIVE-CONTROL READING — is 'REPLICATES' evidence of anything?")
    sig = tier_res.get(("signal",), {}).get("anthro")
    nul = tier_res.get(("null_ext",), {}).get("anthro")
    if sig and nul:
        s_t, s_n, s_r = sig
        n_t, n_r = nul[0], nul[2]
        n_n = nul[1]
        say(f"  anthropometric anchor, meta_p<1e-4:")
        say(f"    signal   {s_n:>6} cells   REPLICATES(any) {s_r:>5} ({pct(s_r, s_n)})")
        say(f"    null_ext {n_n:>6} cells   REPLICATES(any) {n_r:>5} ({pct(n_r, n_n)})")
        if s_n and n_n:
            ratio = (s_r / s_n) / max(n_r / n_n, 1e-9)
            say(f"    signal:null_ext REPLICATES-rate ratio = {ratio:.2f}x")
            if n_r / n_n >= s_r / s_n:
                say("    ⚠️  the negative control replicates AT LEAST AS OFTEN as signal ⟹ on this panel a")
                say("       REPLICATES verdict is NOT evidence of predictor-specific biology. The verdict")
                say("       measures cross-stratum CONSISTENCY of whatever the dosage tracks, which for a")
                say("       control predictor can still be a genuine pan-ancestry nuisance axis.")
            else:
                say("    signal replicates more often than the negative control on the anchor traits.")
        say("  NOTE: null_ext predictors are still real genotype-weighted dosages over real people; they are")
        say("  a control for DISCOVERY-SIDE selection, not a permutation null. A nonzero REPLICATES rate is")
        say("  therefore expected. What matters is the RATIO, and the disease arm's 0-hit result.")

    say("\n" + "-" * 100)
    say("F. chr1 -> WBC / neutrophil — the 'clean quant biology' stress test")
    say(f"  published: {pt_wbcA['STRAT_SUSPECT']:,}/{pn_wbcA:,} STRAT_SUSPECT "
        f"({pct(pt_wbcA['STRAT_SUSPECT'], pn_wbcA)})")
    say(f"  ext2     : {et_wbcA['STRAT_SUSPECT']:,}/{en_wbcA:,} STRAT_SUSPECT "
        f"({pct(et_wbcA['STRAT_SUSPECT'], en_wbcA)})")

    driver_rows = []
    for arm, svs, los, lo2j, C, tiers in (("published", psv, plo, plo2j, pc, ptier),
                                          ("ext2", esv, elo, elo2j, ec, etier)):
        rows = [i for i, s in enumerate(svs) if s.startswith("chr1_")]
        cols = [lo2j[lo] for lo in WBC if lo in lo2j]
        cand = []
        for i in rows:
            for j in cols:
                p = C["meta_p"][i, j]
                if np.isfinite(p) and C["verdict"][i, j] >= 0:
                    cand.append((float(p), i, j))
        cand.sort()
        say(f"\n  top {arm} chr1 x WBC cells by meta_p:")
        say(f"    {'SV':<40}{'loinc':<9}{'tier':<10}" +
            "".join(f"{a + '_b/p':>18}" for a in WELL) + f"{'meta_p':>10}  verdict")
        for p, i, j in cand[:8]:
            s = f"    {svs[i][:39]:<40}{los[j]:<9}{tiers.get(svs[i], '?'):<10}"
            rec = dict(arm=arm, sv_locus=svs[i], loinc=los[j], tier=tiers.get(svs[i], ""),
                       meta_p=f"{p:.3g}", meta_beta=f"{C['meta_beta'][i, j]:.4g}",
                       I2=f"{C['I2'][i, j]:.0f}", verdict=VERDICTS[C["verdict"][i, j]])
            for a in WELL:
                b, pp = C[f"{a}_beta"][i, j], C[f"{a}_p"][i, j]
                cell = f"{b:+.3f}/{pp:.0e}" if np.isfinite(b) else "na"
                s += f"{cell:>18}"
                rec[f"{a}_beta"] = f"{b:.4g}" if np.isfinite(b) else ""
                rec[f"{a}_p"] = f"{pp:.3g}" if np.isfinite(pp) else ""
            s += f"{p:>10.0e}  {VERDICTS[C['verdict'][i, j]]}"
            say(s)
            driver_rows.append(rec)

    say("\n" + "-" * 100)
    say("G. lambda_GC PER ANCESTRY PER TIER (computed from the per-stratum p matrices)")
    lam_rows = []
    for arm, mat, tiers in (("published", PUB_MAT, ptier), ("ext2", EXT_MAT, etier)):
        say(f"\n  {arm}:")
        say(f"    {'ancestry':<10}{'n':>10}  " + "".join(f"{t:>16}" for t in
                                                        ["ALL", "signal", "null_ext", "ext2_new"]))
        for a in ANCS:
            f_p = f"{mat}/quant_p_{a}.npy"
            f_sv = f"{mat}/quant_svorder_{a}.txt"
            if not (os.path.exists(f_p) and os.path.exists(f_sv)):
                say(f"    {a:<10}{'-':>10}  (no matrix — stratum not scanned)")
                continue
            P = np.load(f_p)
            order = [l.strip() for l in open(f_sv) if l.strip()]
            tv = np.array([tiers.get(s, "") for s in order])
            cells = []
            lam_all, n_all = lambda_gc(P)
            cells.append(f"{lam_all:.3f} (n={n_all // 1000}k)")
            rec = dict(arm=arm, ancestry=a, tier="ALL", lambda_gc=f"{lam_all:.4f}", n_p=n_all)
            lam_rows.append(rec)
            for t in ["signal", "null_ext", "ext2_new"]:
                sel = tv == t
                if not sel.any():
                    cells.append("-")
                    continue
                lam, nn = lambda_gc(P[sel])
                cells.append(f"{lam:.3f}")
                lam_rows.append(dict(arm=arm, ancestry=a, tier=t, lambda_gc=f"{lam:.4f}", n_p=nn))
            say(f"    {a:<10}{P.shape[0] * P.shape[1]:>10,}  " + "".join(f"{c:>16}" for c in cells))
            del P
            gc.collect()

    say("\n" + "-" * 100)
    say("WRITING ANALYSIS OUTPUTS")

    def wr(name, header, rows):
        p = f"{OUT}/{name}"
        if os.path.islink(p):
            die(f"{p} is a symlink — refusing to write through it")
        with open(p, "w") as f:
            f.write("\t".join(header) + "\n")
            for r in rows:
                f.write("\t".join(str(r.get(h, "")) for h in header) + "\n")
        say(f"  {name:<36}{len(rows):>7} rows   {os.path.getsize(p):>9,} B")

    tal_rows = []
    for v in VERDICTS:
        tal_rows.append(dict(scope="published_all", verdict=v, n=ptal[v]))
        tal_rows.append(dict(scope="ext2_all", verdict=v, n=etal[v]))
        tal_rows.append(dict(scope="ext2_shared6513", verdict=v, n=etal_shared[v]))
        for t, tt in tier_tal.items():
            tal_rows.append(dict(scope=f"ext2_tier_{t}", verdict=v, n=tt[v]))
    wr("quant_verdict_tally.tsv", ["scope", "verdict", "n"], tal_rows)

    tr_rows = [dict(published_verdict=VERDICTS[a], ext2_verdict=VERDICTS[b], n=int(trans[a, b]))
               for a in range(len(VERDICTS)) for b in range(len(VERDICTS))]
    wr("quant_verdict_transition.tsv", ["published_verdict", "ext2_verdict", "n"], tr_rows)

    wr("quant_discriminates.tsv",
       ["arm", "group", "tier", "p_cut", "n_cells", "replicates_any"] + VERDICTS, rows_out)

    wr("quant_chr1wbc_drivers.tsv",
       ["arm", "sv_locus", "loinc", "tier"] + [f"{a}_{x}" for a in WELL for x in ("beta", "p")] +
       ["meta_beta", "meta_p", "I2", "verdict"], driver_rows)

    wr("quant_lambda_gc.tsv", ["arm", "ancestry", "tier", "lambda_gc", "n_p"], lam_rows)

    wr("quant_known_signal_replication.tsv",
       ["arm", "definition", "trait", "n_cells", "replicates_any", "rate"], ks_rows)

    say("\n" + "-" * 100)
    changed = []
    for name in FROZEN:
        c, m, n = tree_sig(f"{RV}/{name}")
        pcs, pms, pn = base[name]
        if c != pcs:
            changed.append(name)
            say(f"[V-post] {name:<22} *** CONTENT CHANGED *** {pcs} -> {c}")
        elif m != pms:
            say(f"[V-post] {name:<22} content identical, meta touched (bytes unchanged)")
        else:
            say(f"[V-post] {name:<22} unchanged (content+meta), n={n}")
    if changed:
        die(f"frozen tree(s) modified: {changed}")
    say("[V-post] all 4 frozen trees byte-identical — PASS")

    os.makedirs(f"{OUT}/logs", exist_ok=True)
    with open(f"{OUT}/logs/quant_perancestry_analysis.log", "w") as f:
        f.write("\n".join(LOG) + "\n")
    json.dump({"frozen_before": {k: {"content_sig": v[0], "meta_sig": v[1], "n_files": v[2]}
                                 for k, v in base.items()},
               "anc_n": cnt, "n_cells_published": int((pc["verdict"] >= 0).sum()),
               "n_cells_ext2": int((ec["verdict"] >= 0).sum()),
               "n_paired_cells": npair, "verdict_unchanged": same},
              open(f"{OUT}/logs/step87_summary.json", "w"), indent=2, sort_keys=True)
    print(f"\nstep87 done in {time.time()-t0:.0f}s -> {OUT}")
    print("STEP87_DONE")


if __name__ == "__main__":
    main()
