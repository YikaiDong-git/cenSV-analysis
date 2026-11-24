#!/usr/bin/env python3
"""Step 88 (ext2 arm) — WHY the per-ancestry quant answer moved, and whether the Duffy refutation survives."""
import os
import sys
import gc
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np                                        # noqa: E402
from scipy.stats import fisher_exact, chi2 as _chi2       # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from frozen_tree_sig import tree_sig                      # noqa: E402

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
RV = f"{BASE}/aou_replication/results_vds"
PUB, EXT = f"{RV}/step30_phewide", f"{RV}/step30_phewide_ext2"
OUT = f"{EXT}/quant/analysis"
PUB_MAT, EXT_MAT = f"{PUB}/quant_matrices", f"{EXT}/quant/matrices"
PUB_TIER = f"{RV}/phewide_carryover/sv_tier.tsv"
EXT_TIER = f"{RV}/step28_phewide_ext2/sv_tier.tsv"
W_OLD = f"{BASE}/phewide_all_handoff"
W_NEW = f"{BASE}/phewide_all_handoff/full_11129_ext2"
W_NULL = f"{BASE}/phewide_all_handoff/discovery_null_ext"
FROZEN = ["step30_phewide", "step28_phewide", "phewide_carryover", "step30_phewide_mcc1"]

ANCS = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
WBC = ["6690-2", "770-8", "751-8"]
KNOWN = {"29463-7": "weight", "8302-2": "height", "39156-5": "BMI", "32623-1": "MPV", "776-5": "PDW"}
ACKR1_MB = 159.2
PUBLISHED_FOOTPRINT_MB = (120.2, 146.5)
LOG = []


def say(m=""):
    print(m, flush=True)
    LOG.append(m)


def die(m):
    say(f"\n*** GATE FAILED: {m}")
    say("*** aborting without writing outputs (fail-closed)")
    sys.exit(1)


def read_tier(p):
    t = {}
    with open(p) as f:
        f.readline()
        for line in f:
            q = line.rstrip("\n").split("\t")
            t[q[0]] = q[1]
    return t


def load_mat(matdir, anc, kind):
    f = f"{matdir}/quant_{kind}_{anc}.npy"
    if not os.path.exists(f):
        return None, None, None
    a = np.load(f)
    sv = [l.strip() for l in open(f"{matdir}/quant_svorder_{anc}.txt") if l.strip()]
    lo = [l.strip() for l in open(f"{matdir}/quant_loinc_{anc}.txt") if l.strip()]
    return a, sv, lo


def nlp(p):
    return -np.log10(np.clip(p.astype(np.float64), 1e-300, 1.0))


def fisher(a, b, c, d, label):
    try:
        orr, p = fisher_exact([[int(a), int(b)], [int(c), int(d)]])
    except ValueError:
        return float("nan"), float("nan")
    say(f"    {label:<58} OR={orr:>8.2f}  p={p:.3g}   [{a}/{a+b} vs {c}/{c+d}]")
    return orr, p


def h1(ptier, etier):
    say("\n" + "=" * 100)
    say("H1. IS THE PERTURBATION ANCESTRY-STRUCTURED?")
    say("    Prediction from the EUR-defined truncation: EUR moves LEAST, by a clear margin.")
    say("    All metrics are built from p and sign(beta) only - both invariant to the [0,2] rescale,")
    say("    whose scale differs between arms and makes raw beta magnitudes non-comparable.")
    rows = []
    say(f"\n  {'anc':<6}{'tier':<10}{'n cells':>10}{'r(-log10 p)':>13}{'med |dlogP|':>13}"
        f"{'frac |d|>1':>12}{'sign concord':>14}{'lamGC pub':>11}{'lamGC ext2':>11}")
    for anc in ANCS:
        Pp, psv, plo = load_mat(PUB_MAT, anc, "p")
        Pe, esv, elo = load_mat(EXT_MAT, anc, "p")
        if Pp is None or Pe is None:
            say(f"  {anc:<6}(matrix missing in one arm — stratum skipped)")
            continue
        Bp, _, _ = load_mat(PUB_MAT, anc, "beta")
        Be, _, _ = load_mat(EXT_MAT, anc, "beta")
        ei = {s: i for i, s in enumerate(esv)}
        ej = {l: j for j, l in enumerate(elo)}
        ridx = np.array([ei[s] for s in psv if s in ei], np.int64)
        keep_rows = np.array([i for i, s in enumerate(psv) if s in ei], np.int64)
        cidx = np.array([ej[l] for l in plo if l in ej], np.int64)
        keep_cols = np.array([j for j, l in enumerate(plo) if l in ej], np.int64)
        A = Pp[np.ix_(keep_rows, keep_cols)]
        B = Pe[np.ix_(ridx, cidx)]
        Ab = Bp[np.ix_(keep_rows, keep_cols)]
        Bb = Be[np.ix_(ridx, cidx)]
        svs_kept = [psv[i] for i in keep_rows]
        tv = np.array([etier.get(s, "") for s in svs_kept])
        for tier in ("signal", "null_ext"):
            sel = (tv == tier)[:, None] & np.isfinite(A) & np.isfinite(B)
            n = int(sel.sum())
            if n < 100:
                continue
            la, lb = nlp(A[sel]), nlp(B[sel])
            r = float(np.corrcoef(la, lb)[0, 1])
            d = np.abs(la - lb)
            med = float(np.median(d))
            frac = float((d > 1).mean())
            sa, sb = np.sign(Ab[sel]), np.sign(Bb[sel])
            nz = (sa != 0) & (sb != 0)
            conc = float((sa[nz] == sb[nz]).mean()) if nz.any() else float("nan")
            lam_a = float(np.median(_chi2.isf(np.clip(A[sel].astype(np.float64), 1e-300, 1), 1)) / 0.4549364231195736)
            lam_b = float(np.median(_chi2.isf(np.clip(B[sel].astype(np.float64), 1e-300, 1), 1)) / 0.4549364231195736)
            say(f"  {anc:<6}{tier:<10}{n:>10,}{r:>13.4f}{med:>13.4f}{frac:>12.4f}"
                f"{conc:>14.3f}{lam_a:>11.4f}{lam_b:>11.4f}")
            rows.append(dict(ancestry=anc, tier=tier, n_cells=n, r_neglog10p=f"{r:.5f}",
                             median_abs_dlogp=f"{med:.5f}", frac_dlogp_gt1=f"{frac:.5f}",
                             sign_concordance=f"{conc:.4f}",
                             lambda_gc_published=f"{lam_a:.4f}", lambda_gc_ext2=f"{lam_b:.4f}"))
        del Pp, Pe, Bp, Be, A, B, Ab, Bb
        gc.collect()

    sig = {r["ancestry"]: float(r["median_abs_dlogp"]) for r in rows if r["tier"] == "signal"}
    if "EUR" in sig and len(sig) > 1:
        others = {a: v for a, v in sig.items() if a != "EUR"}
        lo = min(others.values())
        say(f"\n  EUR median |dlog10 p| = {sig['EUR']:.4f}; smallest non-EUR = {lo:.4f} "
            f"({min(others, key=others.get)})")
        if sig["EUR"] < lo:
            say(f"  ⟹ H1 SUPPORTED: EUR is perturbed least of all {len(sig)} strata "
                f"({lo / max(sig['EUR'], 1e-12):.2f}x less than the next-smallest).")
            say("     Consistent with a truncation defined on a 1000G-EUR MAF>=0.01 panel: the restored")
            say("     SNPs are the ones rare or absent in Europeans, so EUR dosages barely change.")
        else:
            say("  ⟹ H1 NOT SUPPORTED: EUR is not the least-perturbed stratum. The EUR-panel mechanism")
            say("     does not explain the lambda_GC pattern; treat that pattern as unexplained.")
    return rows


def h2h3():
    say("\n" + "=" * 100)
    say("H2/H3. TESTING THE CALIBRATION CLAIMS (step87 reported the percentages; here they get a test)")
    d = {}
    p = f"{OUT}/quant_discriminates.tsv"
    if not os.path.exists(p):
        die(f"{p} missing — run step87 first")
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            q = dict(zip(hdr, line.rstrip("\n").split("\t")))
            d[(q["arm"], q["group"], q["tier"], q["p_cut"])] = q

    def get(arm, grp, tier, pcut):
        r = d.get((arm, grp, tier, pcut))
        if not r:
            return None
        n = int(r["n_cells"])
        rep = int(r["replicates_any"])
        return rep, n - rep, n

    rows = []

    def T(label, x, y):
        if not x or not y:
            say(f"    {label:<58} (a cell group is missing — skipped)")
            return
        orr, pv = fisher(x[0], x[1], y[0], y[1], label)
        rows.append(dict(test=label, a_repl=x[0], a_n=x[2], b_repl=y[0], b_n=y[2],
                         odds_ratio=f"{orr:.4g}", p_fisher=f"{pv:.4g}"))

    say("\n  H2a — anchor vs chr1->WBC (does the machinery separate them?), WITHIN each arm:")
    T("published: anthro(p<1e-4) vs chr1xWBC(p<1e-6)",
      get("published", "anthro", "ALL", "0.0001"), get("published", "chr1wbc", "ALL", "1e-06"))
    T("ext2 ALL predictors: anthro vs chr1xWBC",
      get("ext2", "anthro", "ALL", "0.0001"), get("ext2", "chr1wbc", "ALL", "1e-06"))
    T("ext2 SHARED 6,513 only: anthro vs chr1xWBC",
      get("ext2_shared", "anthro", "ALL", "0.0001"), get("ext2_shared", "chr1wbc", "ALL", "1e-06"))

    say("\n  H2b — same comparison ACROSS arms (did the separation itself change?):")
    T("anthro: published vs ext2(shared) REPLICATES rate",
      get("published", "anthro", "ALL", "0.0001"), get("ext2_shared", "anthro", "ALL", "0.0001"))
    T("chr1xWBC: published vs ext2(shared) REPLICATES rate",
      get("published", "chr1wbc", "ALL", "1e-06"), get("ext2_shared", "chr1wbc", "ALL", "1e-06"))

    say("\n  H2c — MATCHED THRESHOLDS. The published contrast held the anchor at p<1e-4 but chr1->WBC")
    say("        at p<1e-6, a stricter bar for the group it concluded against. Repeat at ONE threshold:")
    for cut, lab in (("0.0001", "both at p<1e-4"), ("1e-06", "both at p<1e-6")):
        T(f"published, {lab}", get("published", "anthro", "ALL", cut),
          get("published", "chr1wbc", "ALL", cut))
        T(f"ext2 shared, {lab}", get("ext2_shared", "anthro", "ALL", cut),
          get("ext2_shared", "chr1wbc", "ALL", cut))

    say("\n  H3 — is discovery's NEGATIVE CONTROL separable from signal on the anchor traits?")
    T("ext2 anchor: signal vs null_ext REPLICATES rate",
      get("ext2", "anthro", "signal", "0.0001"), get("ext2", "anthro", "null_ext", "0.0001"))
    T("ext2 anchor: signal vs ext2_new REPLICATES rate",
      get("ext2", "anthro", "signal", "0.0001"), get("ext2", "anthro", "ext2_new", "0.0001"))
    return rows


def read_weight_span(path, prefix_filter=None):
    lo, hi, n, svs = None, None, 0, set()
    if not os.path.exists(path):
        return None
    with open(path) as f:
        hdr = f.readline().rstrip("\r\n").split("\t")
        try:
            ip, isv = hdr.index("pos"), hdr.index("sv_locus")
        except ValueError:
            return None
        for line in f:
            q = line.rstrip("\r\n").split("\t")
            if prefix_filter and not q[isv].startswith(prefix_filter):
                continue
            try:
                v = int(q[ip])
            except (ValueError, IndexError):
                continue
            lo = v if lo is None or v < lo else lo
            hi = v if hi is None or v > hi else hi
            n += 1
            svs.add(q[isv])
    if lo is None:
        return None
    return dict(lo=lo, hi=hi, n_rows=n, n_svs=len(svs))


def c1():
    say("\n" + "=" * 100)
    say("C1. DOES THE PUBLISHED DUFFY REFUTATION SURVIVE THE CORRECTED WEIGHTS?")
    say(f"    Published claim: chr1 tag-SNP footprint {PUBLISHED_FOOTPRINT_MB[0]}-{PUBLISHED_FOOTPRINT_MB[1]} Mb;")
    say(f"    ACKR1/Duffy at ~{ACKR1_MB} Mb lies OUTSIDE it, therefore chr1->WBC is not Duffy.")
    rows = []
    cands_old = [f"{W_OLD}/sv_weights_long_chr1.tsv", f"{W_NULL}/sv_weights_long_chr1.tsv"]
    parts = []
    for p in cands_old:
        s = read_weight_span(p, "chr1_")
        if s:
            parts.append((p, s))
            say(f"    old bundle part: {p.split('phewide_all_handoff/')[-1]:<45} "
                f"{s['lo']/1e6:8.2f}-{s['hi']/1e6:8.2f} Mb  ({s['n_svs']} SVs, {s['n_rows']:,} rows)")
    if not parts:
        die("no OLD chr1 weight file found — cannot check the published footprint")
    o_lo = min(s["lo"] for _, s in parts)
    o_hi = max(s["hi"] for _, s in parts)
    o_sv = sum(s["n_svs"] for _, s in parts)
    o_n = sum(s["n_rows"] for _, s in parts)
    say(f"    OLD (truncated)  chr1 footprint: {o_lo/1e6:8.2f} - {o_hi/1e6:8.2f} Mb "
        f"({o_sv} SVs, {o_n:,} weight rows)")

    new = read_weight_span(f"{W_NEW}/sv_weights_long_ext2_chr1.tsv", "chr1_")
    if not new:
        die("ext2 chr1 weight file unreadable")
    say(f"    NEW (corrected)  chr1 footprint: {new['lo']/1e6:8.2f} - {new['hi']/1e6:8.2f} Mb "
        f"({new['n_svs']} SVs, {new['n_rows']:,} weight rows)")

    tol = 0.6
    ok_lo = abs(o_lo / 1e6 - PUBLISHED_FOOTPRINT_MB[0]) <= tol
    ok_hi = abs(o_hi / 1e6 - PUBLISHED_FOOTPRINT_MB[1]) <= tol
    if ok_lo and ok_hi:
        say(f"    [V-fp] OLD span reproduces the published {PUBLISHED_FOOTPRINT_MB[0]}-"
            f"{PUBLISHED_FOOTPRINT_MB[1]} Mb (within {tol} Mb) — recipe matches — PASS")
        licensed = True
    else:
        say(f"    [V-fp] ⚠️  OLD span does NOT reproduce the published "
            f"{PUBLISHED_FOOTPRINT_MB[0]}-{PUBLISHED_FOOTPRINT_MB[1]} Mb.")
        say("           The published figure was computed over a DIFFERENT SV subset (it described the")
        say("           chr1 SVs driving WBC, not every chr1 predictor). The comparison below is still")
        say("           valid as an OLD-vs-NEW contrast under one consistent recipe, but it does not")
        say("           re-derive the published number — stated, not glossed.")
        licensed = False

    inside_old = o_lo / 1e6 <= ACKR1_MB <= o_hi / 1e6
    inside_new = new["lo"] / 1e6 <= ACKR1_MB <= new["hi"] / 1e6
    say(f"\n    ACKR1/Duffy ({ACKR1_MB} Mb) inside OLD footprint: {inside_old}")
    say(f"    ACKR1/Duffy ({ACKR1_MB} Mb) inside NEW footprint: {inside_new}")
    if inside_new and not inside_old:
        say("    ⟹ THE REFUTATION DOES NOT SURVIVE. The corrected weights extend the chr1 tag-SNP span")
        say("       across ACKR1. 'Duffy is outside the footprint' was an artifact of the truncation and")
        say("       must be withdrawn as a refutation. This does NOT make it Duffy — it restores Duffy to")
        say("       the list of hypotheses the footprint cannot exclude.")
    elif inside_new and inside_old:
        say("    ⟹ ACKR1 is inside BOTH spans — the published footprint refutation was never supported by")
        say("       an all-chr1-predictor span; see the [V-fp] note on which subset it described.")
    else:
        say("    ⟹ THE REFUTATION SURVIVES. Even with every restored SNP, no chr1 predictor carries weight")
        say("       at ACKR1, so the Duffy explanation stays excluded on footprint grounds.")

    rows.append(dict(bundle="old_truncated", chrom="chr1", lo_mb=f"{o_lo/1e6:.3f}",
                     hi_mb=f"{o_hi/1e6:.3f}", n_svs=o_sv, n_weight_rows=o_n,
                     ackr1_inside=inside_old, reproduces_published_span=licensed))
    rows.append(dict(bundle="ext2_corrected", chrom="chr1", lo_mb=f"{new['lo']/1e6:.3f}",
                     hi_mb=f"{new['hi']/1e6:.3f}", n_svs=new["n_svs"], n_weight_rows=new["n_rows"],
                     ackr1_inside=inside_new, reproduces_published_span=""))
    return rows


def main():
    t0 = time.time()
    say("=" * 100)
    say("step88 — per-ancestry perturbation structure, calibration tests, chr1 footprint")
    say("=" * 100)
    if os.path.islink(OUT):
        die(f"{OUT} is a symlink")
    os.makedirs(OUT, exist_ok=True)
    base = {n: tree_sig(f"{RV}/{n}") for n in FROZEN}
    for n, (c, _m, k) in base.items():
        say(f"[V-pre] {n:<22} n={k:<5} content={c}")

    ptier, etier = read_tier(PUB_TIER), read_tier(EXT_TIER)
    h1rows = h1(ptier, etier)
    h2rows = h2h3()
    c1rows = c1()

    say("\n" + "=" * 100)
    say("WRITING OUTPUTS")

    def wr(name, header, rows):
        p = f"{OUT}/{name}"
        if os.path.islink(p):
            die(f"{p} is a symlink")
        with open(p, "w") as f:
            f.write("\t".join(header) + "\n")
            for r in rows:
                f.write("\t".join(str(r.get(h, "")) for h in header) + "\n")
        say(f"  {name:<40}{len(rows):>5} rows  {os.path.getsize(p):>8,} B")

    wr("quant_perturbation_by_ancestry.tsv",
       ["ancestry", "tier", "n_cells", "r_neglog10p", "median_abs_dlogp", "frac_dlogp_gt1",
        "sign_concordance", "lambda_gc_published", "lambda_gc_ext2"], h1rows)
    wr("quant_calibration_tests.tsv",
       ["test", "a_repl", "a_n", "b_repl", "b_n", "odds_ratio", "p_fisher"], h2rows)
    wr("chr1_tagsnp_footprint.tsv",
       ["bundle", "chrom", "lo_mb", "hi_mb", "n_svs", "n_weight_rows", "ackr1_inside",
        "reproduces_published_span"], c1rows)

    changed = [n for n in FROZEN if tree_sig(f"{RV}/{n}")[0] != base[n][0]]
    if changed:
        die(f"frozen tree(s) modified: {changed}")
    say("[V-post] all 4 frozen trees byte-identical — PASS")
    os.makedirs(f"{OUT}/logs", exist_ok=True)
    with open(f"{OUT}/logs/quant_perturbation.log", "w") as f:
        f.write("\n".join(LOG) + "\n")
    print(f"\nstep88 done in {time.time()-t0:.0f}s")
    print("STEP88_DONE")


if __name__ == "__main__":
    main()
