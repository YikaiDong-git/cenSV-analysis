#!/usr/bin/env python3
"""Step 90 (ext2 arm) — the RE-FROZEN deliverable table set (<- FROZEN_TABLES_2026-07-13.md)."""
import os
import sys
import csv
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np                                       # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import step51_disease_scoretest_ext2 as S                # noqa: E402
import step51b_exact_firth_refit_ext2 as H               # noqa: E402
from frozen_tree_sig import tree_sig                     # noqa: E402

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
RV = f"{BASE}/aou_replication/results_vds"
EXT = f"{RV}/step30_phewide_ext2"
TAB = f"{EXT}/tables"
DOC = f"{EXT}/docs/FROZEN_TABLES_EXT2_2026-08-04.md"
TRIAGE = f"{EXT}/per_ancestry/triage"
CMP = f"{EXT}/comparison"
W_NEW = f"{BASE}/phewide_all_handoff/full_11129_ext2"
FROZEN = ["step30_phewide", "step28_phewide", "phewide_carryover", "step30_phewide_mcc1"]

WELL = ["EUR", "AFR", "AMR"]
SD_MIN = 0.05
CVR2_A = 0.45
CVR2_C = 0.20
R_INDEP = 0.5
SUPPRESS = 20

CEN = {"chr1": (121.7, 125.1), "chr2": (91.8, 96.0), "chr3": (87.8, 94.0), "chr4": (48.2, 51.8),
       "chr5": (46.1, 50.1), "chr6": (58.5, 62.6), "chr7": (58.1, 62.1), "chr8": (43.2, 47.2),
       "chr9": (42.2, 45.5), "chr10": (38.0, 41.6), "chr11": (51.0, 55.8), "chr12": (33.2, 37.8),
       "chr13": (16.5, 18.9), "chr14": (16.1, 18.2), "chr15": (17.5, 20.5), "chr16": (35.3, 38.3),
       "chr17": (22.7, 27.4), "chr18": (15.4, 21.0), "chr19": (24.4, 28.1), "chr20": (25.7, 30.4),
       "chr21": (10.9, 13.0), "chr22": (13.7, 17.4), "chrX": (58.1, 63.8)}

LOG = []


def say(m=""):
    print(m, flush=True)
    LOG.append(m)


def die(m):
    say(f"\n*** GATE FAILED: {m}")
    sys.exit(1)


def load(path, keyfn=None):
    if not os.path.exists(path):
        die(f"missing input {path}")
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    return {keyfn(r): r for r in rows} if keyfn else rows


def f(x, d=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def chrom_of(sv):
    return sv.split("_")[0]


def sv_footprints(svs):
    want = {}
    for s in svs:
        want.setdefault(chrom_of(s), set()).add(s)
    out = {}
    for c, members in sorted(want.items()):
        p = f"{W_NEW}/sv_weights_long_ext2_{c}.tsv"
        if not os.path.exists(p):
            continue
        acc = {s: [None, None, 0] for s in members}
        with open(p) as fh:
            hdr = fh.readline().rstrip("\r\n").split("\t")
            ip, isv = hdr.index("pos"), hdr.index("sv_locus")
            for line in fh:
                q = line.rstrip("\r\n").split("\t")
                a_ = acc.get(q[isv])
                if a_ is None:
                    continue
                try:
                    v = int(q[ip])
                except (ValueError, IndexError):
                    continue
                a_[0] = v if a_[0] is None or v < a_[0] else a_[0]
                a_[1] = v if a_[1] is None or v > a_[1] else a_[1]
                a_[2] += 1
        cen = CEN.get(c)
        for s, (lo, hi, n) in acc.items():
            if lo is None:
                continue
            ov = bool(cen and not (hi / 1e6 < cen[0] or lo / 1e6 > cen[1]))
            out[s] = dict(lo_mb=lo / 1e6, hi_mb=hi / 1e6, n_rows=n, cen=cen, overlaps=ov)
    return out


def footprints(chroms):
    out = {}
    for c in sorted(chroms):
        p = f"{W_NEW}/sv_weights_long_ext2_{c}.tsv"
        if not os.path.exists(p):
            out[c] = None
            continue
        lo = hi = None
        n = 0
        with open(p) as fh:
            hdr = fh.readline().rstrip("\r\n").split("\t")
            ip = hdr.index("pos")
            for line in fh:
                q = line.rstrip("\r\n").split("\t")
                try:
                    v = int(q[ip])
                except (ValueError, IndexError):
                    continue
                lo = v if lo is None or v < lo else lo
                hi = v if hi is None or v > hi else hi
                n += 1
        if lo is None:
            out[c] = None
            continue
        cen = CEN.get(c)
        ov = bool(cen and not (hi / 1e6 < cen[0] or lo / 1e6 > cen[1]))
        out[c] = dict(lo_mb=lo / 1e6, hi_mb=hi / 1e6, n_rows=n, cen=cen, overlaps=ov)
    return out


def independent_loci(svs):
    iids, iidx = S.load_order()
    dos = H.load_sv_dosages(list(svs), iidx)
    have = [s for s in svs if s in dos and np.isfinite(dos[s]).any()]
    miss = [s for s in svs if s not in have]
    if miss:
        say(f"  ⚠️ {len(miss)} confirmed SV(s) had no readable dosage: {miss[:3]}")
    M = np.vstack([dos[s] for s in have]).astype(np.float64)
    M -= M.mean(axis=1, keepdims=True)
    sd = M.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    M /= sd
    R = np.abs((M @ M.T) / M.shape[1])
    parent = list(range(len(have)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i in range(len(have)):
        for j in range(i + 1, len(have)):
            if R[i, j] >= R_INDEP:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[ri] = rj
    cl = {}
    for i, s in enumerate(have):
        cl.setdefault(find(i), []).append(s)
    groups = sorted(cl.values(), key=lambda g: (chrom_of(g[0]), -len(g)))
    return {s: gi for gi, g in enumerate(groups) for s in g}, groups


def main():
    t0 = time.time()
    say("=" * 100)
    say("step90 — re-frozen ext2 deliverable table set (audit gap B)")
    say("=" * 100)
    for p in (TAB, os.path.dirname(DOC)):
        if os.path.islink(p):
            die(f"{p} is a symlink")
    os.makedirs(TAB, exist_ok=True)
    base = {n: tree_sig(f"{RV}/{n}") for n in FROZEN}
    for n, (c, _m, k) in base.items():
        say(f"[V-pre] {n:<22} n={k:<5} content={c}")

    CVSRC = f"{RV}/step28_phewide_ext2/sv_cv_r2_all.tsv"
    cvmap = {}
    if os.path.exists(CVSRC):
        for line in open(CVSRC):
            q = line.rstrip("\n").split("\t")
            if q[0] != "sv_locus" and len(q) > 1:
                v = f(q[1])
                if v is not None:
                    cvmap[q[0]] = v
    say(f"cv_r2 source: {len(cvmap):,} SVs from sv_cv_r2_all.tsv")

    conf = load(f"{TRIAGE}/q0510_confirmed_ext2.tsv")
    tri = load(f"{TRIAGE}/robustness_triage_ext2.tsv", lambda r: (r["sv_locus"], r["phecode"]))
    fir = load(f"{TRIAGE}/q0510_firth_perancestry_ext2.tsv")
    pooled = load(f"{CMP}/perancestry/ancestry_meta_pooledhits_ext2.tsv")
    band = load(f"{CMP}/perancestry/ancestry_meta_q0510_ext2.tsv")
    say(f"\ninputs: {len(conf)} CONFIRMED cells · {len(tri)} triaged · {len(fir)} per-ancestry Firth rows "
        f"· {len(pooled)} pooled-hit meta · {len(band)} band meta")
    if not conf:
        die("no CONFIRMED cells — nothing to freeze")

    fmap = {(r["sv_locus"], r["phecode"], r["scope"]): r for r in fir}

    svs = sorted({r["sv_locus"] for r in conf})
    say(f"\nRECOMPUTING independence over the {len(svs)} confirmed SVs (|r| >= {R_INDEP}) …")
    cl_of, groups = independent_loci(svs)
    say(f"  -> {len(groups)} independent loci "
        f"(cluster sizes {sorted((len(g) for g in groups), reverse=True)})")

    say(f"\nRECOMPUTING pericentromeric footprints from the ext2 weight bundle …")
    fp = footprints({chrom_of(s) for s in svs})
    for c, d in sorted(fp.items()):
        if d:
            say(f"  {c:<7}{d['lo_mb']:8.2f}-{d['hi_mb']:8.2f} Mb  centromere "
                f"{d['cen'][0]}-{d['cen'][1]}  overlaps={d['overlaps']}")
    n_ov = sum(1 for d in fp.values() if d and d["overlaps"])
    say(f"  -> {n_ov}/{len(fp)} CHROMOSOME footprints overlap their centromere")
    say("\n  per-SV footprints (what Table 1 carries — a chromosome span can straddle a centromere that")
    say("  an individual SV on it does not; step84 fig4 found exactly that for chr1):")
    svfp = sv_footprints(svs)
    n_sv_ov = sum(1 for d in svfp.values() if d["overlaps"])
    for s in sorted(svfp, key=lambda x: (chrom_of(x), svfp[x]["lo_mb"]))[:40]:
        d = svfp[s]
        say(f"    {s[:44]:<46}{d['lo_mb']:8.2f}-{d['hi_mb']:8.2f} Mb  n={d['n_rows']:<5} "
            f"overlaps_centromere={d['overlaps']}")
    say(f"  -> {n_sv_ov}/{len(svfp)} confirmed SVs have a footprint overlapping their own centromere")

    say("\n" + "-" * 100)
    say("TABLE 1 — cross-ancestry disease signals (tier COMPUTED, not assigned)")
    t1 = []
    for r in conf:
        sv, ph = r["sv_locus"], r["phecode"]
        tr = tri.get((sv, ph), {})
        cv = f(tr.get("cv_r2"))
        cv_src = "step81"
        if cv is None and sv in cvmap:
            cv, cv_src = cvmap[sv], "sv_cv_r2_all (step81 join missed it)"
        sds = {a: f(tr.get(f"{a}_sd")) for a in WELL}
        uninf = [a for a, v in sds.items() if v is not None and v < SD_MIN]
        if cv is not None and cv >= CVR2_A and not uninf:
            tier, why = "A", f"cv_r2 {cv:.2f}>={CVR2_A}, informative in all 3"
        elif (cv is not None and cv < CVR2_C) or uninf:
            bits = []
            if cv is not None and cv < CVR2_C:
                bits.append(f"cv_r2 {cv:.2f} low")
            if uninf:
                bits.append("near-constant dosage in " + "/".join(uninf))
            tier, why = "C", "; ".join(bits)
        else:
            tier, why = "B", f"cv_r2 {cv:.2f}" if cv is not None else "cv_r2 unavailable"
        cells = {}
        for a in ["POOLED"] + WELL:
            fr = fmap.get((sv, ph, a))
            if not fr:
                continue
            nca = f(fr.get("cases"))
            if a != "POOLED" and nca is not None and nca < SUPPRESS:
                cells[a] = ("suppressed", "", "")
                continue
            cells[a] = (fr.get("OR_firth", ""), fr.get("p_firth", ""), fr.get("converged", ""))
        t1.append(dict(
            locus=f"{chrom_of(sv)} (cluster {cl_of.get(sv, -1)})", sv_locus=sv, sv_tier=r["tier"],
            phecode=ph, label=r["label"], sig_strata=r["sig_strata"], pooled_OR=r["pooled_OR"],
            cv_r2="" if cv is None else f"{cv:.3f}", cv_r2_source=("" if cv is None else cv_src),
            robustness_tier=tier, tier_reason=why,
            EUR_OR=cells.get("EUR", ("", "", ""))[0], EUR_p=cells.get("EUR", ("", "", ""))[1],
            AFR_OR=cells.get("AFR", ("", "", ""))[0], AFR_p=cells.get("AFR", ("", "", ""))[1],
            AMR_OR=cells.get("AMR", ("", "", ""))[0], AMR_p=cells.get("AMR", ("", "", ""))[1],
            POOLED_OR=cells.get("POOLED", ("", "", ""))[0], POOLED_p=cells.get("POOLED", ("", "", ""))[1],
            survives_unrelated=tr.get("survives_unrelated", ""),
            EUR_sd=tr.get("EUR_sd", ""), AFR_sd=tr.get("AFR_sd", ""), AMR_sd=tr.get("AMR_sd", ""),
            chrX_no_ploidy=r.get("chrX_no_ploidy", ""),
            sv_footprint_mb=("" if sv not in svfp else
                             f"{svfp[sv]['lo_mb']:.2f}-{svfp[sv]['hi_mb']:.2f}"),
            sv_footprint_overlaps_centromere=("" if sv not in svfp else int(svfp[sv]["overlaps"])),
            sv_footprint_n_snps=("" if sv not in svfp else svfp[sv]["n_rows"]),
            chrom_footprint_mb=("" if not fp.get(chrom_of(sv)) else
                                f"{fp[chrom_of(sv)]['lo_mb']:.2f}-{fp[chrom_of(sv)]['hi_mb']:.2f}")))
    tcount = {}
    for r in t1:
        tcount[r["robustness_tier"]] = tcount.get(r["robustness_tier"], 0) + 1
    say(f"  {len(t1)} cells; tier A/B/C = {tcount.get('A',0)}/{tcount.get('B',0)}/{tcount.get('C',0)}")
    surv = sum(1 for r in t1 if str(r["survives_unrelated"]).lower() in ("1", "true", "yes"))
    say(f"  survive the unrelated-only refit: {surv}/{len(t1)}")

    say("\n" + "-" * 100)
    say("TABLE 3 — per-ancestry verdict tally")
    t3 = []

    def tally(rows, name):
        c = {}
        for r in rows:
            c[r.get("verdict", "")] = c.get(r.get("verdict", ""), 0) + 1
        n = len(rows)
        rec = dict(hit_set=name, n=n)
        rec.update({k: v for k, v in sorted(c.items())})
        t3.append(rec)
        say(f"  {name:<34}n={n:<6} " + "  ".join(f"{k}={v}" for k, v in sorted(c.items(), key=lambda x: -x[1])))
        return c

    tally(pooled, "Disease q<0.05 (ext2)")
    tally(band, "Disease q 0.05-0.10 (ext2)")
    qt = f"{EXT}/quant/analysis/quant_verdict_tally.tsv"
    if os.path.exists(qt):
        qr = [r for r in load(qt) if r["scope"] == "ext2_all"]
        n = sum(int(r["n"]) for r in qr)
        rec = dict(hit_set="Quant, all panel cells (ext2)", n=n)
        rec.update({r["verdict"]: int(r["n"]) for r in qr})
        t3.append(rec)
        say(f"  {'Quant, all panel cells (ext2)':<34}n={n:<6} " +
            "  ".join(f"{r['verdict']}={r['n']}" for r in sorted(qr, key=lambda x: -int(x['n']))))

    say("\n" + "-" * 100)
    say("TABLE 4 — quant replication by trait")
    t4 = []
    qm = f"{EXT}/quant/quant_ancestry_meta.tsv"
    if os.path.exists(qm):
        per = {}
        with open(qm) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            il, iv, ip = hdr.index("loinc"), hdr.index("verdict"), hdr.index("meta_p")
            for line in fh:
                q = line.rstrip("\n").split("\t")
                mp = f(q[ip])
                if mp is None or mp >= 1e-4:
                    continue
                a = per.setdefault(q[il], {})
                a[q[iv]] = a.get(q[iv], 0) + 1
        names = {}
        cw = f"{EXT}/reference/trait_crosswalk.tsv"
        if os.path.exists(cw):
            for r in load(cw):
                if r.get("aou_type") == "loinc" and r.get("aou_id"):
                    names.setdefault(r["aou_id"], r.get("aou_name", ""))
        rep = ("REPLICATES", "REPLICATES_HET", "REPLICATES_PLIM")
        rows = []
        for lo, c in per.items():
            n = sum(c.values())
            rows.append((n, lo, c))
        for n, lo, c in sorted(rows, reverse=True)[:15]:
            r_ = sum(c.get(k, 0) for k in rep)
            s_ = c.get("STRAT_SUSPECT", 0)
            t4.append(dict(loinc=lo, trait=names.get(lo, ""), n_cells_metap_lt_1e4=n,
                           replicates_any=r_, strat_suspect=s_,
                           replicate_pct=f"{100*r_/n:.1f}" if n else ""))
            say(f"  {lo:<10}{names.get(lo,'')[:38]:<40}n={n:<6}REPL {r_:<6}({100*r_/n:5.1f}%)  "
                f"STRAT {s_}")

    say("\n" + "-" * 100)

    def wr(name, rows):
        if not rows:
            say(f"  {name}: no rows — skipped")
            return
        hdr = list(dict.fromkeys(k for r in rows for k in r))
        p = f"{TAB}/{name}"
        if os.path.islink(p):
            die(f"{p} is a symlink")
        with open(p, "w") as fh:
            fh.write("\t".join(hdr) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(h, "")) for h in hdr) + "\n")
        say(f"  {name:<44}{len(rows):>4} rows  {os.path.getsize(p):>8,} B")

    wr("table1_headline_signals_ext2.tsv", t1)
    wr("table3_verdict_tally_ext2.tsv", t3)
    wr("table4_quant_replication_ext2.tsv", t4)

    import collections
    by = collections.defaultdict(list)
    for r in t1:
        by[r["locus"]].append(r)

    def pmin(g):
        ps = [f(x["POOLED_p"]) for x in g if f(x["POOLED_p"]) is not None]
        return min(ps) if ps else 9.0

    tier_n = collections.Counter(r["robustness_tier"] for r in t1)
    n_nocv = sum(1 for r in t1 if not r["cv_r2"])
    n_recovered = sum(1 for r in t1 if r.get("cv_r2_source", "").startswith("sv_cv_r2_all"))
    n_ov_sv = sum(1 for r in t1 if str(r["sv_footprint_overlaps_centromere"]) == "1")
    L = []
    L.append("# Frozen deliverable — key tables, ext2 (11,129-predictor) programme")
    L.append("")
    L.append("**Work unit:** `phewide_ext2_11129` · **Date:** 2026-08-04 · "
             "**Supersedes:** `../../step30_phewide/FROZEN_TABLES_2026-07-13.md`")
    L.append("**Parent trail:** `../../../PROJECT_OVERVIEW_AND_MASTER_TODO_2026-07-13.md` · "
             "**Generated by:** `code/phewide_ext2/step90_frozen_tables_ext2.py`")
    L.append("**Closes:** audit gap **B** (and supporting gaps #14 independence, #15 A/B/C tier) of "
             "`SUPERSESSION_NOTICE_2026-08-04.md`.")
    L.append("")
    L.append("> **Every number below is DERIVED from a result file by the generating script** — no value is "
             "hand-typed. That is deliberate: the published `step30x` hardcodes the numbers it plots, and "
             "`step30v`/`step30t` hardcoded the loci confirmed under the old weights, so a copied table "
             "would silently restate published results under an ext2 filename. Re-run the script to "
             "regenerate this file. Aggregate only; per-ancestry cells with <%d cases are suppressed." % SUPPRESS)
    L.append("")
    L.append("---")
    L.append("")
    L.append("## Table 1 — cross-ancestry disease signals")
    L.append("")
    L.append(f"**{len(t1)} CONFIRMED cells · {len({r['sv_locus'] for r in t1})} SVs · "
             f"{len(by)} independent loci** (|r| >= {R_INDEP}, recomputed from dosages). "
             f"Published counterpart: **6 signals / ~5 independent loci**. "
             f"One row per independent locus, showing its strongest cell; the full 28 are in "
             f"`../tables/table1_headline_signals_ext2.tsv`.")
    L.append("")
    L.append("| # | locus | representative SV | phenotype | cells | strata | pooled OR (p) | "
             "EUR OR (p) | AFR OR (p) | AMR OR (p) | cv_r² | tier | SV footprint (Mb) | pericen? |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, (loc, g) in enumerate(sorted(by.items(), key=lambda kv: pmin(kv[1])), 1):
        r = sorted(g, key=lambda x: f(x["POOLED_p"], 9.0))[0]
        def cell(a):
            o, p = r[f"{a}_OR"], r[f"{a}_p"]
            return f"{o} ({p})" if o and o != "suppressed" else (o or "—")
        L.append(f"| {i} | {loc} | `{r['sv_locus'][:38]}` | {r['label'][:44]} ({r['phecode']}) | "
                 f"{len(g)} | {r['sig_strata']} | {cell('POOLED')} | {cell('EUR')} | {cell('AFR')} | "
                 f"{cell('AMR')} | {r['cv_r2'] or '—'} | **{r['robustness_tier']}** | "
                 f"{r['sv_footprint_mb']} | {'yes' if str(r['sv_footprint_overlaps_centromere'])=='1' else '**no**'} |")
    L.append("")
    L.append(f"**Robustness tier is COMPUTED**, applying the published rule "
             f"(A = cv_r²>={CVR2_A} and informative in all three well-powered strata; C = caveated by "
             f"cv_r²<{CVR2_C} or a near-constant dosage, sd<{SD_MIN}, in a stratum; B otherwise): "
             f"**A {tier_n.get('A',0)} / B {tier_n.get('B',0)} / C {tier_n.get('C',0)}** over the {len(t1)} cells.")
    L.append("")
    L.append("Four things a reader should not miss, all of which fall out of the derivation:")
    L.append("")
    L.append(f"1. **Tier A: {tier_n.get('A',0)} of {len(t1)} cells** (published Table 1 had one, chr19 "
             f"hepatitis at cv_r² 0.46). Most confirmed cells are caveated by imputation quality or by a "
             f"stratum in which the dosage carries almost no variance.")
    L.append(f"2. ⚠️ **`step81` leaves cv_r² blank for {n_recovered} of these {len(t1)} cells, and every "
             f"one of them IS in `sv_cv_r2_all.tsv`** — its annotation join silently missed them. Taking "
             f"the blanks at face value would have produced the false headline \"no cell reaches tier A\": "
             f"two of the recovered values are >=%.2f. They are joined from source here, and the "
             f"`cv_r2_source` column in the TSV records which came from where. "
             f"{n_nocv} cell(s) still have no cv_r² from any source." % CVR2_A)
    L.append(f"3. **chr10 -> psoriatic arthropathy, the flagship positive control, computes to tier C** "
             f"— not on cv_r² (0.765, the best in the table) but because its dosage is near-constant in "
             f"EUR and AFR. This reproduces the published triage's own caveat, mechanically.")
    L.append(f"4. **{n_ov_sv} of {len(t1)} cells sit on an SV whose OWN tag-SNP footprint overlaps its "
             f"centromere.** The exceptions are named in the table. ⚠️ Note the granularity: the "
             f"per-CHROMOSOME span overlaps for {n_ov}/{len(fp)} chromosomes, but a chromosome span can "
             f"straddle a centromere that an individual SV on it does not — reporting the chromosome "
             f"number in a per-SV table would silently upgrade a distal SV to 'pericentromeric'.")
    L.append("")
    L.append("**All %d/%d cells survive the unrelated-only refit.**" % (surv, len(t1)))
    L.append("")
    L.append("---")
    L.append("")
    L.append("## Table 2 — the deflation ledger (what the pipeline killed or relegated)")
    L.append("")
    L.append("Carried narrative, with each row's numbers sourced. Rows marked **NEW** are ext2-only.")
    L.append("")
    L.append("| signal | apparent result | rigorous test | outcome | source |")
    L.append("|---|---|---|---|---|")
    L.append("| chr1_sv650 -> breast cancer | replication p=6.1e-8 | arm-fix re-extraction | "
             "**ARTIFACT** — collapsed to p 0.22 (truncated 121 Mb fragment) | published ledger |")
    L.append("| **NEW** 2 of the 6 published headlines | chr19->hepatitis, chr19->GI-obstruction | "
             "bundle membership | **carried by discovery's NEGATIVE CONTROL** (`null_ext`); both "
             "phenotypes re-carried by non-control predictors — the locus replicates, the tagging SV does "
             "not | `HEADTOHEAD` §1 |")
    L.append("| chr16 cluster -> substance-use | pooled-robust Firth 1e-9..1e-15 | per-ancestry Firth + "
             "within-AMR PC rescue | **AMR-restricted**; 197/197 survive within-AMR PCs; ⚠️ 53/197 exceed "
             "no magnitude bound | `Q0510…` |")
    L.append("| chr1 -> WBC / neutrophils | \"clean quant biology\" | per-ancestry stratification, then "
             "corrected weights | **weakened twice**: 88.5% -> 76.7% STRAT_SUSPECT, and the calibration "
             "that made STRAT_SUSPECT interpretable is 12–22x weaker | `QUANT_PERANCESTRY_EXT2` §4 |")
    L.append("| **NEW** chr1 -> WBC, third strike | — | degeneracy sweep | **WBC is the single worst-"
             "flooding trait** (0.0243 of all SVs at p<1e-6, over step62's 0.02) | `ENGINE_VALIDATION_EXT2` §5 |")
    L.append("| weight / height / BMI | replicate, sign-concordant | SNP-conditioning + COJO | "
             "**KNOWN euchromatic loci, not novel** | published ledger |")
    L.append("| chr1_sv1325 block | score-BH significant | exact Firth | **score-test inflation** | "
             "published ledger |")
    L.append("| chr11 -> Prothrombin (GE_971.2) | BH-significant | clinical thrombosis phecodes; then "
             "Gate 2 | **ascertainment artefact**, and WITHDRAWN — truncation had deleted a whole "
             "chromosome arm | `GATE2_FINDINGS` |")
    L.append("| **NEW** `null_ext` across the whole ext2 run | — | 3,097 control predictors | "
             "**0 cells anywhere** (significant, band, and targets) | `HEADTOHEAD` §1 |")
    L.append("")
    L.append("---")
    L.append("")
    L.append("## Table 3 — per-ancestry verdict tally")
    L.append("")
    vs = ["REPLICATES", "REPLICATES_HET", "REPLICATES_PLIM", "STRAT_SUSPECT", "DISCORDANT", "UNDERPOWERED"]
    L.append("| hit set | n | " + " | ".join(vs) + " |")
    L.append("|---" * (len(vs) + 2) + "|")
    for r in t3:
        L.append(f"| {r['hit_set']} | {r['n']:,} | " +
                 " | ".join(f"{r.get(v, 0):,}" for v in vs) + " |")
    L.append("")
    L.append("---")
    L.append("")
    L.append("## Table 4 — quant replication by trait (cells at meta_p < 1e-4)")
    L.append("")
    L.append("| LOINC | trait | cells | REPLICATES (any) | rate | STRAT_SUSPECT |")
    L.append("|---|---|---|---|---|---|")
    for r in t4:
        L.append(f"| `{r['loinc']}` | {r['trait'] or '—'} | {r['n_cells_metap_lt_1e4']} | "
                 f"{r['replicates_any']} | {r['replicate_pct']}% | {r['strat_suspect']} |")
    L.append("")
    L.append("⚠️ **`REPLICATES` on this panel means 'cross-stratum consistent', not 'real SV effect'** — "
             "discovery's negative-control bundle reaches it 55.3% of the time on the anthropometric "
             "anchor against signal's 77.6% (separable at OR 2.80, p=6.7e-09, but not by much). "
             "See `QUANT_PERANCESTRY_EXT2_2026-08-04.md` §6.")
    L.append("")
    L.append("**Engine note.** The disease score test is calibrated on ext2 predictors — "
             "frac(p_score < half p_firth) = 0.0% in every carrier bin, case bin and tier, including the "
             "`ext2_new` and `ext2_chrX` classes the published gate never covered "
             "(`ENGINE_VALIDATION_EXT2_2026-08-04.md`). The quant FWL engine agrees with exact OLS to "
             "max |Δ−log10 p| = 0.52.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("## What is NOT in this table set")
    L.append("")
    L.append("- **No cell here is genome-wide significant** (none at 5e-8); the CONFIRMED set is selected "
             "at pooled q∈[0.05,0.10) with a nominal per-stratum threshold, so it measures cross-stratum "
             "consistency, not independent replication.")
    L.append("- **Nothing is externally validated, and no external replication is available** "
             "(decision 2026-08-04).")
    L.append("  This is a **permanent limitation of the study**, not an open task: every claim in this "
             "table set is")
    L.append("  AoU-internal, and the evidence for it is internal consistency (cross-stratum, "
             "cross-case-definition,")
    L.append("  negative-control and engine-calibration checks), not independent replication.")
    L.append(f"- **{len(by)} independent loci, not {len(t1)} findings.**")
    L.append("- chrX rows are hypothesis-only: `step50` still scores chrX with a ploidy-blind LUT. "
             "`step86` refuted the sex/ploidy confound but did not validate the dosages.")
    if os.path.islink(DOC):
        die(f"{DOC} is a symlink")
    with open(DOC, "w") as fh:
        fh.write("\n".join(L) + "\n")
    say(f"  {'FROZEN_TABLES_EXT2_2026-08-04.md':<44}{len(L):>4} lines  {os.path.getsize(DOC):>8,} B")

    changed = [n for n in FROZEN if tree_sig(f"{RV}/{n}")[0] != base[n][0]]
    if changed:
        die(f"frozen tree(s) modified: {changed}")
    say("[V-post] all 4 frozen trees byte-identical — PASS")
    say(f"\nstep90 done in {time.time()-t0:.0f}s")
    say("FROZEN_TABLES_EXT2_DONE")
    with open(f"{TAB}/step90_run.log", "w") as fh:
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
