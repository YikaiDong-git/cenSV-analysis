#!/usr/bin/env python3
"""Step 80 (ext2 arm) — the q in [0.05,0.10) SUGGESTIVE-BAND chain: band -> per-ancestry meta -> exact Firth."""
import os, sys, csv
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S
import step51b_exact_firth_refit_ext2 as H
import step56_ancestry_disease_ext2 as P
import step58_ancestry_meta_ext2 as M

WELL = ["EUR", "AFR", "AMR"]
SCOPES = ["POOLED", "EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
MIN_CASES = 20
Q_LO, Q_HI = 0.05, 0.10

HITS = f"{S.OUTA}/disease_matrices_firth/disease_hits_final.tsv"
CMP = f"{S.OUTA}/comparison/perancestry"
PANC = f"{S.OUTA}/per_ancestry"
TRIAGE = f"{PANC}/triage"
BAND = f"{CMP}/phewide_disease_hits_q0510_forMETA_ext2.tsv"
META = f"{CMP}/ancestry_meta_q0510_ext2.tsv"
OUT = f"{TRIAGE}/q0510_firth_perancestry_ext2.tsv"
CONF = f"{TRIAGE}/q0510_confirmed_ext2.tsv"

BAND_COLS = ["sv_locus", "tier", "phecode", "label", "cases", "n",
             "beta_final_logOR", "p_final", "p_source", "q_bh_tier", "separation"]


def stage_a():
    os.makedirs(CMP, exist_ok=True); os.makedirs(PANC, exist_ok=True); os.makedirs(TRIAGE, exist_ok=True)
    rows, n_all, n_sep = [], 0, 0
    for r in csv.DictReader(open(HITS), delimiter="\t"):
        n_all += 1
        try:
            q = float(r["q_bh_tier"])
        except (ValueError, KeyError):
            continue
        if not (Q_LO <= q < Q_HI):
            continue
        if int(r.get("separation") or 0):
            n_sep += 1
        rows.append(r)
    with open(BAND, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=BAND_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: float(r["p_final"])):
            w.writerow({c: r.get(c, "") for c in BAND_COLS})
    print(f"[A] band 0.05<=q<0.10 : {len(rows)} cells of {n_all} in the hit table "
          f"({len({r['sv_locus'] for r in rows})} SVs x {len({r['phecode'] for r in rows})} phecodes; "
          f"{n_sep} separation-flagged)\n    -> {BAND}", flush=True)
    if not rows:
        raise SystemExit("ABORT: the ext2 suggestive band is empty — nothing to adjudicate.")
    return rows


def stage_b():
    print(f"\n[B] per-ancestry IVW meta over the band (step58 argv override)", flush=True)
    argv = sys.argv
    sys.argv = ["step58", BAND, META]
    try:
        M.main()
    finally:
        sys.argv = argv
    meta = list(csv.DictReader(open(META), delimiter="\t"))
    print(f"    -> {META} ({len(meta)} rows)", flush=True)
    return meta


def stage_c(meta, band_rows):
    sep = {(r["sv_locus"], r["phecode"]) for r in band_rows if int(r.get("separation") or 0)}
    tgt = [m for m in meta if m["verdict"].startswith("REPLICATES")]
    from collections import Counter
    print(f"\n[C] band verdict tally: {dict(Counter(m['verdict'] for m in meta))}", flush=True)
    print(f"    band tier composition:    {dict(Counter(m.get('tier', '') for m in meta))}", flush=True)
    print(f"    Firth targets (REPLICATES*): {len(tgt)} cells "
          f"({len({m['sv_locus'] for m in tgt})} SVs)  tiers "
          f"{dict(Counter(m.get('tier', '') for m in tgt))}", flush=True)
    for m in tgt:
        f_sep = "  [separation — cannot be CONFIRMED]" if (m["sv_locus"], m["phecode"]) in sep else ""
        f_null = "  ⚠️ NEGATIVE-CONTROL TIER" if m.get("tier") == "null_ext" else ""
        f_x = "  ⚠️ chrX: no ploidy handling — hypothesis only" if m["sv_locus"].startswith("chrX") else ""
        print(f"      {m.get('tier', ''):<9} {m['sv_locus'].split('_best')[0][:34]:<34} {m['phecode']:<10} "
              f"{m['disease'][:28]:<29} {m['verdict']}{f_sep}{f_null}{f_x}", flush=True)
    n_null = sum(1 for m in tgt if m.get("tier") == "null_ext")
    if n_null:
        print(f"    ⚠️ {n_null} target(s) sit on NEGATIVE-CONTROL predictors — excluded from CONFIRMED below.",
              flush=True)
    return tgt, sep


def stage_d(tgt, sep):
    if not tgt:
        print("\n[D] no REPLICATES* targets — nothing to refit.", flush=True)
        return []
    svs = sorted({t["sv_locus"] for t in tgt}); phes = sorted({t["phecode"] for t in tgt})
    tiermap = {(t["sv_locus"], t["phecode"]): t.get("tier", "") for t in tgt}
    print(f"\n[D] exact per-ancestry Firth: {len(tgt)} cells, {len(svs)} SVs x {len(phes)} phecodes", flush=True)

    iids, iidx = S.load_order(); N = len(iids)
    X, names = S.load_X(iidx)
    Xcov = X[:, 1:]
    lab = P.ancestry_labels(X, names, N).astype(str)
    dos = H.load_sv_dosages(svs, iidx)
    Y = H.load_pheno_cols(phes, iidx)
    print(f"    loaded {len(dos)} dosages, {len(Y)} pheno cols", flush=True)

    rows = []
    for t in tgt:
        sv, ph, label = t["sv_locus"], t["phecode"], t["disease"]
        g = dos.get(sv); y = Y.get(ph)
        if g is None or y is None:
            print(f"    WARN missing {sv} / {ph}", flush=True); continue
        for scope in SCOPES:
            present = (y >= 0) if scope == "POOLED" else ((y >= 0) & (lab == scope))
            nc = int((y[present] == 1).sum()); nn = int((y[present] == 0).sum())
            if nc < MIN_CASES or nn < MIN_CASES:
                rows.append((sv, ph, label, scope, nc, nn, np.nan, np.nan, False)); continue
            bf, pf, conv = H.firth(Xcov, g, y, present)
            orf = float(np.exp(bf)) if np.isfinite(bf) else np.nan
            rows.append((sv, ph, label, scope, nc, nn, orf, pf, conv))

    with open(OUT, "w") as f:
        f.write("sv_locus\tphecode\tlabel\tscope\tcases\tcontrols\tOR_firth\tp_firth\tconverged\n")
        for sv, ph, label, sc, nc, nn, orf, pf, conv in rows:
            f.write(f"{sv}\t{ph}\t{label}\t{sc}\t{nc}\t{nn}\t{'' if not np.isfinite(orf) else f'{orf:.4g}'}\t"
                    f"{'' if not np.isfinite(pf) else f'{pf:.3g}'}\t{int(bool(conv))}\n")
    print(f"    -> {OUT}", flush=True)

    cur = None
    for sv, ph, label, sc, nc, nn, orf, pf, conv in rows:
        if (sv, ph) != cur:
            cur = (sv, ph)
            print(f"\n{sv}  {ph} ({label[:38]})", flush=True)
            print(f"  {'scope':<8}{'cases':>7}{'OR':>9}{'p_firth':>11}{'conv':>6}", flush=True)
        o = f"{orf:.3f}" if np.isfinite(orf) else "na"
        p = f"{pf:.2e}" if np.isfinite(pf) else "na"
        print(f"  {sc:<8}{nc:>7}{o:>9}{p:>11}{int(bool(conv)):>6}", flush=True)

    print("\n=== per-(SV,phecode) WELL verdict (exact Firth, nominal p<0.05 within stratum) ===", flush=True)
    bycell = {}
    for sv, ph, label, sc, nc, nn, orf, pf, conv in rows:
        bycell.setdefault((sv, ph, label), {})[sc] = (orf, pf, nc, conv)
    confirmed = []
    for (sv, ph, label), d in bycell.items():
        sig = [a for a in WELL if a in d and np.isfinite(d[a][1]) and d[a][1] < 0.05
               and np.isfinite(d[a][0]) and d[a][3]]
        signs = {np.sign(np.log(d[a][0])) for a in sig}
        pooled = d.get("POOLED", (np.nan,))[0]
        conc = len(sig) >= 2 and len(signs) == 1
        blocked = (sv, ph) in sep
        tr = tiermap.get((sv, ph), "")
        if conc and tr == "null_ext":
            tag = "BLOCKED (negative-control tier — reports FPR, not replication)"
        elif conc and blocked:
            tag = "blocked (separation — Firth non-converged in the pooled fit)"
        elif conc:
            tag = "CONFIRMED(>=2 well, same sign)"; confirmed.append((sv, ph, label))
            if sv.startswith("chrX"):
                tag += "  ⚠️ chrX: no ploidy handling — hypothesis only"
        else:
            tag = "1-well-only" if len(sig) == 1 else "not-confirmed"
        po = f"{pooled:.2f}" if np.isfinite(pooled) else "na"
        print(f"  {tr:<9}{sv.split('_best')[0][:24]:<24}{ph:<10}{label[:24]:<25} sig={sig} "
              f"pooledOR={po} -> {tag}", flush=True)
    print(f"\nCONFIRMED cross-ancestry suggestive-band signals: {len(confirmed)} "
          f"({len({c[0] for c in confirmed})} distinct SVs)", flush=True)
    for sv, ph, label in confirmed:
        print(f"    {sv.split('_best')[0]:<30} {ph:<10} {label}", flush=True)

    with open(CONF, "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tsig_strata\tpooled_OR\tchrX_no_ploidy\n")
        for sv, ph, label in confirmed:
            d = bycell[(sv, ph, label)]
            sig = [a for a in WELL if a in d and np.isfinite(d[a][1]) and d[a][1] < 0.05
                   and np.isfinite(d[a][0]) and d[a][3]]
            po = d.get("POOLED", (np.nan,))[0]
            f.write(f"{sv}\t{tiermap.get((sv, ph), '')}\t{ph}\t{label}\t{','.join(sig)}\t"
                    f"{'' if not np.isfinite(po) else f'{po:.4g}'}\t{int(sv.startswith('chrX'))}\n")
    print(f"    -> {CONF}", flush=True)
    return confirmed


def main():
    print("=" * 96, flush=True)
    print("STEP 80 — ext2 suggestive band (0.05 <= q < 0.10) cross-ancestry chain", flush=True)
    print("=" * 96, flush=True)
    band = stage_a()
    meta = stage_b()
    tgt, sep = stage_c(meta, band)
    stage_d(tgt, sep)
    print("\nQ0510_EXT2_DONE", flush=True)


if __name__ == "__main__":
    main()
