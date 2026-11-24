#!/usr/bin/env python3
"""Step 93 (ext2 arm) — chrX PLOIDY: establish the right model BEFORE changing step50."""
import os
import sys
import csv
import json
import time
import argparse

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phewide_ext"))
import step42_gate2_headline_rescore as G          # noqa: E402
import step51_disease_scoretest_ext2 as S          # noqa: E402
import step51b_exact_firth_refit_ext2 as F         # noqa: E402

BASE = G.BASE
RV = G.RV
BED = G.BED
CARRY = G.CARRY
XPFX = f"{BED}/chrX/chrX.pericentromeric_snps.min_ac1000"
WX = f"{BASE}/phewide_all_handoff/full_11129_ext2/sv_weights_long_ext2_chrX.tsv"
OUT = f"{RV}/step30_phewide_ext2/chrx_ploidy"

PAR1_END, PAR2_START = 2_781_479, 155_701_383
MODELS = ("XCI", "HAP", "SEXSTD")

LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def infer_sex(fam, bim, bpv, Nfam, n_probe=400):
    HET_MIN, SCAN = 0.05, 3000
    scan = np.linspace(0, len(bim) - 1, min(SCAN, len(bim))).astype(int)
    informative = []
    with open(XPFX + ".bed", "rb") as bf:
        assert bf.read(3) == b"\x6c\x1b\x01", "bad .bed magic"
        for vi in scan:
            bf.seek(3 + int(vi) * bpv)
            a1c = G._LUT[np.frombuffer(bf.read(bpv), dtype=np.uint8)].reshape(-1)[:Nfam]
            ok = ~np.isnan(a1c)
            if ok.sum() and (a1c[ok] == 1.0).mean() >= HET_MIN:
                informative.append(int(vi))
            if len(informative) >= n_probe:
                break
        assert informative, f"no chrX SNP reaches pooled het >= {HET_MIN}; cannot infer sex"

        het = np.zeros(Nfam, np.float64)
        obs = np.zeros(Nfam, np.float64)
        for vi in informative:
            bf.seek(3 + vi * bpv)
            a1c = G._LUT[np.frombuffer(bf.read(bpv), dtype=np.uint8)].reshape(-1)[:Nfam]
            ok = ~np.isnan(a1c)
            obs += ok
            het += (a1c == 1.0)
    rate = np.divide(het, obs, out=np.full(Nfam, np.nan), where=obs > 0)
    print(f"  probe SNPs with pooled het >= {HET_MIN}: {len(informative)}", flush=True)
    fin = rate[np.isfinite(rate)]
    grid = np.linspace(0.005, max(0.06, float(np.percentile(fin, 60))), 60)
    counts = np.array([((fin >= g - 0.004) & (fin < g + 0.004)).sum() for g in grid])
    thr = float(grid[int(np.argmin(counts))])
    print(f"  data-derived male/female het threshold: {thr:.4f}", flush=True)
    inferred = np.where(rate < thr, 1.0, 0.0)
    inferred[np.isnan(rate)] = np.nan
    return inferred, rate, obs


def validate_sex(fam, inferred, rate):
    known = {}
    with open(f"{CARRY}/inputs/covar.txt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        si, ii = hdr.index("sex"), hdr.index("IID")
        for line in f:
            p = line.rstrip("\n").split("\t")
            known[p[ii]] = float(p[si])
    pos = [(i, known[x]) for i, x in enumerate(fam) if x in known]
    a = np.array([inferred[i] for i, _ in pos])
    b = np.array([k for _, k in pos])
    rr = np.array([rate[i] for i, _ in pos])
    ok = np.isfinite(a)
    agree = float((a[ok] == b[ok]).mean())
    say(f"  sex inference vs covariate on {int(ok.sum()):,} cohort samples: agreement {agree * 100:.4f}%")
    say(f"    inferred male {int((a[ok] == 1).sum()):,} / covariate male {int((b[ok] == 1).sum()):,}")

    disc = ok & (a != b)
    nd = int(disc.sum())
    say(f"    discordant: {nd:,} ({100.0 * nd / int(ok.sum()):.3f}%)")
    if nd:
        mmode = float(np.median(rr[ok & (a == 1)]))
        fmode = float(np.median(rr[ok & (a == 0)]))
        dmed = float(np.median(rr[disc]))
        between = float(((rr[disc] > mmode) & (rr[disc] < fmode)).mean())
        say(f"    het-rate modes: male {mmode:.4f}  female {fmode:.4f}  |  discordant median {dmed:.4f}")
        say(f"    discordant samples lying BETWEEN the two modes: {between * 100:.1f}% "
            f"({'consistent with aneuploidy/mosaicism' if between > 0.5 else 'NOT mode-intermediate'})")
    assert agree > 0.995, (f"chrX-heterozygosity sex inference agrees only {agree:.4%} with the covariate "
                           f"table — too low to be aneuploidy; the inference itself is suspect")
    return agree, nd


def score_models(rows_by_sv, bim_pos, bpv, Nfam, male, cols):
    pos_map = {}
    for vi, (p, a1, a2) in enumerate(bim_pos):
        pos_map.setdefault(p, []).append((vi, a1, a2))

    isM = male == 1.0
    isF = male == 0.0
    out = {m: {} for m in MODELS}
    qc = []
    bf = open(XPFX + ".bed", "rb")
    assert bf.read(3) == b"\x6c\x1b\x01"

    for sv, rws in rows_by_sv.items():
        acc = {m: np.zeros(Nfam, np.float64) for m in MODELS}
        n_used = 0
        tot = sum(abs(float(r["beta_std"])) for r in rws)
        miss = 0.0
        for r in rws:
            b = abs(float(r["beta_std"]))
            hit = None
            for (vi, a1, a2) in pos_map.get(int(r["pos"]), []):
                if a1 == r["alt"] and a2 == r["ref"]:
                    hit = (vi, False); break
                if a1 == r["ref"] and a2 == r["alt"]:
                    hit = (vi, True); break
            if hit is None:
                miss += b; continue
            vi, sw = hit
            bf.seek(3 + vi * bpv)
            a1c = G._LUT[np.frombuffer(bf.read(bpv), dtype=np.uint8)].reshape(-1)[:Nfam]
            g0 = (2.0 - a1c) if sw else a1c
            beta = float(r["beta_std"])
            used = False

            m_, s_ = np.nanmean(g0), np.nanstd(g0)
            if s_ > 0:
                g = np.where(np.isnan(g0), m_, g0)
                acc["XCI"] += beta * (g - m_) / s_
                used = True

            gh = g0.copy()
            gh[isM] = gh[isM] / 2.0
            m_, s_ = np.nanmean(gh), np.nanstd(gh)
            if s_ > 0:
                g = np.where(np.isnan(gh), m_, gh)
                acc["HAP"] += beta * (g - m_) / s_

            z = np.full(Nfam, np.nan)
            good = True
            for grp in (isM, isF):
                x = g0[grp]
                mm, ss = np.nanmean(x), np.nanstd(x)
                if not ss > 0:
                    good = False; break
                xf = np.where(np.isnan(x), mm, x)
                z[grp] = (xf - mm) / ss
            if good:
                acc["SEXSTD"] += beta * np.nan_to_num(z)

            if used:
                n_used += 1

        rec = dict(sv_locus=sv, n_total=len(rws), n_used=n_used,
                   unscoreable_mass_frac=(miss / tot) if tot else 1.0)
        for m in MODELS:
            d = acc[m][cols]
            lo, hi = float(d.min()), float(d.max())
            if n_used == 0 or not (hi > lo):
                rec[f"{m}_ok"] = 0
                continue
            out[m][sv] = 2.0 * (d - lo) / (hi - lo)
            rec[f"{m}_ok"] = 1
        qc.append(rec)
    bf.close()
    return out, qc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", type=int, default=400)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()

    say("=" * 100)
    say("STEP 93 — chrX ploidy models (XCI / HAP / SEXSTD)")
    say("=" * 100)

    fam = [l.split()[1] for l in open(XPFX + ".fam")]
    Nfam = len(fam)
    bpv = (Nfam + 3) // 4
    bim = []
    for l in open(XPFX + ".bim"):
        p = l.rstrip("\n").split("\t")
        bim.append((int(p[3]), p[4], p[5]))
    say(f"chrX .bed: {Nfam:,} samples x {len(bim):,} variants")
    par = [p for p, _, _ in bim if p <= PAR1_END or p >= PAR2_START]
    say(f"  PAR variants in the .bim: {len(par)}  ({'none — all non-PAR' if not par else 'PRESENT'})")

    say("\n[1/4] sex inference from chrX heterozygosity")
    male, rate, obs = infer_sex(fam, bim, bpv, Nfam, a.probe)
    validate_sex(fam, male, rate)
    nm, nf = int((male == 1).sum()), int((male == 0).sum())
    say(f"  full .bed: male {nm:,}  female {nf:,}  undetermined {int(np.isnan(male).sum()):,}")
    with open(f"{OUT}/sex_inference.tsv", "w") as f:
        f.write("group\tn\tmedian_het_rate\n")
        for lab, m in (("male", male == 1), ("female", male == 0)):
            f.write(f"{lab}\t{int(m.sum())}\t{np.nanmedian(rate[m]):.6f}\n")

    say("\n[2/4] scoring 79 chrX predictors under 3 models")
    rows_by_sv = {}
    for r in csv.DictReader(open(WX), delimiter="\t"):
        rows_by_sv.setdefault(r["sv_locus"], []).append(r)
    say(f"  {len(rows_by_sv)} SVs, {sum(len(v) for v in rows_by_sv.values()):,} weight rows")

    iids, iidx = S.load_order()
    famidx = {x: i for i, x in enumerate(fam)}
    cols = np.array([famidx[i] for i in iids])
    say(f"  cohort {len(iids):,} mapped into .bed order")

    dos, qc = score_models(rows_by_sv, bim, bpv, Nfam, male, cols)
    for m in MODELS:
        say(f"  {m:<8} scored {len(dos[m])}/{len(rows_by_sv)} SVs")

    sexc = np.array([1.0 if male[c] == 1 else 0.0 for c in cols])
    with open(f"{OUT}/dosage_compare.tsv", "w") as f:
        f.write("sv_locus\tn_used\tunscoreable_mass_frac\tr_XCI_HAP\tr_XCI_SEXSTD\t"
                "corr_sex_XCI\tcorr_sex_HAP\tcorr_sex_SEXSTD\tsd_XCI\tsd_HAP\tsd_SEXSTD\n")
        rr = {m: [] for m in MODELS}
        pair = {"HAP": [], "SEXSTD": []}
        for rec in qc:
            sv = rec["sv_locus"]
            if sv not in dos["XCI"]:
                continue
            x = dos["XCI"][sv]
            row = [sv, str(rec["n_used"]), f"{rec['unscoreable_mass_frac']:.4f}"]
            for m in ("HAP", "SEXSTD"):
                v = float(np.corrcoef(x, dos[m][sv])[0, 1]) if sv in dos[m] else float("nan")
                row.append(f"{v:.4f}"); pair[m].append(v)
            for m in MODELS:
                v = float(np.corrcoef(dos[m][sv], sexc)[0, 1]) if sv in dos[m] else float("nan")
                row.append(f"{v:.4f}"); rr[m].append(abs(v))
            row += [f"{float(np.std(dos[m][sv])):.4f}" if sv in dos[m] else "NA" for m in MODELS]
            f.write("\t".join(row) + "\n")
    say("\n  agreement with the current (XCI) dosage:")
    for m in ("HAP", "SEXSTD"):
        v = np.array(pair[m]); v = v[np.isfinite(v)]
        say(f"    r(XCI, {m:<6}) median {np.median(v):.4f}   min {v.min():.4f}   "
            f"n<0.95: {int((v < 0.95).sum())}/{len(v)}")
    say("\n  |correlation with sex| of the SV dosage:")
    for m in MODELS:
        v = np.array(rr[m]); v = v[np.isfinite(v)]
        say(f"    {m:<8} median {np.median(v):.4f}   max {v.max():.4f}   "
            f"n>0.30: {int((v > 0.30).sum())}/{len(v)}")

    say("\n[3/4] re-fitting the chrX cells at q<0.10 under each model")
    H = list(csv.DictReader(open(f"{RV}/step30_phewide_ext2/disease_matrices_firth/"
                                 "disease_hits_final.tsv"), delimiter="\t"))
    cells = []
    for r in H:
        q = r.get("q_bh_tier", "")
        if r["tier"] == "ext2_chrX" and q not in ("", "NA") and float(q) < 0.10:
            cells.append(r)
    phes = sorted({r["phecode"] for r in cells})
    say(f"  {len(cells)} cells · {len({r['sv_locus'] for r in cells})} SVs · phecodes {phes}")

    X, _covnames = S.load_X(iidx)
    Y = F.load_pheno_cols(phes, iidx)
    nulls = {}
    for ph in phes:
        y = Y[ph]
        present = y >= 0
        yy = np.where(present, y, 0).astype(float)
        p, w, r_, conv = S.fit_null(X, yy, present)
        nulls[ph] = (yy, present, p, w, r_, conv)
        say(f"    null {ph:<12} cases {int(yy[present].sum()):>7,}  n {int(present.sum()):>7,}  conv {conv}")

    res = []
    for c in cells:
        sv, ph = c["sv_locus"], c["phecode"]
        yy, present, p, w, r_, _ = nulls[ph]
        row = dict(sv_locus=sv, phecode=ph, label=c["label"], cases=c["cases"],
                   p_published=c["p_final"], q_published=c["q_bh_tier"])
        for m in MODELS:
            if sv not in dos[m]:
                row[f"p_{m}"] = "NA"; row[f"or_{m}"] = "NA"; continue
            g = dos[m][sv]
            gm = g - g.mean()
            U = float(gm @ r_)
            Xw = X * w[:, None]
            A = X.T @ Xw
            bcol = Xw.T @ gm
            V = float((gm * w) @ gm - bcol @ np.linalg.solve(A, bcol))
            z = U / np.sqrt(V) if V > 0 else float("nan")
            from math import erfc, sqrt
            ps = erfc(abs(z) / sqrt(2)) if np.isfinite(z) else float("nan")
            row[f"p_score_{m}"] = f"{ps:.3e}"
            beta, pf, conv = F.firth(X, g, yy, present)
            row[f"p_{m}"] = f"{pf:.3e}"
            row[f"or_{m}"] = f"{np.exp(beta):.4f}"
            row[f"conv_{m}"] = int(conv)
        res.append(row)

    fields = (["sv_locus", "phecode", "label", "cases", "p_published", "q_published"] +
              [f"{k}_{m}" for m in MODELS for k in ("p_score", "p", "or", "conv")])
    with open(f"{OUT}/cell_refits.tsv", "w") as f:
        wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        wr.writeheader()
        for r in res:
            wr.writerow(r)

    say("\n[4/4] VERDICT")
    say(f"\n{'sv_locus':<40}{'phecode':<11}{'p_XCI':>11}{'p_HAP':>11}{'p_SEXSTD':>11}  verdict")
    nrob = 0
    for r in sorted(res, key=lambda x: float(x.get("p_XCI", "1") or 1)):
        try:
            ps = [float(r[f"p_{m}"]) for m in MODELS]
        except (ValueError, KeyError):
            say(f"{r['sv_locus'][:39]:<40}{r['phecode']:<11}{'':>33}  INCOMPLETE"); continue
        rob = all(p < 5e-5 for p in ps)
        nrob += rob
        say(f"{r['sv_locus'][:39]:<40}{r['phecode']:<11}" +
            "".join(f"{p:>11.2e}" for p in ps) +
            f"  {'MODEL-INDEPENDENT' if rob else 'DEPENDS ON MODEL'}")
    say(f"\n{nrob}/{len(res)} cells hold under all three ploidy models.")
    say("A cell that holds under all three does not depend on the ploidy convention. One that appears")
    say("only under XCI is a coding artefact and must not be reported as a chrX association.")

    with open(f"{OUT}/PLOIDY_VERDICT.txt", "w") as f:
        f.write("\n".join(LOG) + "\n")
    say(f"\nwrote {OUT}/  ({time.time() - t0:.0f}s)")
    say("CHRX_PLOIDY_DONE")


if __name__ == "__main__":
    main()
