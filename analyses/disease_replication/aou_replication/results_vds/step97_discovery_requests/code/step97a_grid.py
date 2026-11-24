#!/usr/bin/env python3
"""step97a_grid.py -- LANE M (matrices)."""
import os, sys, gc, csv, math, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from scipy.stats import chi2 as chi2dist

ROOT = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
        "disease_replication/aou_replication")
ARM  = f"{ROOT}/results_vds/step30_phewide_ext2"
SCORE_ALL = f"{ARM}/disease_matrices"
SCORE_ANC = f"{ARM}/disease_matrices/per_ancestry"
FIRTH_ALL = f"{ARM}/disease_matrices_firth"
FIRTH_ANC = f"{ARM}/disease_matrices_firth/per_ancestry"
CVR2      = f"{ROOT}/results_vds/step28_phewide_ext2/sv_cv_r2_all.tsv"
OUT       = f"{ROOT}/results_vds/step97_discovery_requests"
ANCS = ["EUR", "AFR", "AMR", "EAS", "MID", "SAS"]
TARGET_SV = "chr8_sv30_DEL_mon_best8_midpoint1000"
SUPPRESS = 20

os.makedirs(f"{OUT}/logs", exist_ok=True)
LOG = open(f"{OUT}/logs/step97a_grid.log", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n"); LOG.flush()


def fmt(x, spec="%.6g"):
    if x is None:
        return "NA"
    x = float(x)
    return "NA" if not np.isfinite(x) else spec % x


def scount(c, n):
    if c is None or not np.isfinite(c):
        return "NA", "NA"
    c = int(round(c))
    nn = None if (n is None or not np.isfinite(n)) else int(round(n))
    ctl = None if nn is None else nn - c
    if c < SUPPRESS:
        return "<20", "SUPPRESSED"
    if ctl is not None and ctl < SUPPRESS:
        return str(c), "<20"
    return str(c), ("NA" if ctl is None else str(ctl))


def load_svorder():
    sv, tier = [], []
    with open(f"{SCORE_ALL}/disease_svorder.txt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            sv.append(p[0]); tier.append(p[1] if len(p) > 1 else "")
    return np.array(sv), np.array(tier)


def load_meta(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append({"phecode": r["phecode"], "label": r["label"],
                         "cases": float(r["cases"]), "n": float(r["n"]),
                         "tested": int(r["tested"])})
    return rows


def load_hits(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def bh_from_sorted(p_sorted, targets):
    n = p_sorted.size
    tg = np.asarray(targets, dtype=np.float32)
    if n == 0:
        return np.full(tg.shape, np.nan)
    safe = np.where(np.isfinite(tg), tg, np.float32(1.0))
    rank = np.clip(np.searchsorted(p_sorted, safe, side="left"), 0, n - 1)
    k = np.arange(1, n + 1, dtype=np.float32)
    p_sorted *= np.float32(n)
    p_sorted /= k
    del k
    np.minimum.accumulate(p_sorted[::-1], out=p_sorted[::-1])
    np.minimum(p_sorted, np.float32(1.0), out=p_sorted)
    q = p_sorted[rank].astype(np.float64)
    q[~np.isfinite(tg)] = np.nan
    return q


def collect_tier_p(pv_path, tier_rows, tested_mask, chunk=256):
    mm = np.load(pv_path, mmap_mode="r")
    total = 0
    for s in range(0, len(tier_rows), chunk):
        blk = np.asarray(mm[tier_rows[s:s + chunk]])
        total += int((np.isfinite(blk) & tested_mask[None, :]).sum())
        del blk
    out = np.empty(total, dtype=np.float32)
    o = 0
    for s in range(0, len(tier_rows), chunk):
        blk = np.asarray(mm[tier_rows[s:s + chunk]])
        v = blk[np.isfinite(blk) & tested_mask[None, :]]
        out[o:o + v.size] = v; o += v.size
        del blk, v
    del mm; gc.collect()
    assert o == total, (o, total)
    out.sort()
    return out


t0 = time.time()
say("=" * 78)
say("LANE M  step97a_grid.py   arm =", ARM)
say("=" * 78)

svloc, svtier = load_svorder()
svidx = {s: i for i, s in enumerate(svloc)}
say("svorder:", len(svloc), "rows; tiers =",
    {t: int((svtier == t).sum()) for t in sorted(set(svtier))})

meta_all = load_meta(f"{SCORE_ALL}/disease_meta.tsv")
phecodes = [r["phecode"] for r in meta_all]
phidx = {p: j for j, p in enumerate(phecodes)}
J = len(phecodes)
labels = {r["phecode"]: r["label"] for r in meta_all}
say(f"meta: {J} phecodes; tested(all-anc) = {sum(r['tested'] for r in meta_all)}")

meta_anc = {}
for a in ANCS:
    m = load_meta(f"{SCORE_ANC}/disease_meta_{a}.tsv")
    assert [r["phecode"] for r in m] == phecodes, f"phecode order differs for {a}"
    meta_anc[a] = m
say("per-ancestry meta phecode order identical to all-ancestry for all 6 strata: OK")

cvr2 = {}
with open(CVR2) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        cvr2[r["sv_locus"]] = float(r["cv_r2"])

hits_all = load_hits(f"{FIRTH_ALL}/disease_hits_final.tsv")
tgt = []
for r in hits_all:
    q = float(r["q_bh_tier"]) if r["q_bh_tier"] not in ("", "nan") else np.nan
    if np.isfinite(q) and q < 0.10:
        tgt.append(r)
n_sig = sum(1 for r in tgt if float(r["q_bh_tier"]) < 0.05)
say(f"target cells: {len(tgt)}  (q<0.05: {n_sig}, 0.05<=q<0.10: {len(tgt) - n_sig}) "
    f"| distinct SV {len({r['sv_locus'] for r in tgt})} "
    f"| distinct phecode {len({r['phecode'] for r in tgt})}")
assert len(tgt) == 255 and len({r['sv_locus'] for r in tgt}) == 123 and \
       len({r['phecode'] for r in tgt}) == 55, "255/123/55 contract broken"

tgt_rows = sorted({svidx[r["sv_locus"]] for r in tgt})
tgt_cols = sorted({phidx[r["phecode"]] for r in tgt})
row_pos = {r: i for i, r in enumerate(tgt_rows)}
col_pos = {c: j for j, c in enumerate(tgt_cols)}

exported, psource, psource_full = {}, {}, {}
for a in ANCS:
    hh = load_hits(f"{FIRTH_ANC}/disease_hits_final_{a}.tsv")
    for r in hh:
        exported[(a, r["sv_locus"], r["phecode"])] = 1
        psource[(a, r["sv_locus"], r["phecode"])] = r["p_source"]
    for r in load_hits(f"{FIRTH_ANC}/disease_firth_cells_{a}.tsv"):
        psource_full[(a, r["sv_locus"], r["phecode"])] = r["p_source"]
    say(f"  shipped disease_hits_final_{a}.tsv: {len(hh)} rows; "
        f"disease_firth_cells_{a}.tsv: {sum(1 for k in psource_full if k[0] == a)} refit cells")

say("")
say("--- deliverable 1: grid_255cells_6strata.tsv ---")
grid = {}
dev_max_global = 0.0
beta_dev_global = 0.0
dev_report = []
for a in ANCS:
    def sub(path):
        return np.asarray(np.load(path, mmap_mode="r")[tgt_rows])[:, tgt_cols].astype(np.float64)
    U  = sub(f"{SCORE_ANC}/disease_U_{a}.npy")
    V  = sub(f"{SCORE_ANC}/disease_V_{a}.npy")
    Be = sub(f"{SCORE_ANC}/disease_Be_{a}.npy")
    Pv = sub(f"{SCORE_ANC}/disease_Pv_{a}.npy")
    Bf = sub(f"{FIRTH_ANC}/disease_Be_final_{a}.npy")
    Pf = sub(f"{FIRTH_ANC}/disease_Pv_final_{a}.npy")
    with np.errstate(divide="ignore", invalid="ignore"):
        SE = 1.0 / np.sqrt(V)
        Z = Be / SE
    ok = np.isfinite(Z) & np.isfinite(Pv) & (Pv > 0) & (Pv < 1)
    if ok.any():
        implied = chi2dist.isf(Pv[ok], 1)
        rel = np.abs(Z[ok] ** 2 - implied) / np.maximum(implied, 1e-12)
        dev_report.append((a, int(ok.sum()), float(rel.max()), float(np.median(rel))))
        dev_max_global = max(dev_max_global, float(rel.max()))
    ok2 = np.isfinite(Be) & np.isfinite(U) & np.isfinite(V) & (V != 0)
    if ok2.any():
        bd = np.abs(Be[ok2] - U[ok2] / V[ok2]) / np.maximum(np.abs(Be[ok2]), 1e-30)
        beta_dev_global = max(beta_dev_global, float(bd.max()))
    grid[a] = (Be, SE, Z, Pv, Bf, Pf)
    say(f"  {a}: 6 matrices subset to {len(tgt_rows)}x{len(tgt_cols)}; "
        f"finite score p = {int(np.isfinite(Pv).sum())} / {Pv.size}")
    del U, V, ok, ok2; gc.collect()

say("  internal check   z_score^2  vs  chi2.isf(p_score,1)   [relative deviation]")
for a, n, mx, md in dev_report:
    say(f"    {a}: n={n:6d}  max_rel_dev={mx:.3e}  median_rel_dev={md:.3e}")
say(f"    GLOBAL max relative deviation = {dev_max_global:.3e}")
say(f"    check  beta_score == U/V : GLOBAL max relative deviation = {beta_dev_global:.3e}")

hdr = ["sv_locus", "tier", "phecode", "label", "stratum", "cases", "controls", "n",
       "beta_score_logOR", "se_score", "z_score", "p_score",
       "beta_firth_logOR", "p_firth", "firth_source", "exported_in_p1e4_release",
       "firth_source_full"]
tgt_sorted = sorted(tgt, key=lambda r: float(r["p_final"]))
nrow = 0
with open(f"{OUT}/grid_255cells_6strata.tsv", "w") as f:
    f.write("\t".join(hdr) + "\n")
    for r in tgt_sorted:
        i = svidx[r["sv_locus"]]; j = phidx[r["phecode"]]
        ri, cj = row_pos[i], col_pos[j]
        for a in ANCS:
            Be, SE, Z, Pv, Bf, Pf = grid[a]
            m = meta_anc[a][j]
            cs, ct = scount(m["cases"], m["n"])
            f.write("\t".join([
                r["sv_locus"], svtier[i], r["phecode"], labels.get(r["phecode"], ""), a,
                cs, ct, "%d" % int(m["n"]),
                fmt(Be[ri, cj]), fmt(SE[ri, cj]), fmt(Z[ri, cj]), fmt(Pv[ri, cj], "%.4g"),
                fmt(Bf[ri, cj]), fmt(Pf[ri, cj], "%.4g"),
                psource.get((a, r["sv_locus"], r["phecode"]), ""),
                str(exported.get((a, r["sv_locus"], r["phecode"]), 0)),
                psource_full.get((a, r["sv_locus"], r["phecode"]), ""),
            ]) + "\n")
            nrow += 1
n_exp = sum(exported.get((a, r["sv_locus"], r["phecode"]), 0) for r in tgt for a in ANCS)
say(f"  wrote grid_255cells_6strata.tsv: {nrow} rows "
    f"({len(tgt_sorted)} cells x {len(ANCS)} strata)")
say(f"  {n_exp} of {nrow} grid rows are present in a shipped per-ancestry release; "
    f"{nrow - n_exp} were never exported")
fin_by_anc = {}
for a in ANCS:
    Pv = grid[a][3]
    c = 0
    for r in tgt:
        c += int(np.isfinite(Pv[row_pos[svidx[r["sv_locus"]]], col_pos[phidx[r["phecode"]]]]))
    fin_by_anc[a] = c
say(f"  target cells with a finite score p, by stratum: {fin_by_anc} (of 255 each)")
n_fs = sum(1 for r in tgt for a in ANCS if psource.get((a, r["sv_locus"], r["phecode"]), ""))
n_ff = sum(1 for r in tgt for a in ANCS if psource_full.get((a, r["sv_locus"], r["phecode"]), ""))
say(f"  firth_source populated (from the shipped hits table, as specified): {n_fs}")
say(f"  firth_source_full populated (from disease_firth_cells_<ANC>.tsv): {n_ff}"
    f"  -> {n_ff - n_fs} grid rows ARE Firth-refit but carry a BLANK firth_source")
del grid; gc.collect()

say("")
say("--- deliverable 2: chr8_sv30_uncensored_allanc.tsv ---")
i8 = svidx[TARGET_SV]
tier8 = str(svtier[i8])
say(f"  {TARGET_SV}: matrix row {i8}, tier '{tier8}', cv_r2 = {cvr2.get(TARGET_SV)}")

tested_all = np.array([r["tested"] == 1 for r in meta_all])
Be8 = np.asarray(np.load(f"{SCORE_ALL}/disease_Be.npy", mmap_mode="r")[i8]).astype(np.float64)
Pv8 = np.asarray(np.load(f"{SCORE_ALL}/disease_Pv.npy", mmap_mode="r")[i8]).astype(np.float64)
Bf8 = np.asarray(np.load(f"{FIRTH_ALL}/disease_Be_final.npy", mmap_mode="r")[i8]).astype(np.float64)
Pf8 = np.asarray(np.load(f"{FIRTH_ALL}/disease_Pv_final.npy", mmap_mode="r")[i8]).astype(np.float64)

with np.errstate(divide="ignore", invalid="ignore"):
    pp = np.where(np.isfinite(Pv8) & (Pv8 > 0) & (Pv8 < 1), Pv8, np.nan)
    chi8 = chi2dist.isf(pp, 1)
    Z8 = np.sign(Be8) * np.sqrt(chi8)
    SE8 = np.abs(Be8) / np.abs(Z8)

say("  validating that back-derivation against the stored 1/sqrt(V) per-ancestry:")
for a in ("EUR", "AFR", "AMR"):
    Vv = np.asarray(np.load(f"{SCORE_ANC}/disease_V_{a}.npy", mmap_mode="r")[tgt_rows]).astype(np.float64)
    Bv = np.asarray(np.load(f"{SCORE_ANC}/disease_Be_{a}.npy", mmap_mode="r")[tgt_rows]).astype(np.float64)
    Pp = np.asarray(np.load(f"{SCORE_ANC}/disease_Pv_{a}.npy", mmap_mode="r")[tgt_rows]).astype(np.float64)
    m = np.isfinite(Vv) & np.isfinite(Bv) & np.isfinite(Pp) & (Pp > 0) & (Pp < 1) & (Vv > 0)
    se_true = 1.0 / np.sqrt(Vv[m])
    se_back = np.abs(Bv[m]) / np.sqrt(chi2dist.isf(Pp[m], 1))
    rel = np.abs(se_back - se_true) / se_true
    say(f"    {a}: n={int(m.sum()):8d}  max_rel={rel.max():.3e}  median_rel={np.median(rel):.3e}")
    del Vv, Bv, Pp, m, se_true, se_back, rel; gc.collect()

tier_rows = np.nonzero(svtier == tier8)[0]
say(f"  BH: recomputing q within tier '{tier8}' "
    f"({len(tier_rows)} SVs x {int(tested_all.sum())} tested phecodes)")
p_sorted = collect_tier_p(f"{FIRTH_ALL}/disease_Pv_final.npy", tier_rows, tested_all)
say(f"    valid cells in tier = {p_sorted.size}")
val_t, val_m = [], []
for r in hits_all:
    if r["tier"] != tier8 or r["q_bh_tier"] in ("", "nan"):
        continue
    val_t.append(float(r["p_final"])); val_m.append(float(r["q_bh_tier"]))
if val_t:
    qq = bh_from_sorted(p_sorted.copy(), np.array(val_t, dtype=np.float32))
    d = np.abs(qq - np.array(val_m))
    say(f"    reproduce shipped q_bh_tier on {len(val_t)} '{tier8}' hit rows: "
        f"max|dq| = {np.nanmax(d):.3e}  median|dq| = {np.nanmedian(d):.3e}")
q8 = bh_from_sorted(p_sorted, Pf8.astype(np.float32))
del p_sorted; gc.collect()

cand = np.nonzero(tested_all & np.isfinite(Pf8))[0]
pf_safe = np.where(np.isfinite(Pf8[cand]), Pf8[cand], 1.0)
order = cand[np.argsort(pf_safe, kind="stable")]
best_j = int(order[0]); best_p = float(Pf8[best_j]); best_q = float(q8[best_j])
say(f"  tested phecodes with a finite p for this SV: {len(cand)}")
say(f"  BEST p across all tested phecodes = {best_p:.6g}  "
    f"({phecodes[best_j]} / {labels[phecodes[best_j]]}); its q_bh_tier = {best_q:.4g}")
say(f"  export rule is (q_bh_tier<0.10) OR (p<1e-4) [step53:115]. best p / 1e-4 = "
    f"{best_p / 1e-4:.2f}  -> {'ABOVE the 1e-4 threshold' if best_p >= 1e-4 else 'BELOW the 1e-4 threshold'}")

n_firth8 = 0
with open(f"{FIRTH_ALL}/disease_firth_cells.tsv") as f:
    for line in f:
        if line.startswith(TARGET_SV + "\t"):
            n_firth8 += 1
say(f"  all-ancestry Firth-overlay cells for this SV = {n_firth8}")

with open(f"{OUT}/chr8_sv30_uncensored_allanc.tsv", "w") as f:
    f.write(f"# SV = {TARGET_SV}\n")
    f.write(f"# matrix row {i8} (0-based, disease_svorder.txt); tier = {tier8}; "
            f"cv_r2 = {cvr2.get(TARGET_SV)} (SINGLE GLOBAL elastic-net number, NOT per-stratum)\n")
    f.write(f"# cohort = ALL-ANCESTRY; {len(cand)} phecodes tested with a finite p "
            f"(of {int(tested_all.sum())} tested / {J} in the panel)\n")
    f.write(f"# BEST p across all tested phecodes = {best_p:.6g}  "
            f"({phecodes[best_j]}: {labels[phecodes[best_j]]}), q_bh_tier = {best_q:.6g}\n")
    f.write(f"# The shipped release exports a cell when (q_bh_tier<0.10) OR (p<1e-4) "
            f"[step53_build_final_firth_matrix_ext2.py:115]. This SV's best all-ancestry p "
            f"is {best_p / 1e-4:.2f}x the 1e-4 threshold, i.e. "
            f"{'ABOVE it' if best_p >= 1e-4 else 'BELOW it'}, and its best q_bh_tier is "
            f"{best_q:.4g}. NO all-ancestry cell for this SV appears in disease_hits_final.tsv. "
            f"This table is the uncensored view.\n")
    f.write(f"# beta_firth_logOR / p_firth: this SV has {n_firth8} all-ancestry Firth-refit "
            f"cells, so where that count is 0 these columns ARE the score beta/p, which is the "
            f"base step53 initialises the final matrices from. They are not independent Firth fits.\n")
    f.write("# se_score is BACK-DERIVED as |beta_score|/sqrt(chi2.isf(p_score,1)): no "
            "all-ancestry U/V matrices were ever written. This is exact algebra on the same "
            "(U,V) pair, validated against the stored per-ancestry 1/sqrt(V) (see log).\n")
    f.write("# A non-finite p means 'not tested / no evidence', not missing data.\n")
    f.write("\t".join(["phecode", "label", "cases", "controls", "n",
                       "beta_score_logOR", "se_score", "z_score", "p_score",
                       "beta_firth_logOR", "p_firth", "q_bh_tier_as_computed", "tier"]) + "\n")
    for j in order[:50]:
        j = int(j); m = meta_all[j]
        cs, ct = scount(m["cases"], m["n"])
        f.write("\t".join([
            phecodes[j], labels[phecodes[j]], cs, ct, "%d" % int(m["n"]),
            fmt(Be8[j]), fmt(SE8[j]), fmt(Z8[j]), fmt(Pv8[j], "%.4g"),
            fmt(Bf8[j]), fmt(Pf8[j], "%.4g"), fmt(q8[j], "%.4g"), tier8]) + "\n")
say(f"  wrote chr8_sv30_uncensored_allanc.tsv (top 50 of {len(cand)} tested phecodes)")
top20 = [int(j) for j in order[:20]]

say("")
say("--- deliverable 3: chr8_sv30_perancestry.tsv + chr8_sv30_heterogeneity.tsv ---")
per = {}
for a in ANCS:
    Ua = np.asarray(np.load(f"{SCORE_ANC}/disease_U_{a}.npy", mmap_mode="r")[i8]).astype(np.float64)
    Va = np.asarray(np.load(f"{SCORE_ANC}/disease_V_{a}.npy", mmap_mode="r")[i8]).astype(np.float64)
    Ba = np.asarray(np.load(f"{SCORE_ANC}/disease_Be_{a}.npy", mmap_mode="r")[i8]).astype(np.float64)
    Pa = np.asarray(np.load(f"{SCORE_ANC}/disease_Pv_{a}.npy", mmap_mode="r")[i8]).astype(np.float64)
    Bfa = np.asarray(np.load(f"{FIRTH_ANC}/disease_Be_final_{a}.npy", mmap_mode="r")[i8]).astype(np.float64)
    Pfa = np.asarray(np.load(f"{FIRTH_ANC}/disease_Pv_final_{a}.npy", mmap_mode="r")[i8]).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        SEa = 1.0 / np.sqrt(Va)
        bmax = np.nanmax(np.abs(Ba - Ua / Va)) if np.isfinite(Va).any() else np.nan
    per[a] = (Ba, SEa, Pa, Bfa, Pfa)
    say(f"  {a}: finite score p for this SV = {int(np.isfinite(Pa).sum())} / {J}; "
        f"max|Be-U/V| = {bmax:.3e}")
    del Ua, Va; gc.collect()

with open(f"{OUT}/chr8_sv30_perancestry.tsv", "w") as f:
    f.write(f"# SV = {TARGET_SV} (matrix row {i8}, tier {tier8}). "
            f"Top 20 phecodes by ALL-ANCESTRY p.\n")
    f.write("# se_score = 1/sqrt(V) read from the stored per-ancestry V matrix (exact here, "
            "not back-derived).\n")
    f.write("# beta_firth_logOR / p_firth come from disease_{Be,Pv}_final_<ANC>.npy, which are "
            "MIXED: exact-Firth where it converged, approx-Firth where not, score p as base. "
            "No SE was ever stored for a Firth fit -- never pair a Firth beta with se_score.\n")
    f.write("# NA = non-finite = 'not tested / no evidence in this stratum', not missing data.\n")
    f.write("\t".join(["phecode", "label", "stratum", "cases", "controls", "n",
                       "beta_score_logOR", "se_score", "p_score",
                       "beta_firth_logOR", "p_firth"]) + "\n")
    for j in top20:
        for a in ANCS:
            Ba, SEa, Pa, Bfa, Pfa = per[a]
            m = meta_anc[a][j]
            cs, ct = scount(m["cases"], m["n"])
            f.write("\t".join([
                phecodes[j], labels[phecodes[j]], a, cs, ct, "%d" % int(m["n"]),
                fmt(Ba[j]), fmt(SEa[j]), fmt(Pa[j], "%.4g"),
                fmt(Bfa[j]), fmt(Pfa[j], "%.4g")]) + "\n")
say(f"  wrote chr8_sv30_perancestry.tsv: {len(top20) * len(ANCS)} rows")

het = []
for j in range(J):
    bs, ss, used = [], [], []
    for a in ANCS:
        Ba, SEa, Pa, _, _ = per[a]
        if not (np.isfinite(Ba[j]) and np.isfinite(SEa[j]) and SEa[j] > 0):
            continue
        if meta_anc[a][j]["cases"] < SUPPRESS:
            continue
        bs.append(Ba[j]); ss.append(SEa[j]); used.append(a)
    k = len(bs)
    if k == 0:
        continue
    b = np.array(bs); s = np.array(ss); w = 1.0 / s ** 2
    sw = w.sum(); bm = float((w * b).sum() / sw); sem = float(1.0 / math.sqrt(sw))
    zm = bm / sem; pm = float(chi2dist.sf(zm * zm, 1))
    Q = float((w * (b - bm) ** 2).sum()); df = k - 1
    if k > 1:
        I2 = float(max(0.0, (Q - df) / Q) * 100.0) if Q > 0 else 0.0
        C = sw - (w ** 2).sum() / sw
        tau2 = float(max(0.0, (Q - df) / C)) if C > 0 else np.nan
    else:
        I2 = np.nan; tau2 = np.nan
    het.append((j, k, used, bm, sem, pm, Q, df, I2, tau2))
het.sort(key=lambda t: t[5])

with open(f"{OUT}/chr8_sv30_heterogeneity.tsv", "w") as f:
    f.write(f"# SV = {TARGET_SV} (matrix row {i8}, tier {tier8}).\n")
    f.write("# Fixed-effect inverse-variance meta-analysis over the 6 ancestry strata using "
            "beta_score and se_score = 1/sqrt(V) ONLY -- the one coherent beta/SE pair. Firth "
            "betas are NOT meta-analysed because no SE was ever stored for a Firth fit.\n")
    f.write("# A stratum contributes only if its cell is finite AND that stratum has >=20 cases "
            "for the phecode. k_strata and strata_used name exactly which contributed.\n")
    f.write("# Q = Cochran's Q, df = k-1, I2_pct = max(0,(Q-df)/Q)*100, tau2 = DerSimonian-Laird. "
            "I2_pct and tau2 are NA when k=1 (heterogeneity is undefined with one stratum).\n")
    f.write("# Sorted by p_meta ascending. All tested phecodes with k>=1 are listed.\n")
    f.write("\t".join(["phecode", "label", "k_strata", "strata_used",
                       "beta_meta", "se_meta", "p_meta", "Q", "df", "I2_pct", "tau2",
                       "p_allanc_score"]) + "\n")
    for (j, k, used, bm, sem, pm, Q, df, I2, tau2) in het:
        f.write("\t".join([phecodes[j], labels[phecodes[j]], str(k), ",".join(used),
                           fmt(bm), fmt(sem), fmt(pm, "%.4g"), fmt(Q), str(df),
                           fmt(I2, "%.4g"), fmt(tau2), fmt(Pv8[j], "%.4g")]) + "\n")
kdist = {}
for t in het:
    kdist[t[1]] = kdist.get(t[1], 0) + 1
say(f"  wrote chr8_sv30_heterogeneity.tsv: {len(het)} phecodes; "
    f"k_strata distribution = {dict(sorted(kdist.items()))}")

say("")
say("file sizes:")
for fn in ("grid_255cells_6strata.tsv", "chr8_sv30_uncensored_allanc.tsv",
           "chr8_sv30_perancestry.tsv", "chr8_sv30_heterogeneity.tsv"):
    say(f"  {fn}: {os.path.getsize(os.path.join(OUT, fn)) / 1e6:.3f} MB")
say("")
say(f"peak RSS = {__import__('resource').getrusage(__import__('resource').RUSAGE_SELF).ru_maxrss / 1024:.0f} MB")
say(f"DONE in {time.time() - t0:.1f}s")
LOG.close()
