#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd, numpy as np, sys, time, hashlib, argparse
from pathlib import Path
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = Path(CENSV_ROOT + "")
SUF  = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
RES  = PROJ/"results"
LEGACY_DIR = RES/"analysis_qtl/final_tables_cv01_maf005"
LEGACY = LEGACY_DIR/"all_genes_final.tsv"
OUT  = RES/"analysis_qtl/gene_tables_cv01_maf005"
DROOT= RES/"imputation/decompose/cosine_nj"
ENROOT = RES/"elastic_net/decompose/cosine_nj"
CHROMS = ["chr%s" % c for c in list(range(1,23))+["X"]]
PEER = {"eqtl":"peer60", "sqtl":"peer15"}
STRATS = ['edge500','edge1000','edge2000','fksize100k','fksize200k','fksize500k',
          'midpoint500','midpoint1000']
SIG = 0.05

def log(m): sys.stdout.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), m)); sys.stdout.flush()
def md5(p): return hashlib.md5(p.read_bytes()).hexdigest()

ap = argparse.ArgumentParser()
ap.add_argument("--stamp", required=True, help="build timestamp, e.g. 2026-08-06")
args = ap.parse_args()

assert OUT.resolve() != LEGACY_DIR.resolve(), "refusing to write into the legacy directory"
LEGACY_MD5 = md5(LEGACY)
log("legacy table md5 (read-only): %s" % LEGACY_MD5)
OUT.mkdir(parents=True, exist_ok=True)

lg = pd.read_csv(LEGACY, sep="\t", dtype=str)
for c in ("qval","cv_r2","svfreq"): lg[c] = pd.to_numeric(lg[c], errors="coerce")
lg["base"] = lg["ENSG"].str.split(".").str[0]
assert not lg["base"].duplicated().any(), "legacy table has duplicate genes"
log("legacy genes: %d" % len(lg))

need = []
for _, r in lg.iterrows():
    mods = [m.strip() for m in str(r["QTL"]).split(",") if m.strip()]
    assert set(mods) <= {"eqtl","sqtl"}, "unexpected QTL modality %r" % r["QTL"]
    for m in mods:
        need.append(dict(base=r["base"], modality=m))
need = pd.DataFrame(need)
want = set(zip(need["base"], need["modality"]))
log("gene x modality pairs: %d" % len(want))

sv_meta = {}
n_meta_files = 0
for strat in STRATS:
    for ch in CHROMS:
        f = ENROOT/("%s_%s" % (strat, SUF))/ch/"sv_metrics.tsv"
        if not f.exists(): continue
        n_meta_files += 1
        t = pd.read_csv(f, sep="\t", dtype=str)
        have = [c for c in ("svtype","svlen","svpos","svend","svfreq") if c in t.columns]
        for _, x in t[["sv_index"]+have].iterrows():
            sv_meta[(strat, ch, int(x["sv_index"]))] = {c: x[c] for c in have}
log("sv_metrics files: %d ; SV metadata records: %d" % (n_meta_files, len(sv_meta)))

sel = {}
for ens in ("best6","best8"):
    for ch in CHROMS:
        f = DROOT/("%s_%s" % (ens, SUF))/ch/("%s_best_strategy_selection.tsv" % ch)
        if not f.exists(): continue
        t = pd.read_csv(f, sep="\t")
        for _, x in t.iterrows():
            sel[(ens, ch, int(x["sv_index"]))] = (float(x["dosage_cv_r2"]), str(x["dosage_strategy"]))
log("imputation selection records: %d" % len(sel))

def gene_base(pid):
    return (pid.split(":")[-1] if (":" in pid and "ENSG" in pid) else pid).split(".")[0]
def sv_idx(vid):
    p = str(vid).split("_")
    return int(p[1]) if len(p) >= 2 and p[1].isdigit() else -1

rec, nfile = {}, 0
for mod in ("eqtl","sqtl"):
    for ens in ("best6","best8"):
        root = RES/("%s_overlap" % mod)/"decompose"/"cosine_nj_cv01_maf005"/("%s_%s" % (ens, SUF))
        for ch in CHROMS:
            d = root/ch/("sexTrue_gpc5_%s" % PEER[mod])
            for w, wl in (("10000000","w10M"), ("0","w0")):
                f = d/("window_%s" % w)/("%s_cis_permutation.tsv" % ch)
                if not f.exists(): continue
                nfile += 1
                t = pd.read_csv(f, sep="\t")
                t["base"] = t["phenotype_id"].map(gene_base)
                for _, x in t[[(b, mod) in want for b in t["base"]]].iterrows():
                    rec.setdefault((x["base"], mod), []).append(dict(
                        scope="cis", ens=ens, window=wl, sv_chrom=ch,
                        variant_id=x["variant_id"], sv_index=sv_idx(x["variant_id"]),
                        phenotype_id=x["phenotype_id"],
                        stat_name="q_storey", stat=float(x["qval"]),
                        pval=float(x["pval_beta"]), qval_bh_screened=float("nan"),
                        slope=float(x["slope"]), af=float(x["af"])))
            f = d/"trans"/("%s_trans_results.tsv" % ch)
            if not f.exists(): continue
            nfile += 1
            t = pd.read_csv(f, sep="\t")
            t["base"] = t["phenotype_id"].map(gene_base)
            for _, x in t[[(b, mod) in want for b in t["base"]]].iterrows():
                rec.setdefault((x["base"], mod), []).append(dict(
                    scope="trans", ens=ens, window="trans", sv_chrom=ch,
                    variant_id=x["variant_id"], sv_index=sv_idx(x["variant_id"]),
                    phenotype_id=x["phenotype_id"],
                    stat_name="p_nominal", stat=float(x["pval"]),
                    pval=float(x["pval"]), qval_bh_screened=float(x["qval_bh"]),
                    slope=float(x["b"]), af=float(x["af"])))
log("upstream tables read: %d ; pairs with >=1 record: %d / %d" % (nfile, len(rec), len(want)))
missing = want - set(rec)
assert not missing, "gene x modality with no upstream record: %s" % sorted(missing)[:5]

for rs in rec.values():
    for x in rs:
        cv, st = sel.get((x["ens"], x["sv_chrom"], x["sv_index"]), (np.nan, ""))
        x["cv_r2"] = cv; x["strategy"] = st
        m = sv_meta.get((st, x["sv_chrom"], x["sv_index"]), {}) if st else {}
        for k in ("svtype","svlen","svpos","svend","svfreq"): x[k] = m.get(k, "")

LEGACY_COLS = ["ENSG","final_symbol","gene_chr","dist_Mb","Wing","final_biotype","desc"]
rows = []
for _, g in lg.iterrows():
    for mod in [m.strip() for m in str(g["QTL"]).split(",") if m.strip()]:
        pool = rec[(g["base"], mod)]
        pool = [x for x in pool
                if (x["stat"] < 0.10 if x["scope"] == "cis" else x["pval"] < 5e-8)]
        assert pool, "no admitted record for %s %s" % (g["final_symbol"], mod)
        cis = [x for x in pool if x["scope"] == "cis"]
        pool2 = cis if cis else [x for x in pool if x["scope"] == "trans"]
        assert pool2, "no cis or trans record for %s %s" % (g["final_symbol"], mod)
        havecv = [x for x in pool2 if x["cv_r2"] == x["cv_r2"]]
        p = max(havecv, key=lambda x: x["cv_r2"]) if havecv else min(pool2, key=lambda x: x["stat"])
        b = min(pool2, key=lambda x: x["stat"])
        n_sig = sum(1 for x in pool2 if x["stat"] <= SIG)
        picks = [("primary", p)] + ([("bestq", b)] if b["variant_id"] != p["variant_id"] else [])
        for tag, b in picks:
            rows.append(dict(
                ENSG=g["ENSG"], symbol=g["final_symbol"], gene_chr=g["gene_chr"],
                dist_Mb=g["dist_Mb"], wing=g["Wing"], biotype=g["final_biotype"],
                modality=mod, scope=b["scope"], pick=tag,
                n_tested_svs=len(pool2), n_significant_svs=n_sig,
                variant_id=b["variant_id"], sv_chrom=b["sv_chrom"], sv_index=b["sv_index"],
                svtype=b["svtype"], svlen=b["svlen"], svpos=b["svpos"], svend=b["svend"],
                sv_carrier_freq=b["svfreq"], af_in_cohort=b["af"],
                phenotype_id=b["phenotype_id"],
                ens=b["ens"], window=b["window"], strategy=b["strategy"],
                cv_r2_dosage=b["cv_r2"],
                stat_name=b["stat_name"], stat=b["stat"], pval=b["pval"],
                qval_bh_screened=b["qval_bh_screened"], slope=b["slope"],
                significant=int(b["stat"] <= SIG if b["scope"] == "cis"
                                else b["pval"] < 5e-8),
                legacy_sv_id=g["sv_id"], legacy_sv_chrom=g["sv_chrom"],
                legacy_stat=g["qval"], legacy_cv_r2=g["cv_r2"],
                legacy_stat_currency=("p_nominal_MISLABELLED_AS_QVAL"
                                      if "trans" in str(g["cis_trans"]) else "q_storey"),
                sv_differs_from_legacy=int(b["variant_id"] != g["sv_id"]
                                           or b["sv_chrom"] != g["sv_chrom"]),
                desc=g["desc"]))
lo = pd.DataFrame(rows)
lo.to_csv(OUT/"genes_by_modality.tsv", sep="\t", index=False)
log("genes_by_modality.tsv: %d rows (primary %d, bestq %d)"
    % (len(lo), int((lo.pick=="primary").sum()), int((lo.pick=="bestq").sum())))

assert set(lo["ENSG"]) == set(lg["ENSG"]), "GENE SET DRIFTED -- this must never happen"
log("gene set identical to legacy: %d genes  [OK]" % lo["ENSG"].nunique())

one = lg[~lg["QTL"].str.contains(",")]
p1 = lo[lo["pick"]=="primary"].set_index("ENSG")
mism, drift = [], []
for _, g in one.iterrows():
    r = p1.loc[g["ENSG"]]
    why = []
    if r["variant_id"] != g["sv_id"] or r["sv_chrom"] != g["sv_chrom"]:
        why.append("variant")
    ours = r["stat"] if r["scope"] == "cis" else r["pval"]
    drift.append(abs(ours - g["qval"])/g["qval"] if g["qval"] > 0 else 0.0)
    if not np.isclose(ours, g["qval"], rtol=1e-4, atol=0):
        why.append("stat %.6e vs %.6e" % (g["qval"], ours))
    if (ours <= SIG) != (g["qval"] <= SIG):
        why.append("SIGNIFICANCE CLASS FLIPS at %.2f" % SIG)
    if why:
        mism.append((g["final_symbol"], g["QTL"], r["scope"], g["sv_id"],
                     r["variant_id"], "; ".join(why)))
log("")
log("ISOLATION TEST -- single-modality genes must reproduce the legacy record EXACTLY")
log("  (variant AND statistic; for trans the legacy p is compared against our pval)")
log("  single-modality genes: %d ; mismatches: %d" % (len(one), len(mism)))
log("  legacy-vs-upstream relative drift: median %.2e  max %.2e  (writer-chain rounding,"
    % (np.median(drift), np.max(drift)))
log("    not introduced here; our value is read straight from the permutation TSV)")
for m in mism[:20]:
    log("   MISMATCH %-14s %-6s %-6s legacy %-14s -> %-14s  %s" % m)
assert not mism, ("%d single-modality genes disagree with legacy -- v2c is changing more "
                  "than D1/D2" % len(mism))
log("  [OK] D1/D2 correction is isolated; variant reproduces, statistic agrees to 1e-4,")
log("       and no single-modality gene crosses the %.2f boundary." % SIG)

w = lg[LEGACY_COLS+["QTL","cis_trans","sv_id","qval","cv_r2"]].rename(
        columns={"final_symbol":"symbol","sv_id":"legacy_sv_id",
                 "qval":"legacy_stat","cv_r2":"legacy_cv_r2"})
for mod in ("eqtl","sqtl"):
    a = lo[(lo["modality"]==mod) & (lo["pick"]=="primary")].set_index("ENSG")
    i = lo[(lo["modality"]==mod) & (lo["pick"]=="bestq")].set_index("ENSG")
    for c in ("scope","variant_id","sv_chrom","sv_index","phenotype_id","stat_name",
              "stat","pval","cv_r2_dosage","strategy","ens","window","slope","af_in_cohort"):
        w["%s_%s" % (mod, c)] = w["ENSG"].map(a[c])
    for c in ("variant_id","stat","cv_r2_dosage"):
        w["%s_bestq_%s" % (mod, c)] = w["ENSG"].map(i[c])
w.to_csv(OUT/"genes_wide.tsv", sep="\t", index=False)
log("genes_wide.tsv: %d rows x %d cols" % (len(w), len(w.columns)))

d = lo[lo["pick"]=="primary"].copy()
d["stat_ratio"] = d["stat"]/d["legacy_stat"]
d[["ENSG","symbol","modality","scope","pick","legacy_sv_id","legacy_stat",
   "legacy_stat_currency","variant_id","stat_name","stat","stat_ratio",
   "legacy_cv_r2","cv_r2_dosage","sv_differs_from_legacy"]].sort_values(
   "stat_ratio", ascending=False).to_csv(OUT/"DIFF_vs_legacy.tsv", sep="\t", index=False)
log("")
log("---- what the correction changes ----")
log("D1  genes typed eqtl,sqtl now carry a separate record per modality : %d genes"
    % int(lg["QTL"].str.contains(",").sum()))
log("D1  rows whose variant differs from the legacy sv_id               : %d / %d"
    % (int(d["sv_differs_from_legacy"].sum()), len(d)))
log("D2  rows whose legacy 'qval' was actually a nominal p              : %d"
    % int((d["legacy_stat_currency"] != "q_storey").sum()))
tr = d[d["stat_name"]=="p_nominal"]
if len(tr):
    log("D2  those rows now carry stat_name='p_nominal' with the SAME value; the fix is")
    log("    naming, not substitution. qval_bh_screened is carried alongside but is BH")
    log("    over a p<1e-5 pre-screened subset (qtl_overlap.py:524), NOT a genome-wide")
    log("    FDR -- it would be anti-conservative to report it as a q.")
    log("    trans rows admitted at p<5e-8: %d / %d"
        % (int((tr["pval"] < 5e-8).sum()), len(tr)))
log("    genes whose significance CLASS changes at 0.05                 : %d"
    % int(((d["legacy_stat"] <= SIG) != (d["stat"] <= SIG)).sum()))
log("alt  gene x modality where the min-q variant differs from primary  : %d"
    % int((lo["pick"]=="bestq").sum()))
log("     median CV R2 -- primary %.2f  vs  bestq %.2f  (primary is chosen ON cv_r2)"
    % (lo[lo.pick=="primary"]["cv_r2_dosage"].median(),
       lo[lo.pick=="bestq"]["cv_r2_dosage"].median()))
same = lo[lo.pick=="primary"].merge(
    lo[lo.pick=="bestq"][["ENSG","modality","stat"]], on=["ENSG","modality"],
    suffixes=("","_bestq"), how="inner")
if len(same):
    log("     outcome-selection inflation avoided: median stat ratio primary/bestq = %.1fx"
        % (same["stat"]/same["stat_bestq"]).median())

(OUT/"PROVENANCE.txt").write_text(
    "v2c corrected gene table\nbuilt: %s\nbuilder: workflow/scripts/postqtl/gene_table/build_gene_table.py\n\n"
    "LEGACY TABLE READ-ONLY, NOT MODIFIED\n  path: %s\n  md5 : %s\n\n"
    "inputs\n  upstream permutation/trans tables : %d files under {eqtl,sqtl}_overlap/"
    "decompose/cosine_nj_cv01_maf005/\n  imputation selection records      : %d\n"
    "  sv_metrics files                  : %d\n\noutputs\n  genes_by_modality.tsv : %d rows\n"
    "  genes_wide.tsv        : %d rows\n  genes                 : %d (identical to legacy)\n"
    % (args.stamp, LEGACY, LEGACY_MD5, nfile, len(sel), n_meta_files, len(lo), len(w),
       lo["ENSG"].nunique()), encoding="utf-8")

assert md5(LEGACY) == LEGACY_MD5, "LEGACY TABLE CHANGED DURING THE BUILD"
log("")
log("legacy table md5 re-verified unchanged: %s  [OK]" % LEGACY_MD5)
log("DONE")
