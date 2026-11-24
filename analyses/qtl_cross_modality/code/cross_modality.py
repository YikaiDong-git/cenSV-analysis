#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd, numpy as np, sys, time
from pathlib import Path
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = Path(CENSV_ROOT + "")
SUF  = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
RES  = PROJ/"results"
V2C  = RES/"analysis_qtl/gene_tables_cv01_maf005/genes_by_modality.tsv"
OUT  = PROJ/"analyses/qtl_cross_modality"
PEER = {"eqtl":"peer60", "sqtl":"peer15"}
WIN  = {"w0":"0", "w10M":"10000000"}

def log(m): sys.stdout.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), m)); sys.stdout.flush()

d = pd.read_csv(V2C, sep="\t")
p = d[(d["pick"]=="primary") & (d["scope"]=="cis")]
w = p.pivot_table(index="ENSG", columns="modality",
                  values=["variant_id","stat","cv_r2_dosage","phenotype_id","ens","window","sv_chrom"],
                  aggfunc="first").dropna()
sym = p.drop_duplicates("ENSG").set_index("ENSG")["symbol"]
log("dual-modality cis genes: %d" % len(w))

cache = {}
def nominal(mod, ens, window, ch):
    k = (mod, ens, window, ch)
    if k not in cache:
        f = (RES/("%s_overlap" % mod)/"decompose"/"cosine_nj_cv01_maf005"/("%s_%s" % (ens, SUF))
             /ch/("sexTrue_gpc5_%s" % PEER[mod])/("window_%s" % WIN[window])
             /("%s_cis_nominal.cis_qtl_pairs.%s.parquet" % (ch, ch)))
        cache[k] = pd.read_parquet(f) if f.exists() else None
        log("  %s %s %s %s -> %s" % (mod, ens, window, ch,
            "%d pairs" % len(cache[k]) if cache[k] is not None else "MISSING"))
    return cache[k]

def lookup(mod, ens, window, ch, variant, pheno_prefix, base):
    t = nominal(mod, ens, window, ch)
    if t is None: return None
    g = t[t["phenotype_id"].astype(str).str.contains(base, regex=False)]
    if not len(g): return None
    n_var = g["variant_id"].nunique()
    hit = g[g["variant_id"] == variant]
    if not len(hit): return dict(found=0, n_var=n_var, n_pheno=g["phenotype_id"].nunique())
    best = hit.loc[hit["pval_nominal"].idxmin()]
    lead = g["pval_nominal"].min()
    rank = int((g.groupby("variant_id")["pval_nominal"].min() < best["pval_nominal"]).sum()) + 1
    return dict(found=1, p=float(best["pval_nominal"]), slope=float(best["slope"]),
                af=float(best["af"]), n_var=n_var, n_pheno=g["phenotype_id"].nunique(),
                rank=rank, lead_p=float(lead), bonf=min(1.0, float(best["pval_nominal"])*n_var),
                pheno=str(best["phenotype_id"]))

rows = []
for ensg in w.index:
    r = w.loc[ensg]; base = ensg.split(".")[0]
    A, B = r[("variant_id","eqtl")], r[("variant_id","sqtl")]
    rec = dict(ENSG=ensg, symbol=sym.get(ensg, ensg), same_sv=int(A == B),
               eqtl_sv=A, sqtl_sv=B,
               eqtl_lead_q=r[("stat","eqtl")], sqtl_lead_q=r[("stat","sqtl")],
               eqtl_cv_r2=r[("cv_r2_dosage","eqtl")], sqtl_cv_r2=r[("cv_r2_dosage","sqtl")])
    x = lookup("sqtl", r[("ens","sqtl")], r[("window","sqtl")], "chr"+str(r[("sv_chrom","sqtl")]).replace("chr",""),
               A, None, base)
    if x and x.get("found"):
        rec.update(eSV_on_splice_p=x["p"], eSV_on_splice_bonf=x["bonf"],
                   eSV_on_splice_rank=x["rank"], splice_n_var=x["n_var"],
                   splice_lead_p=x["lead_p"])
    else:
        rec.update(eSV_on_splice_p=np.nan, eSV_on_splice_bonf=np.nan,
                   eSV_on_splice_rank=np.nan, splice_n_var=(x or {}).get("n_var", np.nan),
                   splice_lead_p=np.nan)
    y = lookup("eqtl", r[("ens","eqtl")], r[("window","eqtl")], "chr"+str(r[("sv_chrom","eqtl")]).replace("chr",""),
               B, None, base)
    if y and y.get("found"):
        rec.update(sSV_on_expr_p=y["p"], sSV_on_expr_bonf=y["bonf"],
                   sSV_on_expr_rank=y["rank"], expr_n_var=y["n_var"],
                   expr_lead_p=y["lead_p"])
    else:
        rec.update(sSV_on_expr_p=np.nan, sSV_on_expr_bonf=np.nan,
                   sSV_on_expr_rank=np.nan, expr_n_var=(y or {}).get("n_var", np.nan),
                   expr_lead_p=np.nan)
    rows.append(rec)

o = pd.DataFrame(rows)
OUT.mkdir(parents=True, exist_ok=True)
(OUT/"output").mkdir(exist_ok=True)
o.to_csv(OUT/"output"/"cross_modality_lead_sv.tsv", sep="\t", index=False)

diff = o[o["same_sv"] == 0]
log("")
log("genes with DIFFERENT lead SVs: %d" % len(diff))
log("")
log("Does the eQTL's lead SV also affect SPLICING?")
log("  tested: %d | Bonferroni p<0.05: %d | nominal p<0.05: %d"
    % (int(diff["eSV_on_splice_p"].notna().sum()),
       int((diff["eSV_on_splice_bonf"] < 0.05).sum()),
       int((diff["eSV_on_splice_p"] < 0.05).sum())))
log("Does the sQTL's lead SV also affect EXPRESSION?")
log("  tested: %d | Bonferroni p<0.05: %d | nominal p<0.05: %d"
    % (int(diff["sSV_on_expr_p"].notna().sum()),
       int((diff["sSV_on_expr_bonf"] < 0.05).sum()),
       int((diff["sSV_on_expr_p"] < 0.05).sum())))
log("")
c = diff[["symbol","eqtl_sv","sqtl_sv","eqtl_lead_q","sqtl_lead_q",
          "eSV_on_splice_p","eSV_on_splice_bonf","sSV_on_expr_p","sSV_on_expr_bonf"]]
log(c.to_string(index=False, float_format=lambda v: "%.2e" % v))
log("")
log("DONE")
