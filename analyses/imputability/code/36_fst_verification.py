#!/usr/bin/env python3
from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

PROJECT = Path(CENSV_ROOT + "")
OUT_DIR = PROJECT / "analyses/imputability/output"
AUDIT_DIR = OUT_DIR / "fst"
AUDIT_DIR.mkdir(parents=True, exist_ok=True)
SUFFIX = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
POP_COLORS = {"AFR": "#e41a1c", "AMR": "#377eb8", "EAS": "#4daf4a",
              "EUR": "#984ea3", "SAS": "#ff7f00"}

def atomic_write(df, path):
    tmp = path.parent / f"{path.name}.tmp.{os.getpid()}"
    df.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)


print("[1/8] Panel composition...", flush=True)
panel = pd.read_csv(PROJECT / "resources/1000G/20130606_g1k_3202_samples_ped_population.txt", sep=r"\s+", engine="python")
panel_pop_map = dict(zip(panel["SampleID"], panel["Superpopulation"]))
print(f"  1KGP panel: {len(panel_pop_map)} samples")
print(f"  super-pop breakdown: {panel['Superpopulation'].value_counts().to_dict()}")

sv_samples_chr1 = (PROJECT / f"sv_extraction/decompose/{SUFFIX}/chr1/chr1_sv_samples.txt").read_text().splitlines()
matched = [s for s in sv_samples_chr1 if s in panel_pop_map]
unmatched = [s for s in sv_samples_chr1 if s not in panel_pop_map]
sv_panel_pops = pd.Series([panel_pop_map[s] for s in matched]).value_counts()
print(f"  HPRC SV samples (chr1): {len(sv_samples_chr1)}")
print(f"  matched to 1KGP: {len(matched)}  unmatched: {len(unmatched)}")
print(f"  super-pop in HPRC matched: {sv_panel_pops.to_dict()}")

panel_stats = pd.DataFrame({
    "super_pop": POPS,
    "n_1KGP_panel": [int((panel["Superpopulation"]==p).sum()) for p in POPS],
    "n_HPRC_matched": [int(sv_panel_pops.get(p, 0)) for p in POPS],
})
panel_stats.loc[len(panel_stats)] = {
    "super_pop": "TOTAL",
    "n_1KGP_panel": int(len(panel)),
    "n_HPRC_matched": int(len(matched)),
}
atomic_write(panel_stats, AUDIT_DIR / "panel_stats.tsv")

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
ax = axes[0]
y = np.arange(len(POPS))
ax.barh(y - 0.2, panel_stats.iloc[:5]["n_1KGP_panel"], 0.4, color=[POP_COLORS[p] for p in POPS], alpha=0.4, label="1KGP panel (3202)")
ax.barh(y + 0.2, panel_stats.iloc[:5]["n_HPRC_matched"], 0.4, color=[POP_COLORS[p] for p in POPS], edgecolor="black", label="HPRC matched (188)")
ax.set_yticks(y); ax.set_yticklabels(POPS)
ax.set_xlabel("sample count")
ax.set_title(f"Panel composition\n(1KGP n={len(panel)}; HPRC matched n={len(matched)})")
ax.legend()
ax.grid(alpha=0.3, axis="x")

ax = axes[1]
sizes = [int(sv_panel_pops.get(p, 0)) for p in POPS]
ax.pie(sizes, labels=[f"{p}\nn={s}" for p, s in zip(POPS, sizes)],
       colors=[POP_COLORS[p] for p in POPS], startangle=90, autopct="%.1f%%")
ax.set_title("HPRC SV-panel super-pop fraction (n=188)")
plt.tight_layout()
fig.savefig(AUDIT_DIR / "panel_composition.jpg", format="jpg", dpi=200, bbox_inches="tight", pil_kwargs={"quality": 92})
plt.close(fig)

print("[2/8] Loading data...", flush=True)
anc_full = pd.read_csv(OUT_DIR / "sv_ancestry_features.tsv", sep="\t")
print(f"  full catalog F_ST table: {len(anc_full)} rows")
v3 = pd.read_csv(OUT_DIR / "per_sv_features_with_genome.tsv", sep="\t")
print(f"  imputable subset (v3): {len(v3)} rows")

print("[3/8] F_ST descriptive stats...", flush=True)
def fst_desc(s, label):
    s = s.dropna()
    return {
        "label": label, "n": len(s),
        "mean": s.mean(), "median": s.median(),
        "p25": s.quantile(0.25), "p75": s.quantile(0.75),
        "p95": s.quantile(0.95), "p99": s.quantile(0.99),
        "max": s.max(),
        "n_above_0.15": int((s > 0.15).sum()),
        "frac_above_0.15": (s > 0.15).mean(),
        "n_neg": int((s < 0).sum()),
        "frac_neg": (s < 0).mean(),
    }
desc_rows = [
    fst_desc(anc_full["max_FST_pair_sv"], "catalog (192k)"),
    fst_desc(v3["max_FST_pair_sv"], "imputable (10,705)"),
]
desc_df = pd.DataFrame(desc_rows)
atomic_write(desc_df, AUDIT_DIR / "fst_descriptive_stats.tsv")
print(desc_df.to_string(index=False))

print("[4/8] F_ST histogram...", flush=True)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
bins = np.linspace(-0.1, 0.6, 71)
ax.hist(anc_full["max_FST_pair_sv"].dropna().clip(-0.1, 0.6), bins=bins, color="#999999", alpha=0.5, label="catalog (192k)", density=True)
ax.hist(v3["max_FST_pair_sv"].dropna().clip(-0.1, 0.6), bins=bins, color="#1f77b4", alpha=0.6, label="imputable (10.7k)", density=True)
ax.axvline(0.15, color="red", lw=1, ls="--", label="F_ST=0.15 (inter-continent)")
ax.axvline(0, color="black", lw=0.5, ls="-")
ax.set_xlabel("max_FST_pair_sv (Hudson, sample-size corrected)")
ax.set_ylabel("density")
ax.set_title("F_ST distribution — catalog vs imputable")
ax.legend(fontsize=8)
ax.set_yscale("log")
ax.grid(alpha=0.3)

ax = axes[1]
ax.hist(v3["max_FST_pair_sv"].dropna().clip(0, 0.6), bins=np.linspace(0, 0.6, 61), color="#1f77b4", alpha=0.6)
ax.axvline(0.15, color="red", lw=1, ls="--", label="F_ST=0.15")
n_high = (v3["max_FST_pair_sv"] > 0.15).sum()
ax.set_xlabel("max_FST_pair_sv")
ax.set_ylabel("count")
ax.set_title(f"Imputable-subset detail\n({n_high}/{len(v3)} = {100*n_high/len(v3):.2f}% above 0.15)")
ax.set_yscale("log")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig(AUDIT_DIR / "fst_distribution.jpg", format="jpg", dpi=200, bbox_inches="tight", pil_kwargs={"quality": 92})
plt.close(fig)

print("[5/8] Per-chrom F_ST stats...", flush=True)
def chrom_int(c):
    s = c.replace("chr", "")
    return 23 if s == "X" else int(s)
per_chrom = v3.groupby("chrom").agg(
    n_imputable=("max_FST_pair_sv", "count"),
    mean_FST=("max_FST_pair_sv", "mean"),
    median_FST=("max_FST_pair_sv", "median"),
    n_high=("max_FST_pair_sv", lambda v: (v > 0.15).sum()),
    pct_high=("max_FST_pair_sv", lambda v: 100 * (v > 0.15).mean()),
    mean_cv_r2_best8=("cv_r2_best8", "mean"),
).reset_index()
per_chrom["chrom_int"] = per_chrom["chrom"].apply(chrom_int)
per_chrom = per_chrom.sort_values("chrom_int").drop(columns="chrom_int")
atomic_write(per_chrom, AUDIT_DIR / "fst_per_chrom.tsv")

fig, ax = plt.subplots(figsize=(11, 4.5))
x = np.arange(len(per_chrom))
ax.bar(x, per_chrom["mean_FST"], color="#1f77b4", edgecolor="black", label="mean F_ST")
ax.bar(x, per_chrom["pct_high"] / 100 * per_chrom["mean_FST"].max() * 2, color="#ff7f0e", alpha=0.5, label="% F_ST>0.15 (scaled)")
ax.set_xticks(x)
ax.set_xticklabels(per_chrom["chrom"], rotation=45, fontsize=8)
ax.set_ylabel("mean max_FST_pair_sv")
ax.set_title("Per-chromosome F_ST in imputable subset (n=10,705)\n"
              "Bars: mean F_ST. Orange: % SVs above F_ST=0.15 (scaled to bar height).")
ax.axhline(v3["max_FST_pair_sv"].mean(), color="grey", lw=0.5, ls="--")
ax.legend(fontsize=8)
ax.grid(alpha=0.3, axis="y")
plt.tight_layout()
fig.savefig(AUDIT_DIR / "fst_per_chrom.jpg", format="jpg", dpi=200, bbox_inches="tight", pil_kwargs={"quality": 92})
plt.close(fig)

print("[6/8] F_ST vs CV R² scatter...", flush=True)
fig, ax = plt.subplots(figsize=(7, 5))
sub = v3.dropna(subset=["max_FST_pair_sv", "cv_r2_best8"])
hb = ax.hexbin(sub["max_FST_pair_sv"], sub["cv_r2_best8"], gridsize=40,
               cmap="Blues", mincnt=1, bins="log")
ax.set_xlabel("max_FST_pair_sv")
ax.set_ylabel("cv_r2_best8")
from scipy.stats import spearmanr
r, p = spearmanr(sub["max_FST_pair_sv"], sub["cv_r2_best8"])
ax.set_title(f"F_ST vs CV R² (imputable, n={len(sub)})\nSpearman ρ = {r:+.4f}  p = {p:.3g}")
ax.axvline(0.15, color="red", lw=1, ls="--", alpha=0.5)
fig.colorbar(hb, label="density (log)")
ax.grid(alpha=0.3)
plt.tight_layout()
fig.savefig(AUDIT_DIR / "fst_vs_cvr2.jpg", format="jpg", dpi=200, bbox_inches="tight", pil_kwargs={"quality": 92})
plt.close(fig)

print("[7/8] eGene/sGene hit enrichment...", flush=True)
tiers = ["cv01_maf005", "cv01_maf010", "cv03_maf005", "cv03_maf010"]
enrich_rows = []
all_hit_dfs = {}
for tier in tiers:
    hits = pd.read_csv(PROJECT / f"results/analysis_qtl/final_tables_{tier}/all_genes_final.tsv", sep="\t")
    hits["chrom"] = hits["sv_chrom"]; hits["pos"] = hits["svpos"]
    merged = hits.merge(v3[["chrom", "pos", "svtype", "svlen", "max_FST_pair_sv",
                              "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
                              "cv_r2_best8"]],
                         on=["chrom", "pos", "svtype", "svlen"], how="left")
    all_hit_dfs[tier] = merged
    m = merged.dropna(subset=["max_FST_pair_sv"])
    n_total = len(merged)
    n_match = len(m)
    n_high_hit = int((m["max_FST_pair_sv"] > 0.15).sum())
    frac_high_hit = (m["max_FST_pair_sv"] > 0.15).mean()
    frac_high_cat = (v3["max_FST_pair_sv"] > 0.15).mean()
    n_imp_v3 = len(v3)
    n_imp_high = int((v3["max_FST_pair_sv"] > 0.15).sum())
    a = n_high_hit
    b = n_match - n_high_hit
    c = n_imp_high - n_high_hit
    d = (n_imp_v3 - n_match) - c
    odds, p_fisher = fisher_exact([[a, b], [c, d]], alternative="greater")
    enrich_rows.append({
        "tier": tier, "n_hits": n_total, "n_matched_to_v3": n_match,
        "mean_FST_hits": m["max_FST_pair_sv"].mean(),
        "median_FST_hits": m["max_FST_pair_sv"].median(),
        "n_high_FST_hits": n_high_hit,
        "frac_high_FST_hits_pct": frac_high_hit * 100,
        "catalog_frac_high_pct": frac_high_cat * 100,
        "enrichment_fold": frac_high_hit / frac_high_cat if frac_high_cat > 0 else np.nan,
        "fisher_p_one_sided": p_fisher,
    })
enrich_df = pd.DataFrame(enrich_rows)
atomic_write(enrich_df, AUDIT_DIR / "hit_enrichment.tsv")
print(enrich_df.to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
x = np.arange(len(tiers))
ax.bar(x, enrich_df["frac_high_FST_hits_pct"], color="#d62728", label="hits % at F_ST > 0.15")
ax.bar(x, enrich_df["catalog_frac_high_pct"], color="#999999", alpha=0.7, width=0.4, label="imputable baseline (1.20%)")
ax.set_xticks(x); ax.set_xticklabels(tiers, rotation=20, fontsize=9)
ax.set_ylabel("% SVs with F_ST > 0.15")
ax.set_title("eGene/sGene hit enrichment for high-F_ST SVs")
for i, row in enrich_df.iterrows():
    ax.text(i, row["frac_high_FST_hits_pct"] + 0.2, f"{row['enrichment_fold']:.1f}×",
            ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.legend(fontsize=8)
ax.grid(alpha=0.3, axis="y")

ax = axes[1]
def cdf(s, label, color):
    s = s.dropna().sort_values()
    y = np.linspace(0, 1, len(s))
    ax.plot(s, y, label=label, color=color, lw=2)
cdf(v3["max_FST_pair_sv"], "imputable (10.7k)", "#999999")
for tier, color in [("cv03_maf005", "#d62728"), ("cv03_maf010", "#1f77b4")]:
    cdf(all_hit_dfs[tier]["max_FST_pair_sv"], f"hits ({tier})", color)
ax.axvline(0.15, color="red", lw=1, ls="--", alpha=0.5)
ax.set_xlabel("max_FST_pair_sv")
ax.set_ylabel("cumulative fraction")
ax.set_title("F_ST CDF: catalog vs cTWAS-tier hits")
ax.legend(fontsize=8)
ax.set_xlim(0, 0.4)
ax.grid(alpha=0.3)
plt.tight_layout()
fig.savefig(AUDIT_DIR / "hit_enrichment.jpg", format="jpg", dpi=200, bbox_inches="tight", pil_kwargs={"quality": 92})
plt.close(fig)

print("[8/8] Top high-F_ST tables...", flush=True)
top20_cat = v3.dropna(subset=["max_FST_pair_sv"]).nlargest(20, "max_FST_pair_sv")[
    ["chrom", "pos", "svtype", "svlen", "max_FST_pair_sv",
     "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
     "total_af_sv", "cv_r2_best8", "sv_category"] if "sv_category" in v3.columns else
    ["chrom", "pos", "svtype", "svlen", "max_FST_pair_sv",
     "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
     "total_af_sv", "cv_r2_best8"]]
atomic_write(top20_cat, AUDIT_DIR / "top20_high_fst_catalog.tsv")

top20_hits = all_hit_dfs["cv03_maf005"].dropna(subset=["max_FST_pair_sv"]).nlargest(20, "max_FST_pair_sv")
cols = ["final_symbol", "QTL", "chrom", "pos", "svtype", "svlen",
        "max_FST_pair_sv", "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
        "cv_r2_best8", "qval"]
cols = [c for c in cols if c in top20_hits.columns]
atomic_write(top20_hits[cols], AUDIT_DIR / "top20_high_fst_hits.tsv")

fig, ax = plt.subplots(figsize=(7, 7))
data = top20_cat[["af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv"]].values
labels = [f"{r.chrom}:{r.pos:,} {r.svtype}{r.svlen}  F_ST={r.max_FST_pair_sv:.2f}" for r in top20_cat.itertuples()]
im = ax.imshow(data, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=1)
ax.set_xticks(range(5)); ax.set_xticklabels(POPS)
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
ax.set_title("Top-20 high-F_ST SVs in v3 catalog\nPer-super-pop AF (red = high, blue = low)")
for i in range(len(labels)):
    for j in range(5):
        ax.text(j, i, f"{data[i,j]:.2f}", ha="center", va="center", fontsize=6,
                 color="white" if data[i,j] > 0.5 else "black")
fig.colorbar(im, label="AF in super-pop")
plt.tight_layout()
fig.savefig(AUDIT_DIR / "top_fst_heatmap.jpg", format="jpg", dpi=200, bbox_inches="tight", pil_kwargs={"quality": 92})
plt.close(fig)

print(f"\n[OK] all artifacts written to {AUDIT_DIR}/")
print("Files:")
for f in sorted(AUDIT_DIR.glob("*")):
    print(f"  {f.name}  ({f.stat().st_size} bytes)")
