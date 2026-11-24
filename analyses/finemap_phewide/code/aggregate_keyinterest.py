#!/usr/bin/env python3
import sys, glob, re
import numpy as np, pandas as pd

diag_glob, svs_tsv, traits_tsv, out_wide, out_digest, out_cand = sys.argv[1:7]

AOU = {
    "chr4_sv398":  "platelet MPV/PDW: nominal, SURVIVES conditioning",
    "chr5_sv213":  "finemap underpowered survivor (chr5 vertigo twin)",
    "chr8_sv2":    "collinear twin of chr8_sv5 r=0.77 (hearing region)",
    "chr16_sv465": "psychoactive substance: robust p=2.1e-15",
    "chr16_sv37":  "psychoactive/nicotine: robust cv0.46",
    "chr16_sv73":  "nicotine dep: robust",
    "chr16_sv1017":"psychoactive: robust_low_impute",
    "chr16_sv559": "psychoactive: robust",
    "chr16_sv396": "alcohol: PIP=1.0 TAGGING ARTIFACT (null marginal p=1.0)",
    "chr10_sv697": "psoriatic arthropathy p=1.4e-6, cv0.765 WELL-IMPUTED",
    "chr10_sv698": "psoriatic arthropathy (twin of sv697)",
    "chr11_sv200": "prothrombin/F2 tag p=7e-9 (F2 ~9Mb away)",
    "chr11_sv400": "prothrombin/F2 tag",
    "chr11_sv148": "prothrombin/F2 tag",
    "chr11_sv389": "congenital hip p=1.2e-6",
    "chr17_sv317": "vesicoureteral reflux (score_inflated n=72)",
    "chr17_sv331": "vesicoureteral reflux (score_inflated)",
    "chr1_sv1325": "congenital genital/spine (score_inflated)",
    "chr11_sv271": "albumin: UNRESOLVED (EUR p=0.14 sign-disc; external GWAS pending)",
    "chr12_sv721": "height: NEGATIVE (discordant, cv0.15)",
    "chr1_sv650":  "breast: NEGATIVE (LD-shadow, 61.5% support dropped)",
    "chr19_sv1949":"SNP-equiv NEGATIVE, cv~1.0 (near-perfect SNP combo = perfect-tag ANCHOR)",
    "chr11_sv121": "hypertension: x-ancestry q=0.009 but fine-map NULL",
    "chr3_sv40":   "PROS1->VTE: positive control",
}
ORDER = ["chr4_sv398","chr5_sv213","chr8_sv2","chr16_sv465","chr16_sv37","chr16_sv73","chr16_sv1017",
         "chr16_sv559","chr16_sv396","chr10_sv697","chr10_sv698","chr11_sv200","chr11_sv400",
         "chr11_sv148","chr11_sv389","chr17_sv317","chr17_sv331","chr1_sv1325","chr11_sv271",
         "chr12_sv721","chr1_sv650","chr19_sv1949","chr11_sv121","chr3_sv40"]
EXPECT_TAG = {"chr16_sv396","chr1_sv650","chr19_sv1949","chr12_sv721"}

files = sorted(glob.glob(diag_glob))
df = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
df["am"] = df["arm"] + "_" + df["mode"]
df["phys"] = df["sv_locus"].str.split("_").str[:2].str.join("_")
key = ["sv_locus", "phys", "trait", "region", "cv_r2"]

pip = df.pivot_table(index=key, columns="am", values="sv_pip", aggfunc="first")
zsv = df.pivot_table(index=key, columns="arm", values="z_SV", aggfunc="first").add_prefix("zSV_")
snpr = df.pivot_table(index=key, columns="mode", values="svSNP_r", aggfunc="first").add_prefix("svSNPr_")
W = pip.join(zsv).join(snpr).reset_index()
for c in ["matched_allwindow", "matched_hitbundle", "impute_allwindow", "impute_hitbundle"]:
    if c not in W: W[c] = np.nan
W = W.rename(columns={"matched_allwindow": "pip_m_aw", "matched_hitbundle": "pip_m_hb",
                      "impute_allwindow": "pip_i_aw", "impute_hitbundle": "pip_i_hb"})
W["zSV_max"] = W[["zSV_matched", "zSV_impute"]].abs().max(axis=1)
W["div_hb"] = W["pip_m_hb"] - W["pip_i_hb"]
W["diverge"] = W["div_hb"].abs() > 0.4
W.to_csv(out_wide, sep="\t", index=False, float_format="%.4f")

svmap = {}
for ln in open(svs_tsv).read().splitlines()[1:]:
    f = ln.split("\t"); svmap[f[1]] = (f[3], f[4])
traits_bgz = dict(l.split("\t")[:2] for l in open(traits_tsv).read().splitlines() if l.strip())


def parse_region(r):
    m = re.match(r"chr(\d+)_super\d+_(\d+)_(\d+)", r); return m.group(1), m.group(2), m.group(3)


is_cand = (W["zSV_max"] > 4) & (W["svSNPr_hitbundle"] < 0.90) & ((W["pip_m_hb"] > 0.3) | (W["pip_i_hb"] > 0.3))
sel = W[is_cand | W["diverge"]].copy()
rows = []
for _, r in sel.iterrows():
    if r["phys"] not in svmap: continue
    rb, col = svmap[r["phys"]]; chrom, start, end = parse_region(r["region"])
    reason = "candidate" if (is_cand.loc[r.name] if r.name in is_cand.index else False) else "divergence"
    if r["diverge"] and not (r["zSV_max"] > 4): reason = "divergence"
    for arm in ("matched", "impute"):
        rows.append([r["phys"], r["sv_locus"], r["region"], rb, chrom, start, end, col, arm, r["trait"],
                     traits_bgz.get(r["trait"], ""), "hitbundle",
                     f"{(r['pip_m_hb'] if arm=='matched' else r['pip_i_hb']):.4f}",
                     f"{r['svSNPr_hitbundle']:.4f}", f"{r['zSV_max']:.3f}", reason])
cc = pd.DataFrame(rows, columns=["phys", "sv_locus", "region", "region_bundle", "chrom", "start", "end",
                                 "col", "arm", "trait", "bgz", "mode", "pip", "svSNP_r", "z_SV", "reason"])
cc.to_csv(out_cand, sep="\t", index=False)

d = []
converged = int(df["converged"].astype(str).eq("True").mean() * 100)
d.append(f"# KEY-INTEREST single-SV sweep  (cells={len(W)} SVxtrait, diag_files={len(files)}, converged={converged}%)")
d.append("# PIP is NON-IDENTIFIABLE at svSNP_r->1 (project theorem); read = svSNP_r (distinct<0.85) + z_SV (burden) + arm-concordance.")
r = W["svSNPr_hitbundle"].dropna()
d.append(f"## Identifiability (hitbundle svSNP_r): median={r.median():.3f} distinct(<0.85)={100*(r<0.85).mean():.0f}% tagged(>=0.95)={100*(r>=0.95).mean():.0f}%")
sub = W.dropna(subset=["pip_m_hb", "pip_i_hb"])
d.append(f"## matched<->imputeG hitbundle PIP: corr={sub['pip_m_hb'].corr(sub['pip_i_hb']):.3f} "
         f"divergence_cells(|d|>0.4)={int(W['diverge'].sum())} "
         f"over-call(m>0.5&i<0.1)={int(((sub['pip_m_hb']>0.5)&(sub['pip_i_hb']<0.1)).sum())} "
         f"reverse(i>0.5&m<0.1)={int(((sub['pip_i_hb']>0.5)&(sub['pip_m_hb']<0.1)).sum())}")

d.append("\n## Per key-interest SV -- strongest-|z| trait | PIP[m_aw/m_hb/i_aw/i_hb] | svSNP_r_hb | read vs AoU")
for phys in ORDER:
    s = W[W["phys"] == phys]
    if s.empty:
        d.append(f"  {phys:13s}  (no cells -- trait absent / filtered)  | {AOU.get(phys,'')}"); continue
    top = s.loc[s["zSV_max"].idxmax()]
    tagged = top["svSNPr_hitbundle"] >= 0.85
    read = "TAGGED " if tagged else "DISTINCT"
    flag = ""
    if phys in EXPECT_TAG:
        flag = "  [OK: tag as expected]" if tagged else "  [!! distinct, AoU expected tag]"
    d.append(f"  {phys:13s} cv={top['cv_r2']:.2f} {str(top['trait'])[:16]:16s} z={top['zSV_max']:+.2f} "
             f"PIP[{top['pip_m_aw']:.2f}/{top['pip_m_hb']:.2f}/{top['pip_i_aw']:.2f}/{top['pip_i_hb']:.2f}] "
             f"snpR={top['svSNPr_hitbundle']:.2f} {read}{flag} | {AOU.get(phys,'')}")

dv = W[W["diverge"]].copy()
dv = dv.sort_values("div_hb", key=lambda s: s.abs(), ascending=False) if len(dv) else dv
d.append(f"\n## matched<->imputeG DIVERGENCE cells (|pip_m_hb - pip_i_hb|>0.4): n={len(dv)}")
for _, x in dv.head(15).iterrows():
    d.append(f"  {str(x['phys']):13s} {str(x['trait'])[:16]:16s} m_hb={x['pip_m_hb']:.2f} i_hb={x['pip_i_hb']:.2f} "
             f"div={x['div_hb']:+.2f} z={x['zSV_max']:+.2f} snpR={x['svSNPr_hitbundle']:.2f}")

cd = W[is_cand].sort_values("zSV_max", ascending=False)
d.append(f"\n## CANDIDATE cells (|z_SV|>4, svSNP_r<0.90 distinct, hitbundle PIP>0.3 in >=1 arm): n={len(cd)}")
for _, x in cd.head(20).iterrows():
    d.append(f"  {str(x['phys']):13s} {str(x['trait'])[:16]:16s} z={x['zSV_max']:+.2f} snpR={x['svSNPr_hitbundle']:.2f} "
             f"m_hb={x['pip_m_hb']:.2f} i_hb={x['pip_i_hb']:.2f} cv={x['cv_r2']:.2f}")

d.append(f"\n[s_rss to compute on {len(cc)//2} candidate/divergence cells] -> {out_cand}")
d.append("per-mode wall_s(mean): " + " ".join(f"{k}={v:.1f}" for k, v in df.groupby('am')['wall_s'].mean().items()))
open(out_digest, "w").write("\n".join(d) + "\n")
print("\n".join(d))
print(f"\n[wrote] {out_wide} ({len(W)} cells) + {out_digest} + {out_cand} ({len(cc)} rows)")
