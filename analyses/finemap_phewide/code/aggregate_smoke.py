#!/usr/bin/env python3
"""Aggregate the single-SV sweep smoke diag TSVs into a Table-4-style side-by-side + a digest."""
import sys, glob
import numpy as np, pandas as pd

diag_glob, out_wide, out_digest = sys.argv[1], sys.argv[2], sys.argv[3]
df = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(diag_glob)], ignore_index=True)
df["am"] = df["arm"] + "_" + df["mode"]
key = ["sv_locus", "trait", "region", "cv_r2"]

pip = df.pivot_table(index=key, columns="am", values="sv_pip", aggfunc="first")
zsv = df.pivot_table(index=key, columns="arm", values="z_SV", aggfunc="first").add_prefix("zSV_")
snpr = df.pivot_table(index=key, columns="mode", values="svSNP_r", aggfunc="first").add_prefix("svSNPr_")
W = pip.join(zsv).join(snpr).reset_index()
for c in ["matched_allwindow", "matched_hitbundle", "impute_allwindow", "impute_hitbundle"]:
    if c not in W: W[c] = np.nan
W = W.rename(columns={"matched_allwindow": "pip_m_aw", "matched_hitbundle": "pip_m_hb",
                      "impute_allwindow": "pip_i_aw", "impute_hitbundle": "pip_i_hb"})
W["zSV_max"] = W[["zSV_matched", "zSV_impute"]].abs().max(axis=1)
W.to_csv(out_wide, sep="\t", index=False, float_format="%.4f")

d = []
d.append(f"# single-SV sweep SMOKE aggregate  (cells={len(W)} SVxtrait, converged={int(df['converged'].astype(str).eq('True').mean()*100)}%)")
d.append(f"per-mode wall_s (mean): " + " ".join(f"{k}={v:.1f}s" for k, v in df.groupby('am')['wall_s'].mean().items()))
r = W["svSNPr_hitbundle"].dropna()
d.append(f"\n## Identifiability (hit_bundle svSNP_r = max|corr(SV,SNP)|):  median={r.median():.3f}  "
         f"distinct(<0.85)={100*(r<0.85).mean():.1f}%  tagged(>=0.95)={100*(r>=0.95).mean():.1f}%")
sub = W.dropna(subset=["pip_m_hb", "pip_i_hb"])
flip = ((sub["pip_m_hb"] > 0.5) & (sub["pip_i_hb"] < 0.1))
d.append(f"## matched<->imputeG hit_bundle PIP:  corr={sub['pip_m_hb'].corr(sub['pip_i_hb']):.3f}  "
         f"matched-drop OVER-CALL cells (m_hb>0.5 & i_hb<0.1)={int(flip.sum())}  "
         f"reverse(i>0.5&m<0.1)={int(((sub['pip_i_hb']>0.5)&(sub['pip_m_hb']<0.1)).sum())}")
d.append(f"## all_window vs hit_bundle (matched):  mean pip_aw={W['pip_m_aw'].mean():.3f}  "
         f"mean pip_hb={W['pip_m_hb'].mean():.3f}  (aw inflated by SV-vs-SV degeneracy if >> hb)")

d.append("\n## Anchor SVs — strongest-|z_SV| trait, 4-way PIP [m_aw / m_hb / i_aw / i_hb], svSNP_r_hb")
anchors = ["chr4_sv349", "chr11_sv263", "chr3_sv20", "chr5_sv324", "chr8_sv5", "chr16_sv1029",
           "chr8_sv575", "chr4_sv572", "chr16_sv71", "chr20_sv87"]
W["phys"] = W["sv_locus"].str.split("_").str[:2].str.join("_")
for a in anchors:
    s = W[W["phys"] == a]
    if s.empty: continue
    top = s.loc[s["zSV_max"].idxmax()]
    d.append(f"  {a:13s} cv={top['cv_r2']:.2f}  {top['trait']:14s} z_SV={top['zSV_max']:+.2f}  "
             f"PIP[{top['pip_m_aw']:.2f}/{top['pip_m_hb']:.2f}/{top['pip_i_aw']:.2f}/{top['pip_i_hb']:.2f}]  "
             f"svSNPr_hb={top['svSNPr_hitbundle']:.2f}")

cand = W[(W["zSV_max"] > 5) & (W["svSNPr_hitbundle"] < 0.85) &
         ((W["pip_m_hb"] > 0.5) | (W["pip_i_hb"] > 0.5))].sort_values("zSV_max", ascending=False)
d.append(f"\n## CANDIDATE cells (|z_SV|>5, svSNP_r<0.85 distinct, hit_bundle PIP>0.5 in >=1 arm): n={len(cand)}")
for _, r in cand.head(12).iterrows():
    d.append(f"  {r['phys']:13s} {r['trait']:14s} z={r['zSV_max']:+.2f} svSNPr={r['svSNPr_hitbundle']:.2f} "
             f"m_hb={r['pip_m_hb']:.2f} i_hb={r['pip_i_hb']:.2f} cv={r['cv_r2']:.2f}")

open(out_digest, "w").write("\n".join(d) + "\n")
print("\n".join(d))
print(f"\n[wrote] {out_wide} ({len(W)} cells) + {out_digest}")
