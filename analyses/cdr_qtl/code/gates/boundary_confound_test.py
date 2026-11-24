#!/usr/bin/env python3
import os, json
import numpy as np
import pandas as pd
import matplotlib
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import struct

ROOT = CENSV_ROOT + ""
DISC = f"{ROOT}/analyses/cdr_qtl"
SUB = f"{ROOT}/subproject_methylation"
PHENO = f"{DISC}/output/per_hap_chrom_phenotypes.tsv"
GRM = f"{SUB}/results/greml/pericen_1m_all"
PED = f"{ROOT}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
SVDROOT = (f"{ROOT}/sv_extraction/decompose/"
           "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200")
OUTD = f"{DISC}/results/figures/05_confound"; os.makedirs(OUTD, exist_ok=True)
LEADS = [(44300221, "DEL", 113, "DEL-113"), (44302147, "INS", 59, "INS-59")]
DOS_COL = {0: "#2c7fb8", 1: "#7fcdbb", 2: "#edf8b1"}


def read_grm_pcs(prefix, npc=10):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]; K[i, j] = v; K[j, i] = v
    ids = [l.split()[1] for l in open(f"{prefix}.grm.id")]
    ev, evec = np.linalg.eigh(K); order = np.argsort(ev)[::-1]
    P = pd.DataFrame(evec[:, order[:npc]], columns=[f"PC{i+1}" for i in range(npc)])
    P["sample"] = ids
    return P.set_index("sample")


def ols(y, Xcols):
    X = np.column_stack([np.ones(len(y))] + Xcols)
    m = ~np.isnan(y) & ~np.isnan(X).any(axis=1)
    Xm, ym, n, k = X[m], y[m], int(m.sum()), X.shape[1]
    beta = np.linalg.lstsq(Xm, ym, rcond=None)[0]
    resid = ym - Xm @ beta; s2 = (resid ** 2).sum() / (n - k)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(Xm.T @ Xm)))
    t = beta / se; p = 2 * (1 - stats.t.cdf(np.abs(t), n - k))
    return dict(beta=beta, se=se, t=t, p=p, n=n)


def rint(a):
    r = stats.rankdata(a); return stats.norm.ppf((r - 0.5) / len(r))


ph = pd.read_csv(PHENO, sep="\t")
ph8 = ph[ph.chm13_chrom == "chr8"]
PCs = read_grm_pcs(GRM)
ped = pd.read_csv(PED, sep=r"\s+")
sc = [c for c in ped.columns if c.lower() in ("sampleid", "sample", "individual", "iid")][0]
pcc = [c for c in ped.columns if "super" in c.lower()][0]
spmap = ped.set_index(sc)[pcc].to_dict()

ind = ph8.groupby("sample").agg(off=("cdr_pos_offset", "mean"), alen=("active_len", "mean"),
                                cwidth=("median_cdr_width", "mean")).reset_index()

npz = np.load(f"{SVDROOT}/chr8/chr8_sv_genotypes.npz")
key = "arr_0" if "arr_0" in npz else list(npz.keys())[0]
G = npz[key]; sv_samp = open(f"{SVDROOT}/chr8/chr8_sv_samples.txt").read().split()
svv = pd.read_csv(f"{SVDROOT}/chr8/chr8_sv_variants.tsv", sep="\t")
if G.shape[0] != len(sv_samp) and G.shape[1] == len(sv_samp):
    G = G.T

results = {}
fig, axes = plt.subplots(2, len(LEADS) * 2, figsize=(7 * len(LEADS), 8))
for li, (pos, svt, svl, lab) in enumerate(LEADS):
    cand = svv[(svv.pos == pos) & (svv.svtype.astype(str).str.startswith(svt)) &
               (svv.svlen.between(svl - 3, svl + 3))]
    dos = pd.Series(np.rint(G[:, cand.index[0]].astype(float)), index=sv_samp).clip(0, 2)
    d = ind.merge(dos.rename("dos"), left_on="sample", right_index=True)
    d = d.merge(PCs, left_on="sample", right_index=True)
    pcm = [d[f"PC{i+1}"].values for i in range(10)]
    af = d["dos"].mean() / 2
    yoff = rint(d["off"].values); alen = d["alen"].values; dv = d["dos"].values
    T1 = ols(rint(alen), [dv] + pcm)
    T2 = ols(rint(d["cwidth"].values), [dv] + pcm)
    T3 = ols(yoff, [dv] + pcm)
    T4 = ols(yoff, [dv, rint(alen)] + pcm)
    atten = 100 * (1 - T4["beta"][1] / T3["beta"][1])
    results[lab] = dict(af=float(af), n=int(T3["n"]),
                        T1_len_beta=float(T1["beta"][1]), T1_len_p=float(T1["p"][1]),
                        T2_width_beta=float(T2["beta"][1]), T2_width_p=float(T2["p"][1]),
                        T3_off_beta=float(T3["beta"][1]), T3_off_t=float(T3["t"][1]), T3_off_p=float(T3["p"][1]),
                        T4_off_beta_adj=float(T4["beta"][1]), T4_off_t_adj=float(T4["t"][1]), T4_off_p_adj=float(T4["p"][1]),
                        attenuation_pct=float(atten),
                        r_dos_len=float(np.corrcoef(dv, alen)[0, 1]))
    axA = axes[0, li * 2]
    for g in [0, 1, 2]:
        m = dv == g
        axA.scatter(alen[m] / 1e6, d["off"].values[m], s=18, color=DOS_COL[g], label=f"{g} (n={m.sum()})", alpha=0.8)
    axA.set_xlabel("active-array length (Mb)"); axA.set_ylabel("cdr_pos_offset")
    axA.set_title(f"{lab}: offset vs array length\nby dosage", fontsize=8); axA.legend(fontsize=6)
    axB = axes[0, li * 2 + 1]
    bp = axB.boxplot([alen[dv == g] / 1e6 for g in [0, 1, 2]], positions=[0, 1, 2], widths=0.6,
                     patch_artist=True, showfliers=False)
    for patch, g in zip(bp["boxes"], [0, 1, 2]):
        patch.set_facecolor(DOS_COL[g])
    axB.set_xticks([0, 1, 2]); axB.set_xlabel(f"{lab} dosage"); axB.set_ylabel("array length (Mb)")
    axB.set_title(f"T1: length~dosage\nbeta={T1['beta'][1]:+.2f} p={T1['p'][1]:.1e}", fontsize=8)
    axC = axes[1, li * 2]
    res_off = yoff - (ols(yoff, [rint(alen)])["beta"][0] + ols(yoff, [rint(alen)])["beta"][1] * rint(alen))
    bp1 = axC.boxplot([d["off"].values[dv == g] for g in [0, 1, 2]], positions=[0, 1, 2], widths=0.35,
                      patch_artist=True, showfliers=False)
    for patch, g in zip(bp1["boxes"], [0, 1, 2]):
        patch.set_facecolor(DOS_COL[g]); patch.set_alpha(0.8)
    axC.set_xticks([0, 1, 2]); axC.set_xlabel(f"{lab} dosage"); axC.set_ylabel("cdr_pos_offset")
    axC.set_title(f"T3 raw beta={T3['beta'][1]:+.2f} t={T3['t'][1]:+.1f}\n"
                  f"T4 +length beta={T4['beta'][1]:+.2f} t={T4['t'][1]:+.1f}  (atten {atten:.0f}%)", fontsize=8)
    axD = axes[1, li * 2 + 1]
    bp2 = axD.boxplot([d["cwidth"].values[dv == g] for g in [0, 1, 2]], positions=[0, 1, 2], widths=0.6,
                      patch_artist=True, showfliers=False)
    for patch, g in zip(bp2["boxes"], [0, 1, 2]):
        patch.set_facecolor(DOS_COL[g])
    axD.set_xticks([0, 1, 2]); axD.set_xlabel(f"{lab} dosage"); axD.set_ylabel("median CDR width (bp)")
    axD.set_title(f"T2: width~dosage\nbeta={T2['beta'][1]:+.2f} p={T2['p'][1]:.1e}", fontsize=8)

fig.suptitle("chr8 boundary/geometry confound test — is the CDR-position QTL real, or array-architecture?",
             fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{OUTD}/chr8_boundary_confound.jpg", format="jpg", dpi=145, pil_kwargs={"quality": 90})
print(f"[fig] {OUTD}/chr8_boundary_confound.jpg")
json.dump(results, open(f"{OUTD}/boundary_confound_results.json", "w"), indent=2)
print(json.dumps(results, indent=2))
for lab, r in results.items():
    length_route = r["T1_len_p"] < 0.05
    survives = r["T4_off_p_adj"] < 0.05 and abs(r["attenuation_pct"]) < 50
    verdict = ("ARCHITECTURE/MECHANICAL (length~dosage sig AND offset effect attenuates)"
               if length_route and not survives else
               "REAL POSITION (offset survives length adjustment)" if survives and not length_route else
               "MIXED/INDETERMINATE")
    print(f"\nVERDICT {lab}: {verdict}")
    print(f"  length~dosage p={r['T1_len_p']:.2e}  offset raw t={r['T3_off_t']:+.2f} -> +length t={r['T4_off_t_adj']:+.2f} (atten {r['attenuation_pct']:.0f}%)")
