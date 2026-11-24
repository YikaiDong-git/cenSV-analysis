#!/usr/bin/env python3
import os, json, struct
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
DISC = f"{ROOT}/analyses/cdr_qtl"
SUB = f"{ROOT}/subproject_methylation"
PHENO = f"{DISC}/output/per_hap_chrom_phenotypes.tsv"
GRM = f"{SUB}/results/greml/pericen_1m_all"
PED = f"{ROOT}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
SVDROOT = (f"{ROOT}/sv_extraction/decompose/"
           "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200")
OUTD = f"{DISC}/results/figures/05_confound"; os.makedirs(OUTD, exist_ok=True)
ANCHORS = [
    ("chr10", 41048073, "DEL", 340, "cdr_pos_offset", "chr10 DEL-340 cdr_pos_offset (p=2e-6)"),
    ("chr10", 41640398, "INS", 339, "cdr_pos_iqr",    "chr10 INS-339 cdr_pos_iqr (p=5e-6)"),
]
SPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
SP_COL = {"AFR": "#E24A33", "AMR": "#348ABD", "EAS": "#988ED5", "EUR": "#777777", "SAS": "#FBC15E"}
DOS_COL = {0: "#2c7fb8", 1: "#7fcdbb", 2: "#edf8b1"}
MIN_MAF, MIN_N, MIN_POPS = 0.05, 12, 3


def read_grm_pcs(prefix, npc=10):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]; K[i, j] = v; K[j, i] = v
    ids = [l.split()[1] for l in open(f"{prefix}.grm.id")]
    ev, evec = np.linalg.eigh(K); o = np.argsort(ev)[::-1]
    P = pd.DataFrame(evec[:, o[:npc]], columns=[f"PC{i+1}" for i in range(npc)]); P["sample"] = ids
    return P.set_index("sample")


def rint(a):
    r = stats.rankdata(a); return stats.norm.ppf((r - 0.5) / len(r))


def ols(y, Xcols):
    X = np.column_stack([np.ones(len(y))] + Xcols)
    m = ~np.isnan(y) & ~np.isnan(X).any(axis=1)
    Xm, ym, n, k = X[m], y[m], int(m.sum()), X.shape[1]
    if n <= k + 1 or Xm[:, 1].std() == 0:
        return None
    b = np.linalg.lstsq(Xm, ym, rcond=None)[0]; r = ym - Xm @ b
    s2 = (r ** 2).sum() / (n - k); se = np.sqrt(np.diag(s2 * np.linalg.inv(Xm.T @ Xm)))
    t = b / se; p = 2 * (1 - stats.t.cdf(np.abs(t), n - k))
    return dict(beta=b, se=se, t=t, p=p, n=n)


def re_meta(betas, ses):
    betas, ses = np.array(betas), np.array(ses)
    w = 1 / ses ** 2; bfe = (w * betas).sum() / w.sum(); Q = (w * (betas - bfe) ** 2).sum()
    df = len(betas) - 1; I2 = max(0, (Q - df) / Q) if Q > 0 else 0
    tau2 = max(0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum())) if w.sum() else 0
    wre = 1 / (ses ** 2 + tau2); bre = (wre * betas).sum() / wre.sum(); sere = np.sqrt(1 / wre.sum())
    return dict(beta=float(bre), z=float(bre / sere), I2=float(I2), pQ=float(1 - stats.chi2.cdf(Q, df)))


ph = pd.read_csv(PHENO, sep="\t")
PCs = read_grm_pcs(GRM)
ped = pd.read_csv(PED, sep=r"\s+")
sc = [c for c in ped.columns if c.lower() in ("sampleid", "sample", "individual", "iid")][0]
pcc = [c for c in ped.columns if "super" in c.lower()][0]
spmap = ped.set_index(sc)[pcc].to_dict()
PHCOL = {"cdr_pos_offset": "off", "cdr_pos_iqr": "iqr"}

results = {}
fig, axes = plt.subplots(len(ANCHORS), 4, figsize=(17, 4.0 * len(ANCHORS)), squeeze=False)
for ai, (chrom, pos, svt, svl, pheno, lab) in enumerate(ANCHORS):
    phc = ph[ph.chm13_chrom == chrom]
    ind = phc.groupby("sample").agg(off=("cdr_pos_offset", "mean"), iqr=("cdr_pos_iqr", "mean"),
                                    alen=("active_len", "mean"), cwidth=("median_cdr_width", "mean")).reset_index()
    npz = np.load(f"{SVDROOT}/{chrom}/{chrom}_sv_genotypes.npz")
    key = "arr_0" if "arr_0" in npz else list(npz.keys())[0]
    G = npz[key]; sv_samp = open(f"{SVDROOT}/{chrom}/{chrom}_sv_samples.txt").read().split()
    svv = pd.read_csv(f"{SVDROOT}/{chrom}/{chrom}_sv_variants.tsv", sep="\t")
    if G.shape[0] != len(sv_samp) and G.shape[1] == len(sv_samp):
        G = G.T
    cand = svv[(svv.pos == pos) & (svv.svtype.astype(str).str.startswith(svt)) & (svv.svlen.between(svl - 3, svl + 3))]
    dos = pd.Series(np.rint(G[:, cand.index[0]].astype(float)), index=sv_samp).clip(0, 2)
    d = ind.merge(dos.rename("dos"), left_on="sample", right_index=True).merge(PCs, left_on="sample", right_index=True)
    d["sp"] = d["sample"].map(spmap)
    pcm = [d[f"PC{i+1}"].values for i in range(10)]
    dv = d["dos"].values
    ycol = PHCOL[pheno]
    is_off = (pheno == "cdr_pos_offset")
    cdr_from_start_kb = (d["alen"].values / 2.0) * (1 + d["off"].values) / 1000.0
    T1 = ols(rint(d["alen"].values), [dv] + pcm)
    T2 = ols(rint(d["cwidth"].values), [dv] + pcm)
    T3 = ols(rint(d[ycol].values), [dv] + pcm)
    T4 = ols(rint(d[ycol].values), [dv, rint(d["alen"].values)] + pcm)
    TA = ols(rint(cdr_from_start_kb), [dv] + pcm) if is_off else None
    atten = 100 * (1 - T4["beta"][1] / T3["beta"][1])
    rows, dropped = [], []
    for s in SPS:
        idx = d.index[d.sp == s]
        if len(idx) < MIN_N:
            dropped.append(f"{s}(n={len(idx)})"); continue
        maf = min(d.loc[idx, "dos"].mean() / 2, 1 - d.loc[idx, "dos"].mean() / 2)
        if maf < MIN_MAF:
            dropped.append(f"{s}(MAF={maf:.3f})"); continue
        rr = ols(rint(d.loc[idx, ycol].values), [d.loc[idx, "dos"].values])
        if rr:
            rows.append(dict(sp=s, beta=float(rr["beta"][1]), se=float(rr["se"][1]), maf=float(maf), n=rr["n"]))
    gate, meta, loso = "REJECT(<3 pops)", None, {}
    if len(rows) >= MIN_POPS:
        betas = [r["beta"] for r in rows]; ses = [r["se"] for r in rows]
        meta = re_meta(betas, ses)
        flip = bool((np.sign(betas) != np.sign(betas[int(np.argmax([1 / s**2 for s in ses]))])).any())
        for k in range(len(rows)):
            sub = [j for j in range(len(rows)) if j != k]
            loso[rows[k]["sp"]] = re_meta([betas[j] for j in sub], [ses[j] for j in sub])["z"]
        ok = (not flip) and meta["I2"] <= 0.5 and abs(meta["z"]) > 2 and all(np.sign(v) == np.sign(meta["z"]) for v in loso.values())
        gate = "PASS" if ok else "CONFOUND-SUSPECT"
    geom_ok = (T1["p"][1] > 0.05) and (T4["p"][1] < 0.05) and abs(atten) < 50
    abs_ok = bool(is_off and TA and (TA["p"][1] < 0.05) and (np.sign(TA["beta"][1]) == np.sign(T3["beta"][1])))
    overall = "ROBUST" if (gate == "PASS" and geom_ok and (abs_ok or not is_off)) else "PARTIAL/REVIEW"
    results[lab] = dict(chrom=chrom, pos=pos, phenotype=pheno, af=float(dv.mean()/2), n=int(T3["n"]),
                        T1_len_p=float(T1["p"][1]), T2_width_p=float(T2["p"][1]),
                        T3_t=float(T3["t"][1]), T3_p=float(T3["p"][1]),
                        T4_t=float(T4["t"][1]), atten_pct=float(atten),
                        ABS_t=(float(TA["t"][1]) if TA else None), ABS_p=(float(TA["p"][1]) if TA else None),
                        gate=gate, meta=meta, loso={k: round(v, 2) for k, v in loso.items()},
                        dropped=dropped, geom_ok=bool(geom_ok), abs_ok=abs_ok, overall=overall)

    def box(ax, vals, ttl, ylab):
        bp = ax.boxplot([vals[dv == g] for g in [0, 1, 2]], positions=[0, 1, 2], widths=0.6,
                        patch_artist=True, showfliers=False)
        for patch, g in zip(bp["boxes"], [0, 1, 2]):
            patch.set_facecolor(DOS_COL[g]); patch.set_alpha(0.85)
        ax.set_xticks([0, 1, 2]); ax.set_xlabel("dosage"); ax.set_ylabel(ylab); ax.set_title(ttl, fontsize=8)
    box(axes[ai, 0], d[ycol].values, f"{lab}\n{pheno}~dos t={T3['t'][1]:+.1f} (+len t={T4['t'][1]:+.1f})", pheno)
    if is_off:
        box(axes[ai, 1], cdr_from_start_kb, f"ABS: CDR-from-start (kb)\n~dos t={TA['t'][1]:+.1f} p={TA['p'][1]:.1e}", "kb from start")
    else:
        box(axes[ai, 1], d["cwidth"].values, f"T2: median CDR width~dos\np={T2['p'][1]:.2f}", "median width")
    box(axes[ai, 2], d["alen"].values / 1e6, f"T1: array length~dos\np={T1['p'][1]:.2f}", "array len (Mb)")
    axf = axes[ai, 3]
    if meta:
        yy = np.arange(len(rows) + 1)[::-1]
        bb = [r["beta"] for r in rows] + [meta["beta"]]; ee = [r["se"] for r in rows] + [abs(meta["beta"]/meta["z"]) if meta["z"] else 0]
        cols = [SP_COL[r["sp"]] for r in rows] + ["#d95f0e"]
        for i in range(len(bb)):
            axf.errorbar(bb[i], yy[i], xerr=1.96 * ee[i], fmt="o", color=cols[i], capsize=3, ms=6)
        axf.axvline(0, color="black", lw=0.8)
        axf.set_yticks(yy); axf.set_yticklabels([f"{r['sp']}(MAF{r['maf']:.2f})" for r in rows] + [f"RE z={meta['z']:+.1f}"], fontsize=7)
        axf.set_title(f"GATE={gate} I2={meta['I2']:.2f}\noverall={overall}", fontsize=8)
        axf.set_xlabel("beta")
    else:
        axf.text(0.5, 0.5, f"{gate}\ndropped {dropped}", ha="center", fontsize=7); axf.axis("off")

fig.suptitle("chr10 CDR-QTL leads (SUGGESTIVE): boundary + ABS-coord + strict ancestry gate", fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(f"{OUTD}/chr10_gate.jpg", format="jpg", dpi=140, pil_kwargs={"quality": 90})
json.dump(results, open(f"{OUTD}/chr10_gate.json", "w"), indent=2, default=str)
print(f"[fig] {OUTD}/chr10_gate.jpg")
for lab, r in results.items():
    print(f"\n=== {lab} === OVERALL={r['overall']}  (n={r['n']}, af={r['af']:.2f}, pheno={r['phenotype']})")
    print(f"  geom: len~dos p={r['T1_len_p']:.2f} width~dos p={r['T2_width_p']:.2f} pheno~dos t={r['T3_t']:+.1f}->+len {r['T4_t']:+.1f} (atten {r['atten_pct']:.0f}%) geom_ok={r['geom_ok']}")
    print(f"  ABS-coord: t={r['ABS_t']} p={r['ABS_p']} abs_ok={r['abs_ok']}")
    print(f"  gate={r['gate']} meta={r['meta']} LOSO={r['loso']} dropped={r['dropped']}")
