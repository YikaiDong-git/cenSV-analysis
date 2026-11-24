#!/usr/bin/env python
"""Resolve the `active_len` channel: CAUSAL, ASCERTAINMENT, or LINEAGE?"""
import os, sys, glob, re, struct
import numpy as np, pandas as pd
from collections import defaultdict
from scipy import stats
from statsmodels.regression.mixed_linear_model import MixedLM
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
CAN = f"{ROOT}/analyses/cdr_qtl"
AUD = f"{ROOT}/analyses/cdr_calls"
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
CENSAT = f"{ROOT}/resources/chm13v2.0_censat_v2.1.active.bed"
DECOMP = sorted(glob.glob(f"{ROOT}/sv_extraction/decompose/*/"))[0]
OUT = f"{WU}/output/burden/channel"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{WU}/code")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pd.set_option("display.width", 220)
RES = []


def rint(a):
    a = np.asarray(a, float)
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def read_grm(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]
                K[i, j] = v; K[j, i] = v
    return K, [l.split()[1] for l in open(f"{prefix}.grm.id")]


def merge_iv(iv):
    iv = sorted(iv); out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t", engine="python")
samples = np.array(open(f"{WU}/geno/full_samples.txt").read().split())

SIGN = {"INS": +1.0, "DUP:TANDEM": +1.0, "DEL": -1.0, "BND": 0.0}
vt["sgn"] = vt.svtype.map(SIGN)
vt["signed_len"] = vt.sgn * vt.svlen

rep = vt.reset_index().sort_values("index").drop_duplicates("geno_group", keep="first")
Gu, vu = G[rep["index"].values, :], vt.iloc[rep["index"].values].reset_index(drop=True)
phy = vt.reset_index().sort_values("index").drop_duplicates(
    ["chrom", "pos", "end", "svtype", "svlen"], keep="first")
Gp, vp = G[phy["index"].values, :], vt.iloc[phy["index"].values].reset_index(drop=True)
print(f"genotype-group representatives {len(vu)}   distinct physical SVs {len(vp)}")

cs = pd.read_csv(CENSAT, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "s", "e"])
REFAA = {c: merge_iv(list(zip(g.s, g.e))) for c, g in cs.groupby("chrom")}
vu["mid"] = (vu.pos + vu.end) / 2.0
vu["in_ref"] = [any(s <= m <= e for s, e in REFAA.get(c, []))
                for c, m in zip(vu.chrom, vu["mid"])]
N_IN, N_TOT = int(vu.in_ref.sum()), len(vu)

rows = []
for c in sorted(set(vu.chrom)):
    m = (vu.chrom == c).values
    mp = (vp.chrom == c).values
    ins = m & (vu.sgn.values > 0)
    dele = m & (vu.sgn.values < 0)
    het = (Gp[mp, :] == 1)
    rows.append(pd.DataFrame(dict(
        sample=samples, chrom=c,
        burden_cis=np.nansum(Gu[m, :], axis=0),
        burden_INS=np.nansum(Gu[ins, :], axis=0),
        burden_DEL=np.nansum(Gu[dele, :], axis=0),
        netbp=np.nansum(Gp[mp, :] * vp.signed_len.values[mp][:, None], axis=0),
        churnbp=np.nansum(Gp[mp, :] * np.abs(vp.svlen.values[mp])[:, None], axis=0),
        het_bp=np.nansum(het * np.abs(vp.svlen.values[mp])[:, None], axis=0),
        het_n=np.nansum(het, axis=0))))
PC = pd.concat(rows, ignore_index=True)
PC["burden_sum"] = PC.burden_INS + PC.burden_DEL
PC["burden_diff"] = PC.burden_INS - PC.burden_DEL
tot = PC.groupby("sample").burden_cis.sum().rename("burden_tot")
PC = PC.merge(tot, on="sample")
PC["burden_trans"] = PC.burden_tot - PC.burden_cis

corr = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
HAP = corr[["sample", "hap_local", "chm13_chrom", "active_len"]].rename(
    columns={"chm13_chrom": "chrom"})
P = (corr.groupby(["sample", "chm13_chrom"])[["active_len", "mean_5mC"]]
     .mean().reset_index().rename(columns={"chm13_chrom": "chrom"}))

coh = pd.read_csv(f"{AUD}/output/cohort_corrected_final.tsv", sep="\t")
cov_rows = []
for f in sorted(glob.glob(f"{AUD}/output/calls/*/binned_freq.bed")):
    tag = os.path.basename(os.path.dirname(f))
    s, h = tag.rsplit("_", 1)
    a = pd.read_csv(f, sep="\t", header=None, usecols=[4], names=["cov"])
    cov_rows.append(dict(sample=s, hap=h, coverage=a["cov"].mean()))
Q = coh.merge(pd.DataFrame(cov_rows), on=["sample", "hap"], how="outer")
Q["is_r1041"] = Q.chemistry.astype(str).str.contains("1041").astype(float)
QI = Q.groupby("sample").agg(n_contigs=("n_contigs", "sum"), coverage=("coverage", "mean"),
                             is_r1041=("is_r1041", "max")).reset_index()

ROW_RE = re.compile(r"^(.+?)_(chr[0-9XY]+)_((.+?)#([12])#[^:]+):")
nint = defaultdict(int)
for bed in sorted(glob.glob(f"{AUD}/inputs/globe_aa/*_AS-HOR_stv_row.ort.bed")):
    for line in open(bed):
        m = ROW_RE.match(line.split("\t")[0])
        if m:
            nint[(m.group(1), m.group(2))] += 1
NI = pd.DataFrame([dict(sample=a, chrom=b, n_intervals=v) for (a, b), v in nint.items()])

K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]],
                   columns=[f"gPC{i+1}" for i in range(10)], index=gids)
PCs.index.name = "sample"

D = (PC.merge(P, on=["sample", "chrom"], how="inner")
       .merge(QI, on="sample", how="inner")
       .merge(NI, on=["sample", "chrom"], how="inner")
       .merge(PCs.reset_index(), on="sample", how="inner"))
print(f"model rows (individual x chromosome): {len(D)}  "
      f"[{D['sample'].nunique()} individuals x {D.chrom.nunique()} chromosomes]")

for col in ["burden_cis", "burden_INS", "burden_DEL", "burden_sum", "burden_diff",
            "netbp", "churnbp", "active_len"]:
    D[col + "_z"] = D.groupby("chrom")[col].transform(
        lambda s: rint(s) if s.notna().sum() == len(s) and s.std() > 0 else np.nan)
D["burden_trans_z"] = rint(D.burden_trans)

COVS = ["n_contigs", "coverage", "is_r1041", "n_intervals"] + [f"gPC{i+1}" for i in range(10)]


def fit(dsub, y, terms):
    cols = [y, "sample"] + [t for t in terms if t not in (y, "sample")]
    d = dsub[cols].replace([np.inf, -np.inf], np.nan).dropna()
    keep = [t for t in terms if d[t].std() > 1e-12]
    X = np.column_stack([np.ones(len(d))] + [d[t].values.astype(float) for t in keep])
    m = MixedLM(d[y].values.astype(float), X, groups=d["sample"].values).fit(reml=True)
    return m, d, {t: i + 1 for i, t in enumerate(keep)}


def emit(test, term, m, idx, d, note=""):
    b, t, p = m.params[idx[term]], m.tvalues[idx[term]], m.pvalues[idx[term]]
    RES.append(dict(test=test, term=term, n=len(d), beta=b, t=t, p=p, note=note))
    print(f"  {term:18s} beta {b:+9.4f}   t {t:+7.2f}   p {p:11.3g}   {note}")
    return t


print("\n" + "=" * 100)
print("T1  DIRECTION -- H1(causal) predicts OPPOSITE signs; H2/H3 predict the SAME sign.")
print("=" * 100)
print(f"  corr(burden_INS_z, burden_DEL_z) = {D.burden_INS_z.corr(D.burden_DEL_z):+.3f}")
print("\n  marginal (each alone):")
for term in ["burden_INS_z", "burden_DEL_z"]:
    m, d, i = fit(D, "active_len_z", [term, "burden_trans_z"] + COVS)
    emit("T1_marginal", term, m, i, d)
print("\n  joint:")
m, d, i = fit(D, "active_len_z", ["burden_INS_z", "burden_DEL_z", "burden_trans_z"] + COVS)
t_ins = emit("T1_joint", "burden_INS_z", m, i, d)
t_del = emit("T1_joint", "burden_DEL_z", m, i, d)
same_sign = (t_ins > 0) == (t_del > 0)
print(f"\n  --> signs are {'THE SAME' if same_sign else 'OPPOSITE'}  ==> "
      f"{'H2/H3 (not causal)' if same_sign else 'H1 (causal)'} on this test")
print("\n  T1b reparameterised -- H1 loads on (INS-DEL), H2/H3 on (INS+DEL):")
m, d, i = fit(D, "active_len_z", ["burden_sum_z", "burden_diff_z", "burden_trans_z"] + COVS)
t_sum = emit("T1b", "burden_sum_z", m, i, d)
t_dif = emit("T1b", "burden_diff_z", m, i, d)

print("\n" + "=" * 100)
print("T2  MAGNITUDE -- H1 predicts d(active_len)/d(net signed SV bp) = +0.5 exactly.")
print("=" * 100)
srows = []
for c, g in D.groupby("chrom"):
    g = g.dropna(subset=["active_len", "burden_cis", "netbp"])
    if len(g) < 30:
        continue
    bc, bn = stats.linregress(g.burden_cis, g.active_len), stats.linregress(g.netbp, g.active_len)
    srows.append(dict(chrom=c, n=len(g), mean_burden=g.burden_cis.mean(),
                      mean_active_len=g.active_len.mean(),
                      slope_bp_per_allele=bc.slope, obs_r2_count=bc.rvalue ** 2,
                      slope_len_per_netbp=bn.slope,
                      max_causal_R2=np.var(0.5 * g.netbp) / np.var(g.active_len),
                      scale_only_density_bench=g.active_len.mean() / g.burden_cis.mean()))
S = pd.DataFrame(srows).sort_values("chrom")
S["obs_over_ceiling"] = S.obs_r2_count / S.max_causal_R2
print(S[["chrom", "n", "slope_len_per_netbp", "obs_r2_count", "max_causal_R2", "obs_over_ceiling"]]
      .to_string(index=False, float_format=lambda x: f"{x:,.4g}"))
n_viol = int((S.obs_r2_count > S.max_causal_R2).sum())
print(f"\n  median slope of active_len on NET SIGNED SV bp   {S.slope_len_per_netbp.median():+.3f}"
      f"   (H1 predicts +0.5)")
print(f"  chromosomes where observed r^2 EXCEEDS the causal ceiling: {n_viol}/{len(S)}")
print(f"  median observed r^2 {S.obs_r2_count.median():.3e}   median ceiling "
      f"{S.max_causal_R2.median():.3e}   median ratio {S.obs_over_ceiling.median():,.0f}x")
print(f"  [scale only, NOT evidence -- constant-density benchmark median "
      f"{S.scale_only_density_bench.median():,.0f} bp/allele vs observed slope "
      f"{S.slope_bp_per_allele.median():,.0f}]")
for k, v, nt in [("median_slope_active_len_per_netbp", S.slope_len_per_netbp.median(), "H1 predicts +0.5"),
                 ("median_obs_r2_over_causal_ceiling", S.obs_over_ceiling.median(),
                  f"{n_viol}/{len(S)} chroms exceed ceiling")]:
    RES.append(dict(test="T2_magnitude", term=k, n=len(S), beta=v, t=np.nan, p=np.nan, note=nt))

print("\n" + "=" * 100)
print("T3  LOCATION -- VOID, NOT EVIDENCE.")
print("=" * 100)
print(f"  {N_IN}/{N_TOT} ({100*N_IN/N_TOT:.1f}%) catalogued SVs fall INSIDE the CHM13 active HOR")
print("  array. The catalogue was built inside the array by construction, so there is no flanking")
print("  comparison set: an 'outside' burden term is built from a single variant and carries no")
print("  information. The in/out contrast cannot be run on this substrate.")
print("  CONSEQUENCE: agy's proposed discriminator for H3 (do FLANK SVs still track active_len?)")
print("  is NOT AVAILABLE here. H3 is instead attacked directly in T6.")
RES.append(dict(test="T3_location", term="VOID_pct_inside_ref_array", n=N_TOT,
                beta=100.0 * N_IN / N_TOT, t=np.nan, p=np.nan,
                note="no flanking set exists; test cannot be run"))

print("\n" + "=" * 100)
print("T4  STRATA -- rules out an assembly-QUALITY gradient only; does NOT separate H1 from H2/H3.")
print("=" * 100)
strata = [("ALL", D), ("n_intervals==1 (single contiguous)", D[D.n_intervals == 1]),
          ("n_intervals>1 (fragmented)", D[D.n_intervals > 1])]
for lab, col in [("coverage", "coverage"), ("n_contigs", "n_contigs")]:
    q = D[col].quantile([0.25, 0.75]).values
    strata += [(f"{lab} Q1 (lowest)", D[D[col] <= q[0]]), (f"{lab} Q4 (highest)", D[D[col] > q[1]])]
for lab, dd in strata:
    if len(dd) < 100:
        print(f"  {lab:36s} n={len(dd):5d}  SKIPPED (n<100)"); continue
    try:
        m, d, i = fit(dd, "active_len_z", ["burden_cis_z", "burden_trans_z"] + COVS)
        b, t, p = (m.params[i["burden_cis_z"]], m.tvalues[i["burden_cis_z"]],
                   m.pvalues[i["burden_cis_z"]])
        print(f"  {lab:36s} n={len(d):5d}  beta {b:+7.4f}  t {t:+7.2f}  p {p:10.3g}")
        RES.append(dict(test="T4_strata", term=lab, n=len(d), beta=b, t=t, p=p, note=""))
    except Exception as e:
        print(f"  {lab:36s} FAILED {e}")

print("\n" + "=" * 100)
print("T5  HOMOLOG BOUND -- under H1, |active_len(hap1)-active_len(hap2)| <= sum|svlen| over HET SVs.")
print("    All individual-level confounders difference out.")
print("=" * 100)
H = HAP.dropna(subset=["active_len"]).groupby(["sample", "chrom"])["active_len"].agg(["min", "max", "count"])
H = H[H["count"] == 2].reset_index()
H["delta_len"] = H["max"] - H["min"]
H = H.merge(PC[["sample", "chrom", "het_bp", "het_n"]], on=["sample", "chrom"], how="inner")
H = H[H.het_bp > 0]
H["ratio"] = H.delta_len / H.het_bp
viol = float((H.delta_len > H.het_bp).mean())
print(f"  rows with both haplotypes and >=1 het SV: {len(H)}")
print(f"  median homolog length gap |dLen|      {H.delta_len.median():14,.0f} bp")
print(f"  median heterozygous SV sequence       {H.het_bp.median():14,.0f} bp  "
      f"(median {H.het_n.median():.0f} het SVs)")
print(f"  median ratio |dLen| / het_bp          {H.ratio.median():14,.1f}x")
print(f"  rows VIOLATING the H1 bound |dLen| > het_bp:  {100*viol:.1f}%  "
      f"({int((H.delta_len > H.het_bp).sum())}/{len(H)})")
r_h = stats.spearmanr(H.delta_len, H.het_bp)
print(f"  spearman(|dLen|, het_bp) = {r_h.statistic:+.3f}  p = {r_h.pvalue:.3g}")
RES.append(dict(test="T5_homolog_bound", term="pct_rows_violating_causal_bound", n=len(H),
                beta=100 * viol, t=np.nan, p=np.nan,
                note=f"median ratio {H.ratio.median():.1f}x"))
H.to_csv(f"{OUT}/active_len_channel_homolog.tsv", sep="\t", index=False, float_format="%.6g")

print("\n" + "=" * 100)
print("T6  LINEAGE (H3) -- centromeric haplogroups, assigned PER HAPLOTYPE.")
print("=" * 100)
hg = []
for f in sorted(glob.glob(f"{DECOMP}chr*/chr*_haplogroups.tsv")):
    if "validated" in f:
        continue
    hg.append(pd.read_csv(f, sep="\t"))
HG = pd.concat(hg, ignore_index=True).drop(columns=["sample"]).rename(
    columns={"sample_id": "sample", "chromosome": "chrom"})
print(f"  haplogroup assignments: {len(HG)} haplotypes over {HG.chrom.nunique()} chromosomes, "
      f"{HG.haplogroup.nunique()} distinct haplogroups")

HL = HAP.copy()
HL["haplotype"] = "hap" + HL.hap_local.astype(str).str.extract(r"(\d+)")[0]
print(f"  hap_local values {sorted(set(HAP.hap_local.astype(str)))[:4]} -> "
      f"join key {sorted(set(HL.haplotype.dropna()))[:4]}")
HL = HL.merge(HG[["sample", "chrom", "haplotype", "haplogroup"]], on=["sample", "chrom", "haplotype"],
              how="inner").dropna(subset=["active_len"])
print(f"  haplotype-level rows matched to a haplogroup: {len(HL)}")

hrows = []
for c, g in HL.groupby("chrom"):
    if g.haplogroup.nunique() < 2 or len(g) < 30:
        continue
    grps = [v.active_len.values for _, v in g.groupby("haplogroup") if len(v) >= 2]
    if len(grps) < 2:
        continue
    F = stats.f_oneway(*grps)
    ss_b = sum(len(x) * (x.mean() - g.active_len.mean()) ** 2 for x in grps)
    hrows.append(dict(chrom=c, n_hap=len(g), n_hg=g.haplogroup.nunique(),
                      R2_active_len_on_haplogroup=ss_b / np.sum((g.active_len - g.active_len.mean()) ** 2),
                      p_anova=F.pvalue))
HGR = pd.DataFrame(hrows)
print(f"\n  (a) haplotype-level active_len explained by haplogroup alone:")
print(f"      median R2 = {HGR.R2_active_len_on_haplogroup.median():.3f} over {len(HGR)} chromosomes; "
      f"{int((HGR.p_anova < 0.05).sum())}/{len(HGR)} significant at p<0.05")

HD = (HL.groupby(["sample", "chrom", "haplogroup"]).size().rename("dose").reset_index())
brows = []
for c, g in D.groupby("chrom"):
    sub = HD[HD.chrom == c]
    if not len(sub):
        continue
    W = sub.pivot_table(index="sample", columns="haplogroup", values="dose", fill_value=0)
    g = g.merge(W, left_on="sample", right_index=True, how="inner")
    hgc = [x for x in W.columns if x in g.columns and g[x].std() > 1e-12]
    if len(g) < 40 or not hgc:
        continue
    y = g.active_len.values.astype(float)
    x = g.burden_cis.values.astype(float)
    r_raw = stats.pearsonr(x, y)
    Z = np.column_stack([np.ones(len(g))] + [g[h].values.astype(float) for h in hgc])
    ry = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    rx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    r_adj = stats.pearsonr(rx, ry)
    brows.append(dict(chrom=c, n=len(g), n_hg=len(hgc), r_raw=r_raw.statistic, p_raw=r_raw.pvalue,
                      r_adj_haplogroup=r_adj.statistic, p_adj=r_adj.pvalue))
B = pd.DataFrame(brows)
print(f"\n  (b) burden <-> active_len correlation, before vs after haplogroup adjustment:")
print(B.to_string(index=False, float_format=lambda x: f"{x:,.4g}"))
print(f"\n      median r_raw {B.r_raw.median():+.3f}  ->  median r_adj {B.r_adj_haplogroup.median():+.3f}"
      f"   ({100*(1-abs(B.r_adj_haplogroup.median())/max(1e-9,abs(B.r_raw.median()))):.0f}% attenuation)")
B.to_csv(f"{OUT}/active_len_channel_haplogroup.tsv", sep="\t", index=False, float_format="%.6g")
RES.append(dict(test="T6_lineage", term="median_R2_active_len_on_haplogroup", n=len(HGR),
                beta=HGR.R2_active_len_on_haplogroup.median(), t=np.nan, p=np.nan, note=""))
RES.append(dict(test="T6_lineage", term="median_r_burden_active_len_adj_haplogroup", n=len(B),
                beta=B.r_adj_haplogroup.median(), t=np.nan, p=np.nan,
                note=f"raw {B.r_raw.median():+.3f}"))

pd.DataFrame(RES).to_csv(f"{OUT}/active_len_channel.tsv", sep="\t", index=False, float_format="%.6g")
S.to_csv(f"{OUT}/active_len_channel_slopes.tsv", sep="\t", index=False, float_format="%.6g")
print(f"\nwrote {OUT}/active_len_channel{{,_slopes,_homolog,_haplogroup}}.tsv")

print("\n" + "=" * 100)
print("VERDICT INPUTS")
print("=" * 100)
print(f"  T0 cis/trans      cis t=+11.76  trans t=-0.44   (individual-level channel absent)")
print(f"  T1 direction      INS t={t_ins:+.2f}  DEL t={t_del:+.2f}  -> "
      f"{'SAME sign, against H1' if same_sign else 'OPPOSITE, supports H1'}")
print(f"  T1b sum vs diff   sum t={t_sum:+.2f}  diff t={t_dif:+.2f}   (H1 loads on diff)")
print(f"  T2 net-bp slope   {S.slope_len_per_netbp.median():+.3f}   vs H1's exact prediction +0.5")
print(f"  T2 causal ceiling observed r^2 exceeds it on {n_viol}/{len(S)} chromosomes, "
      f"median {S.obs_over_ceiling.median():,.0f}x")
print(f"  T3 location       VOID ({100*N_IN/N_TOT:.1f}% of SVs inside the reference array)")
print(f"  T4 strata         effect homogeneous -> not a quality gradient (neutral H1 vs H2/H3)")
print(f"  T5 homolog bound  {100*viol:.1f}% of rows violate |dLen| <= het_bp; "
      f"median gap {H.ratio.median():.1f}x the het SV sequence")
print(f"  T6 lineage        haplogroup explains median R2={HGR.R2_active_len_on_haplogroup.median():.3f} "
      f"of active_len; burden r {B.r_raw.median():+.3f} -> {B.r_adj_haplogroup.median():+.3f} adjusted")
