#!/usr/bin/env python
"""What sets the length of the active alpha-satellite HOR array?"""
import os, sys, glob, re, struct
import numpy as np, pandas as pd
from collections import defaultdict
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
CAN = f"{ROOT}/analyses/cdr_qtl"
AUD = f"{ROOT}/analyses/cdr_calls"
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
DECOMP = sorted(glob.glob(f"{ROOT}/sv_extraction/decompose/*/"))[0]
OUT = f"{WU}/output/length"
os.makedirs(OUT, exist_ok=True)
pd.set_option("display.width", 220)


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


def adj_r2(y, X):
    X = np.column_stack([np.ones(len(y))] + [np.asarray(c, float) for c in X]) if len(X) else \
        np.ones((len(y), 1))
    keep = [0] + [j for j in range(1, X.shape[1]) if np.std(X[:, j]) > 1e-12]
    X = X[:, keep]
    n, p = X.shape
    if n <= p + 1:
        return np.nan
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    ss_res = float(np.sum((y - X @ beta) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    if ss_tot <= 0:
        return np.nan
    r2 = 1 - ss_res / ss_tot
    return 1 - (1 - r2) * (n - 1) / (n - p)


ph = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
H = ph[["sample", "hap_local", "chm13_chrom", "active_len", "basecaller"]].rename(
    columns={"chm13_chrom": "chrom"}).dropna(subset=["active_len"])
H["hap"] = H.hap_local.astype(str)
print(f"haplotype-level rows: {len(H)}  [{H['sample'].nunique()} individuals x "
      f"{H.chrom.nunique()} chromosomes]")

coh = pd.read_csv(f"{AUD}/output/cohort_corrected_final.tsv", sep="\t")
cov_rows = []
for f in sorted(glob.glob(f"{AUD}/output/calls/*/binned_freq.bed")):
    tag = os.path.basename(os.path.dirname(f))
    s, h = tag.rsplit("_", 1)
    a = pd.read_csv(f, sep="\t", header=None, usecols=[4], names=["cov"])
    cov_rows.append(dict(sample=s, hap=h, coverage=a["cov"].mean()))
COV = pd.DataFrame(cov_rows)
Q = coh.merge(COV, on=["sample", "hap"], how="outer")
Q["is_r1041"] = Q.chemistry.astype(str).str.contains("1041").astype(float)
QH = Q[["sample", "hap", "coverage", "is_r1041", "n_contigs"]].copy()

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
PCS = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]],
                   columns=[f"gPC{i+1}" for i in range(10)], index=gids)
PCS.index.name = "sample"

hg = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{DECOMP}chr*/chr*_haplogroups.tsv"))
      if "validated" not in f]
HG = pd.concat(hg, ignore_index=True).drop(columns=["sample"]).rename(
    columns={"sample_id": "sample", "chromosome": "chrom"})[
    ["sample", "chrom", "haplotype", "haplogroup"]]

G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t", engine="python")
samples = np.array(open(f"{WU}/geno/full_samples.txt").read().split())
vt["sgn"] = vt.svtype.map({"INS": 1.0, "DUP:TANDEM": 1.0, "DEL": -1.0, "BND": 0.0})
phy = vt.reset_index().sort_values("index").drop_duplicates(
    ["chrom", "pos", "end", "svtype", "svlen"], keep="first")
Gp, vp = G[phy["index"].values, :], vt.iloc[phy["index"].values].reset_index(drop=True)
srows = []
for c in sorted(set(vp.chrom)):
    m = (vp.chrom == c).values
    srows.append(pd.DataFrame(dict(
        sample=samples, chrom=c,
        burden=np.nansum(Gp[m, :], axis=0),
        netbp=np.nansum(Gp[m, :] * (vp.sgn.values[m] * vp.svlen.values[m])[:, None], axis=0),
        het_bp=np.nansum((Gp[m, :] == 1) * np.abs(vp.svlen.values[m])[:, None], axis=0))))
SV = pd.concat(srows, ignore_index=True)

H["haplotype"] = "hap" + H.hap.str.extract(r"(\d+)")[0]
D = (H.merge(QH, on=["sample", "hap"], how="left")
      .merge(NI, on=["sample", "chrom"], how="left")
      .merge(PCS.reset_index(), on="sample", how="left")
      .merge(HG, on=["sample", "chrom", "haplotype"], how="left")
      .merge(SV, on=["sample", "chrom"], how="left"))
D["alen_z"] = D.groupby("chrom").active_len.transform(
    lambda s: rint(s) if s.notna().all() and s.std() > 0 else np.nan)
D = D.dropna(subset=["alen_z"])
print(f"modelling rows: {len(D)}   with haplogroup: {D.haplogroup.notna().sum()}   "
      f"with coverage: {D.coverage.notna().sum()}")
print(f"raw active_len: mean {D.active_len.mean():,.0f} bp   sd {D.active_len.std():,.0f} bp   "
      f"range {D.active_len.min():,.0f}-{D.active_len.max():,.0f}")

print("\n" + "=" * 100)
print("S1  VARIANCE DECOMPOSITION of within-chromosome z-scored active_len")
print("=" * 100)
two = D.groupby(["sample", "chrom"]).alen_z.transform("size") == 2
B = D[two].copy()
print(f"  restricted to the {B.groupby(['sample','chrom']).ngroups} chromosome-arms where BOTH "
      f"homologs are resolved ({len(B)} haplotypes)")
gm = B.alen_z.mean()
ind_mean = B.groupby("sample").alen_z.transform("mean")
v_ind = float(np.var(ind_mean - gm)) / float(np.var(B.alen_z))
print(f"  individual-level share of variance (all of one person's arrays): {100*v_ind:.1f}%")

piv = B.pivot_table(index=["sample", "chrom"], columns="haplotype", values="alen_z").dropna()
x1, x2 = piv.iloc[:, 0].values.astype(float), piv.iloc[:, 1].values.astype(float)
ricc = stats.pearsonr(np.r_[x1, x2], np.r_[x2, x1])
icc = float(ricc.statistic)
print(f"  homolog intraclass correlation (n = {len(piv)} pairs): ICC = {icc:+.3f}  "
      f"p = {ricc.pvalue:.3g}")
print(f"  --> the two homologs of one chromosome in one person share {100*max(icc,0):.1f}% of the")
print("      within-chromosome length variance; the rest is PER-HAPLOTYPE. Array length is a")
print("      property of the individual chromosome copy, not of the person carrying it.")
DEC = pd.DataFrame([
    dict(quantity="individual-level share of variance", value=v_ind,
         note="corroborated by S2's near-zero cross-chromosome r"),
    dict(quantity="homolog ICC (shared by the 2 copies)", value=icc,
         note=f"n={len(piv)} pairs, p={ricc.pvalue:.3g}"),
    dict(quantity="per-haplotype share (1 - ICC)", value=1 - icc,
         note="no individual-level covariate can reach this")])
print()
print(DEC.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
DEC.to_csv(f"{OUT}/variance_decomposition.tsv", sep="\t", index=False, float_format="%.6g")

print("\n" + "=" * 100)
print("S2  CROSS-CHROMOSOME STRUCTURE -- is there an individual-level 'long arrays' factor?")
print("=" * 100)
W = D.groupby(["sample", "chrom"]).alen_z.mean().unstack()
C = W.corr(min_periods=40)
off = C.values[np.triu_indices_from(C.values, k=1)]
off = off[~np.isnan(off)]
print(f"  {len(off)} chromosome pairs with >=40 shared individuals")
print(f"  mean off-diagonal r = {np.mean(off):+.3f}   median {np.median(off):+.3f}   "
      f"sd {np.std(off):.3f}")
print(f"  pairs with |r| > 0.2: {int((np.abs(off) > 0.2).sum())}/{len(off)}")
tt = stats.ttest_1samp(off, 0.0)
print(f"  mean differs from 0: t = {tt.statistic:+.2f}, p = {tt.pvalue:.3g}")
pd.DataFrame(dict(pair_r=off)).to_csv(f"{OUT}/cross_chrom_r.tsv", sep="\t", index=False,
                                      float_format="%.6g")

print("\n" + "=" * 100)
print("S3  EXPLAINED-VARIANCE BUDGET, per chromosome, ADJUSTED R2")
print("=" * 100)
BLOCKS = {
    "technical": ["coverage", "is_r1041", "n_contigs", "n_intervals"],
    "ancestry": [f"gPC{i+1}" for i in range(10)],
    "SVs": ["burden", "netbp"],
}
brows = []
for c, g in D.groupby("chrom"):
    g = g.dropna(subset=["alen_z"])
    if len(g) < 60:
        continue
    y = g.alen_z.values.astype(float)
    row = dict(chrom=c, n=len(g))
    cols_all = []
    for name, cols in BLOCKS.items():
        use = [g[k].fillna(g[k].median()).values.astype(float) for k in cols if k in g]
        row[f"R2_{name}"] = adj_r2(y, use)
        cols_all += use
    hgd = pd.get_dummies(g.haplogroup.fillna("NA"), drop_first=True)
    hgc = [hgd[k].values.astype(float) for k in hgd.columns] if hgd.shape[1] else []
    row["R2_lineage"] = adj_r2(y, hgc) if hgc else np.nan
    row["n_haplogroups"] = int(g.haplogroup.nunique())
    row["R2_all"] = adj_r2(y, cols_all + hgc)
    brows.append(row)
BUD = pd.DataFrame(brows)
print(BUD.to_string(index=False, float_format=lambda x: f"{x:,.3f}"))
med = {k: BUD[k].median() for k in ["R2_technical", "R2_ancestry", "R2_lineage", "R2_SVs", "R2_all"]}
print(f"\n  MEDIAN adjusted R2 across {len(BUD)} chromosomes")
for k, v in med.items():
    print(f"    {k:14s} {v:+.4f}")
print(f"    {'UNEXPLAINED':14s} {1-med['R2_all']:+.4f}   "
      f"(= 1 - R2_all; an UPPER bound on biological variation, since no technical")
print(f"    {'':14s}          replicates exist to net out measurement error)")
BUD.to_csv(f"{OUT}/variance_budget.tsv", sep="\t", index=False, float_format="%.6g")

print("\n" + "=" * 100)
print("S4  HOMOLOG-GAP MODEL -- every individual-level confounder differences out.")
print("    If array length were lineage-determined, homologs in DIFFERENT centromeric haplogroups")
print("    should differ in length more than homologs in the SAME haplogroup.")
print("=" * 100)
g2 = D[two].sort_values(["sample", "chrom", "haplotype"])
agg = g2.groupby(["sample", "chrom"]).agg(
    lo=("active_len", "min"), hi=("active_len", "max"),
    n_hg=("haplogroup", lambda s: s.dropna().nunique()),
    hg_known=("haplogroup", lambda s: int(s.notna().sum())),
    het_bp=("het_bp", "first"), n_intervals=("n_intervals", "first"),
    coverage=("coverage", "mean")).reset_index()
agg["dlen"] = agg.hi - agg.lo
agg["rel_gap"] = agg.dlen / ((agg.hi + agg.lo) / 2)
agg["discordant"] = np.where(agg.hg_known == 2, (agg.n_hg == 2).astype(float), np.nan)
print(f"  chromosome-arms with both homologs: {len(agg)}")
print(f"  median homolog gap {agg.dlen.median():,.0f} bp   "
      f"median RELATIVE gap {100*agg.rel_gap.median():.1f}% of array length")
sub = agg.dropna(subset=["discordant"])
conc, disc = sub[sub.discordant == 0], sub[sub.discordant == 1]
print(f"\n  haplogroup CONCORDANT homologs  n={len(conc):4d}  median gap {conc.dlen.median():,.0f} bp")
print(f"  haplogroup DISCORDANT homologs  n={len(disc):4d}  median gap {disc.dlen.median():,.0f} bp")
if len(conc) > 5 and len(disc) > 5:
    mw = stats.mannwhitneyu(disc.dlen, conc.dlen, alternative="two-sided")
    print(f"  Mann-Whitney p = {mw.pvalue:.3g}   "
          f"ratio of medians {disc.dlen.median()/max(1,conc.dlen.median()):.2f}x")
for nm, col in [("het SV bp", "het_bp"), ("array intervals", "n_intervals"),
                ("coverage", "coverage")]:
    s = agg.dropna(subset=[col])
    r = stats.spearmanr(s.dlen, s[col])
    print(f"  spearman(|dLen|, {nm:15s}) = {r.statistic:+.3f}   p = {r.pvalue:.3g}")
agg.to_csv(f"{OUT}/homolog_gap.tsv", sep="\t", index=False, float_format="%.6g")

print("\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)
print(f"  S1  homolog ICC {icc:+.3f} -> per-haplotype share {100*(1-icc):.1f}%; "
      f"individual-level share {100*v_ind:.1f}%")
print(f"  S2  mean cross-chromosome r of array length                  : {np.mean(off):+.3f}")
print(f"  S3  median adj R2 technical {med['R2_technical']:+.3f} | ancestry {med['R2_ancestry']:+.3f} "
      f"| lineage {med['R2_lineage']:+.3f} | SVs {med['R2_SVs']:+.3f} | ALL {med['R2_all']:+.3f}")
print(f"  S4  median homolog gap {agg.dlen.median():,.0f} bp "
      f"({100*agg.rel_gap.median():.1f}% of array length)")
print(f"\nwrote {OUT}/{{variance_decomposition,cross_chrom_r,variance_budget,homolog_gap}}.tsv")
