import numpy as np, pandas as pd, struct
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = ROOT + "/analyses/cdr_qtl_full_panel"
CAN = ROOT + "/analyses/cdr_qtl"
RD = ROOT + "/analyses/cdr_phenotypes"
GRM = ROOT + "/subproject_methylation/results/greml/pericen_1m_all"

def rint(x):
    x = np.asarray(x, float); ok = ~np.isnan(x); o = np.full(x.shape, np.nan)
    from scipy.stats import rankdata, norm
    o[ok] = norm.ppf((rankdata(x[ok]) - .5) / ok.sum()); return o

def read_grm(p):
    ids = [l.split()[1] for l in open(p + ".grm.id")]
    n = len(ids); K = np.zeros((n, n))
    raw = open(p + ".grm.bin", "rb").read()
    v = struct.unpack("f" * (len(raw) // 4), raw); k = 0
    for i in range(n):
        for j in range(i + 1):
            K[i, j] = K[j, i] = v[k]; k += 1
    return K, ids

def slope(y, d, COV):
    ok = ~np.isnan(y) & ~np.isnan(d)
    y, d, COV = y[ok], d[ok], COV[ok]
    X = np.column_stack([np.ones(len(y)), d, COV])
    b = np.linalg.pinv(X) @ y
    r = y - X @ b
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X) * (r**2).sum() / (len(y) - X.shape[1])))
    t = b[1] / se[1]
    return float(t), float(2 * (1 - stats.t.cdf(abs(t), len(y) - X.shape[1]))), int(len(y))

v1 = pd.read_csv(CAN + "/output/per_hap_chrom_phenotypes.tsv", sep="\t", engine="python")
v2 = pd.read_csv(RD + "/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv", sep="\t", engine="python")
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], index=gids)
G = np.load(WU + "/geno/full_genotypes.npz"); G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
samp = np.array(open(WU + "/geno/full_samples.txt").read().split())
loc2row = {l: i for i, l in enumerate(vt.sv_locus)}

SVS = [("DEL-678@44982379",  "chr9_sv20_DEL_mon_best6_alt_edge500"),
       ("DEL-7310@44984820", "chr9_sv25_DEL_mon_best6_alt_edge2000")]

ARMS = [
    ("v1 as shipped (46.44 Mb row IN)",      v1, "active_len"),
    ("v1 minus HG00146 only",                v1[v1["sample"] != "HG00146"], "active_len"),
    ("v2 corrected (SUM, artifact dropped)", v2, "active_len_sum"),
]
print("chr9 array-length test: does the manuscript's t = -4.65 survive?\n")
print("%-38s %-20s %8s %12s %6s" % ("arm", "SV", "t(len)", "p(len)", "n"))
for name, df, col in ARMS:
    s = df[df.chm13_chrom == "chr9"]
    b = s.groupby("sample").agg(aa_l=(col, "mean"))
    bc = s.groupby("sample")["basecaller"].agg(
        lambda v: 1.0 if v.astype(str).str.contains("1041").any() else 0.0)
    b = b.join(bc.rename("bc"))
    for label, locus in SVS:
        dos = pd.Series(G[loc2row[locus], :], index=samp).replace(-1, np.nan)
        keep = sorted(set(b.dropna().index) & set(dos.dropna().index) & set(PCs.index))
        COV = np.column_stack([PCs.loc[keep].values, b.loc[keep, "bc"].values[:, None]])
        t, p, n = slope(rint(b.loc[keep, "aa_l"].values), dos.loc[keep].values, COV)
        print("%-38s %-20s %+8.2f %12.2e %6d" % (name, label, t, p, n))
    print()

s9 = v1[v1.chm13_chrom == "chr9"]
print("chr9 active_len in v1: n=%d, median %.2f Mb, max %.2f Mb (%s)"
      % (len(s9), s9.active_len.median()/1e6, s9.active_len.max()/1e6,
         s9.loc[s9.active_len.idxmax(), "sample"]))
