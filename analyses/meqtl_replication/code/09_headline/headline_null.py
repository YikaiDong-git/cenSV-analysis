#!/usr/bin/env python3
import glob
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
RES = PROJ + "/resources"
UNIV = UNIT + "/output/null_design/cpg_universe"
OUT = UNIT + "/output/headline_svs"

K_NULL = 5
MAX_SV_DIST = 100_000_000
EXCLUDE_BUFFER = 10_000
EXCLUDE_REL = 0.25
CALIPERS = [
    dict(name="tight",  dist_rel=0.20, meth_abs=0.05, nobs_sd=0.25, require_censat=True),
    dict(name="medium", dist_rel=0.35, meth_abs=0.10, nobs_sd=0.50, require_censat=True),
    dict(name="loose",  dist_rel=0.50, meth_abs=0.15, nobs_sd=1.00, require_censat=True),
]
SVS = [("chr12_DEL333", "chr12", 37183837, 37184170),
       ("chr8_DEL113", "chr8", 44300221, 44300334)]


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


D = []
for f in sorted(glob.glob(OUT + "/inputs/*_cpgs.tsv")):
    d = pd.read_csv(f, sep="\t")
    d["chrom"] = d.feature_id.astype(str).str.split(":").str[0]
    d["pos"] = d.feature_id.astype(str).str.split(":").str[1].str.split("-").str[0].astype(np.int64)
    D.append(d)
D = pd.concat(D, ignore_index=True).drop_duplicates(["sv_name", "feature_id"])
print("headline discovery CpGs: " + str(len(D)) + "  " + str(dict(D.sv_name.value_counts())))

var = pd.read_csv(SUB + "/geno/full_variants.tsv", sep="\t", engine="python")
rows = []
for name, ch, pos, end in SVS:
    hit = var[(var.chrom == ch) & (var.pos == pos) & (var.end == end) & (var.svtype == "DEL")]
    for i, r in hit.iterrows():
        rows.append(dict(sv_name=name, sv_locus=r.sv_locus, geno_group=int(r.geno_group),
                         row_index=int(i), af=float(r.af), strategy=r.strategy,
                         status=r.status))
B = pd.DataFrame(rows).drop_duplicates(["sv_name", "geno_group"])
print("\ngenotype branches carried forward (all, per the run-all-branches decision):")
print(B[["sv_name", "geno_group", "af", "strategy", "status"]].to_string(index=False))
write_atomic(B, OUT + "/inputs/branches.tsv")

def load_bed(p):
    b = pd.read_csv(p, sep="\t", header=None,
                    names=["chrom", "start", "end", "name", "score", "strand", "ts", "te", "rgb"])
    b["cls"] = b.name.astype(str).str.split("_").str[0]
    return b[["chrom", "start", "end", "cls"]]


cen = load_bed(RES + "/chm13v2.0_censat_v2.1.bed")


def annotate(df):
    out = np.full(len(df), "none", dtype=object)
    for c, sub in df.groupby("chrom"):
        cc = cen[cen.chrom == c]
        if not len(cc):
            continue
        s = cc.start.values.astype(np.int64)
        e = cc.end.values.astype(np.int64)
        o = np.argsort(s)
        s, e, cl = s[o], e[o], cc.cls.values[o]
        p = sub.pos.values.astype(np.int64)
        j = np.searchsorted(s, p, side="right") - 1
        ok = (j >= 0) & (p < e[np.clip(j, 0, len(e) - 1)])
        v = np.full(len(p), "none", dtype=object)
        v[ok] = cl[j[ok]]
        out[sub.index.values] = v
    return out


keep = []
for f in sorted(glob.glob(UNIV + "/*/*.cpg.tsv.gz")):
    ch = os.path.basename(os.path.dirname(f))
    if ch not in {"chr12", "chr8"}:
        continue
    d = pd.read_csv(f, sep="\t", header=None,
                    names=["chrom", "start", "end", "feature_id", "n_observed",
                           "n_valid", "mean_meth", "sd_meth"])
    d = d[d.mean_meth.notna()]
    keep.append(d)
U = pd.concat(keep, ignore_index=True).reset_index(drop=True)
U = U.rename(columns={"start": "pos"})
U["cls"] = annotate(U)
SD_NVALID = float(U.n_valid.std())
SD_METH = float(U.mean_meth.std())
print("\nuniverse (chr8+chr12): " + str(len(U)) + " CpGs   sd(n_valid)="
      + format(SD_NVALID, ".3f") + "  sd(mean_meth)=" + format(SD_METH, ".4f"))

Uk = U.set_index("feature_id")
disc = D.merge(U[["feature_id", "n_valid", "mean_meth", "cls"]], on="feature_id", how="left")
print("discovery CpGs located in the universe: "
      + str(int(disc.mean_meth.notna().sum())) + " / " + str(len(disc)))

rows, diag = [], []
for name, ch, pos, end in SVS:
    mid = (pos + end) // 2
    cand = U[U.chrom == ch].copy()
    cand["d"] = np.abs(cand.pos.values.astype(np.int64) - mid)
    ex = np.sort(disc[disc.sv_name == name].pos.values.astype(np.int64))
    p = cand.pos.values.astype(np.int64)
    k = np.searchsorted(ex, p)
    lo = ex[np.clip(k - 1, 0, len(ex) - 1)]
    hi = ex[np.clip(k, 0, len(ex) - 1)]
    cand["excl_dist"] = np.minimum(np.abs(p - lo), np.abs(p - hi))
    cand = cand[~cand.feature_id.isin(set(disc.feature_id))].reset_index(drop=True)
    print("\n" + name + ": candidate pool " + str(len(cand)))

    sub = disc[(disc.sv_name == name) & disc.mean_meth.notna()]
    for t in sub.itertuples(index=False):
        d_i = abs(int(t.pos) - mid)
        if d_i > MAX_SV_DIST:
            diag.append(dict(sv_name=name, feature_id=t.feature_id, caliper="NONE",
                             pool=0, reason="beyond_window"))
            continue
        chosen, used, last = None, None, 0
        for cal in CALIPERS:
            buf = min(EXCLUDE_BUFFER, EXCLUDE_REL * max(d_i, 1))
            m = (np.abs(cand.d.values - d_i) <= cal["dist_rel"] * max(d_i, 1)) \
                & (cand.excl_dist.values > buf) \
                & (np.abs(cand.mean_meth.values - float(t.mean_meth)) <= cal["meth_abs"]) \
                & (np.abs(cand.n_valid.values - float(t.n_valid)) <= cal["nobs_sd"] * SD_NVALID)
            if cal["require_censat"]:
                m = m & (cand.cls.values == t.cls)
            n = int(m.sum())
            last = max(last, n)
            if n >= K_NULL:
                chosen, used = cand[m], cal["name"]
                break
        if chosen is None:
            diag.append(dict(sv_name=name, feature_id=t.feature_id, caliper="FAILED",
                             pool=last, reason="pool_below_K"))
            continue
        z = (np.abs(np.log10(chosen.d.clip(lower=1)) - np.log10(max(d_i, 1)))
             + np.abs(chosen.mean_meth - float(t.mean_meth)) / SD_METH
             + np.abs(chosen.n_valid - float(t.n_valid)) / SD_NVALID)
        take = chosen.assign(_z=z.values).sort_values(["_z", "pos"]).head(K_NULL)
        diag.append(dict(sv_name=name, feature_id=t.feature_id, caliper=used,
                         pool=int(len(chosen)), reason="ok"))
        for _, c in take.iterrows():
            rows.append(dict(sv_name=name, matched_to=t.feature_id,
                             null_feature_id=c.feature_id, chrom=c.chrom,
                             pos=int(c.pos), d=int(c.d), disc_d=int(d_i),
                             cls=c.cls, mean_meth=float(c.mean_meth),
                             n_valid=int(c.n_valid), caliper=used))
N = pd.DataFrame(rows)
G = pd.DataFrame(diag)
print("\n=== positivity for the headline SVs ===")
print(G.groupby(["sv_name", "caliper"]).size().to_string())
ok = G[G.caliper.isin([c["name"] for c in CALIPERS])]
print("\nCpGs with an adequate matched pool: " + str(len(ok)) + " / " + str(len(G)))
print(G.assign(ok=G.caliper.isin([c["name"] for c in CALIPERS]))
      .groupby("sv_name").ok.agg(["size", "sum", "mean"]).round(4).to_string())
write_atomic(N, OUT + "/null/headline_null_pairs.tsv")
write_atomic(G, OUT + "/null/headline_positivity.tsv")

d1 = disc[["feature_id", "chrom", "pos", "sv_name"]].copy()
d1["target_class"] = "discovery"
n1 = N[["null_feature_id", "chrom", "pos", "sv_name"]].rename(
    columns={"null_feature_id": "feature_id"}).drop_duplicates("feature_id")
n1["target_class"] = "matched_null"
T = pd.concat([d1, n1], ignore_index=True).drop_duplicates("feature_id")
print("\ncombined headline targets: " + str(len(T))
      + "   " + str(dict(T.target_class.value_counts())))
write_atomic(T, OUT + "/null/headline_targets_raw.tsv")
print("\nnext: normalise to the plus-strand C and assert CG, then lift.")
