#!/usr/bin/env python3
"""B2: build the covariate-matched null CpG set, and in doing so RUN the positivity check."""
import os
import sys
import glob
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

SEED = 20260830
rng = np.random.default_rng(SEED)

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
RES = PROJ + "/resources"
UNIV = UNIT + "/output/null_design/cpg_universe"
OUT = UNIT + "/output/null_design"

K_NULL = 20
MAX_SV_DIST = 8_000_000
EXCLUDE_BUFFER = 10_000
EXCLUDE_REL = 0.25

CALIPERS = [
    dict(name="tight",  dist_rel=0.20, meth_abs=0.05, nobs_sd=0.25, require_censat=True),
    dict(name="medium", dist_rel=0.35, meth_abs=0.10, nobs_sd=0.50, require_censat=True),
    dict(name="loose",  dist_rel=0.50, meth_abs=0.15, nobs_sd=1.00, require_censat=True),
]


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
tg = pd.read_csv(UNIT + "/output/targets/targets_cpg.tsv", sep="\t")
canon = pd.read_csv(OUT + "/b0_sv_units/b0_sv_unit_map.tsv", sep="\t")
print("discovery pairs " + str(len(pairs)) + ", target CpGs " + str(len(tg))
      + ", canonical SV units " + str(len(canon)))

sv = pairs.drop_duplicates("sv_id")[["sv_id", "source_chrom", "source_pos",
                                     "source_end", "source_type"]].copy()
sv["sv_mid"] = ((sv.source_pos + sv.source_end) // 2).astype(np.int64)
print("source SVs " + str(len(sv)) + " on " + str(sv.source_chrom.nunique()) + " chromosomes")

def load_censat(path, class_from_name=True):
    b = pd.read_csv(path, sep="\t", header=None,
                    names=["chrom", "start", "end", "name", "score", "strand",
                           "ts", "te", "rgb"])
    if class_from_name:
        b["cls"] = b.name.astype(str).str.split("_").str[0]
    return b[["chrom", "start", "end", "name", "cls"]]


cen = load_censat(RES + "/chm13v2.0_censat_v2.1.bed")
act = load_censat(RES + "/chm13v2.0_censat_v2.1.active.bed")
print("censat intervals " + str(len(cen)) + ", classes "
      + str(sorted(cen.cls.unique())[:10]))
print("active HOR arrays " + str(len(act)))


def annotate(df, chrom_col, pos_col):
    out_cls = np.full(len(df), "none", dtype=object)
    out_dist = np.full(len(df), np.nan)
    for c, sub in df.groupby(chrom_col):
        idx = sub.index.values
        p = sub[pos_col].values.astype(np.int64)
        cc = cen[cen.chrom == c]
        if len(cc):
            s = cc.start.values.astype(np.int64)
            e = cc.end.values.astype(np.int64)
            order = np.argsort(s)
            s, e = s[order], e[order]
            cls = cc.cls.values[order]
            j = np.searchsorted(s, p, side="right") - 1
            ok = (j >= 0) & (p < e[np.clip(j, 0, len(e) - 1)])
            v = np.full(len(p), "none", dtype=object)
            v[ok] = cls[j[ok]]
            out_cls[np.searchsorted(df.index.values, idx)] = v
        aa = act[act.chrom == c]
        if len(aa):
            edges = np.sort(np.concatenate([aa.start.values, aa.end.values]).astype(np.int64))
            k = np.searchsorted(edges, p)
            lo = edges[np.clip(k - 1, 0, len(edges) - 1)]
            hi = edges[np.clip(k, 0, len(edges) - 1)]
            out_dist[np.searchsorted(df.index.values, idx)] = np.minimum(
                np.abs(p - lo), np.abs(p - hi))
    return out_cls, out_dist


files = sorted(glob.glob(UNIV + "/*/*.cpg.tsv.gz"))
tasklist = [l.strip() for l in open(UNIV + "/_tasklist.txt") if l.strip()]
print("\nuniverse chunk files present: " + str(len(files)) + " / " + str(len(tasklist)))
if len(files) < len(tasklist):
    print("INCOMPLETE -- the array has not finished. Aborting rather than building a")
    print("null from a partial universe, which would silently bias the matched pool.")
    sys.exit(2)

keep = []
for f in files:
    d = pd.read_csv(f, sep="\t", header=None,
                    names=["chrom", "start", "end", "feature_id", "n_observed",
                           "n_valid", "mean_meth", "sd_meth"])
    d = d[d.mean_meth.notna()]
    sv_c = sv[sv.source_chrom == d.chrom.iloc[0]] if len(d) else sv.iloc[:0]
    if not len(sv_c):
        continue
    pos = d.start.values.astype(np.int64)
    mind = np.full(len(d), np.inf)
    for _, r in sv_c.iterrows():
        mind = np.minimum(mind, np.abs(pos - r.sv_mid))
    d = d[mind <= MAX_SV_DIST]
    if len(d):
        keep.append(d)
    print("  " + os.path.basename(f) + "  kept " + str(len(d)))
U = pd.concat(keep, ignore_index=True)
print("universe within " + str(MAX_SV_DIST // 1_000_000) + " Mb of a source SV: " + str(len(U)))

U = U.sort_values(["chrom", "start"]).reset_index(drop=True)
U["cls"], U["act_dist"] = annotate(U, "chrom", "start")
print("universe cenSat classes: " + str(dict(pd.Series(U.cls).value_counts().head(8))))

SD_NVALID = float(U.n_valid.std())
SD_METH = float(U.mean_meth.std())
print("empirical scales: sd(n_valid)=" + format(SD_NVALID, ".3f")
      + "  sd(mean_meth)=" + format(SD_METH, ".4f"))
assert SD_NVALID > 0 and SD_METH > 0

disc = tg.merge(U[["chrom", "start", "feature_id", "n_observed", "n_valid", "mean_meth",
                   "sd_meth", "cls", "act_dist"]],
                left_on=["chm13_chrom", "chm13_pos"], right_on=["chrom", "start"],
                how="left", suffixes=("", "_u"))
n_missing = int(disc.mean_meth.isna().sum())
print("\ndiscovery targets located in the universe: "
      + str(len(disc) - n_missing) + " / " + str(len(disc)))
assert n_missing < len(disc), "no discovery target found in the universe -- key mismatch"
if n_missing:
    print("  " + str(n_missing) + " not found (outside the +/-" + str(MAX_SV_DIST // 10**6)
          + " Mb window or a different matrix chunk); those pairs cannot be matched")

excl = pairs.merge(tg, on="target_id")[["sv_id", "chm13_chrom", "chm13_pos"]].drop_duplicates()

disc_by_tid = disc.set_index("target_id")

print("\nprecomputing candidate pools per source SV ...")
POOLS = {}
for _, s in sv.iterrows():
    cand = U[U.chrom == s.source_chrom]
    if not len(cand):
        POOLS[s.sv_id] = cand.iloc[:0].assign(d=np.array([], dtype=np.int64))
        continue
    cand = cand.copy()
    cand["d"] = np.abs(cand.start.values.astype(np.int64) - int(s.sv_mid))
    ex = excl[(excl.sv_id == s.sv_id) & (excl.chm13_chrom == s.source_chrom)]
    if len(ex):
        e = np.sort(ex.chm13_pos.values.astype(np.int64))
        p = cand.start.values.astype(np.int64)
        k = np.searchsorted(e, p)
        lo = e[np.clip(k - 1, 0, len(e) - 1)]
        hi = e[np.clip(k, 0, len(e) - 1)]
        cand["excl_dist"] = np.minimum(np.abs(p - lo), np.abs(p - hi))
    else:
        cand["excl_dist"] = np.inf
    cand = cand[~cand.feature_id.isin(set(tg.feature_id))]
    POOLS[s.sv_id] = cand.reset_index(drop=True)
    print("  " + s.sv_id + "  candidates=" + str(len(cand)))

SVROW = sv.set_index("sv_id")
rows, diag = [], []
for _, pr in pairs.iterrows():
    tid = pr.target_id
    if tid not in disc_by_tid.index or pd.isna(disc_by_tid.loc[tid, "mean_meth"]):
        diag.append(dict(pair_id=pr.pair_id, sv_id=pr.sv_id, target_id=tid,
                         caliper="NONE", pool=0, reason="target_not_in_universe"))
        continue
    t = disc_by_tid.loc[tid]
    s = SVROW.loc[pr.sv_id]
    d_i = abs(int(t.start) - int(s.sv_mid))
    if d_i > MAX_SV_DIST:
        diag.append(dict(pair_id=pr.pair_id, sv_id=pr.sv_id, target_id=tid,
                         caliper="NONE", pool=0, reason="distance_beyond_window"))
        continue

    cand = POOLS[pr.sv_id]
    chosen, used, last_pool = None, None, 0
    for cal in CALIPERS:
        buf = min(EXCLUDE_BUFFER, EXCLUDE_REL * max(d_i, 1))
        m = (np.abs(cand.d.values - d_i) <= cal["dist_rel"] * max(d_i, 1)) \
            & (cand.excl_dist.values > buf) \
            & (np.abs(cand.mean_meth.values - float(t.mean_meth)) <= cal["meth_abs"]) \
            & (np.abs(cand.n_valid.values - float(t.n_valid)) <= cal["nobs_sd"] * SD_NVALID)
        if cal["require_censat"]:
            m = m & (cand.cls.values == t.cls)
        n_pool = int(m.sum())
        last_pool = max(last_pool, n_pool)
        if n_pool >= K_NULL:
            chosen, used = cand[m], cal["name"]
            break
    if chosen is None:
        diag.append(dict(pair_id=pr.pair_id, sv_id=pr.sv_id, target_id=tid,
                         caliper="FAILED", pool=last_pool,
                         reason="pool_below_K_at_loosest_caliper"))
        continue

    z = (np.abs(np.log10(chosen.d.clip(lower=1)) - np.log10(max(d_i, 1)))
         + np.abs(chosen.mean_meth - float(t.mean_meth)) / SD_METH
         + np.abs(chosen.n_valid - float(t.n_valid)) / SD_NVALID)
    take = chosen.assign(_z=z.values).sort_values(["_z", "start"]).head(K_NULL)
    diag.append(dict(pair_id=pr.pair_id, sv_id=pr.sv_id, target_id=tid,
                     caliper=used, pool=int(len(chosen)), reason="ok",
                     excl_buffer_bp=int(min(EXCLUDE_BUFFER, EXCLUDE_REL * max(d_i, 1)))))
    for _, c in take.iterrows():
        rows.append(dict(
            null_pair_id=pr.pair_id + "__n" + str(len(rows)).zfill(6),
            matched_to_pair_id=pr.pair_id, sv_id=pr.sv_id,
            null_chrom=c.chrom, null_pos=int(c.start), null_feature_id=c.feature_id,
            null_dist_to_sv=int(c.d), disc_dist_to_sv=int(d_i),
            null_mean_meth=float(c.mean_meth), disc_mean_meth=float(t.mean_meth),
            null_n_valid=int(c.n_valid), disc_n_valid=int(t.n_valid),
            null_cls=c.cls, disc_cls=t.cls, caliper=used))

nullp = pd.DataFrame(rows)
dg = pd.DataFrame(diag)

print("\n=== POSITIVITY / OVERLAP ===")
print(dg.caliper.value_counts().to_string())
ok = dg[dg.caliper.isin([c["name"] for c in CALIPERS])]
print("\npairs with an adequate matched pool: " + str(len(ok)) + " / " + str(len(dg))
      + "  (" + str(round(100.0 * len(ok) / max(len(dg), 1), 1)) + "%)")
print("pairs OFF the null's support        : " + str(len(dg) - len(ok)))
print("\nreason breakdown for failures:")
print(dg[~dg.caliper.isin([c["name"] for c in CALIPERS])].reason.value_counts().to_string())

if len(ok):
    st = pd.read_csv(UNIT + "/output/targets/prereg_strata.tsv", sep="\t")
    okj = ok.merge(st[["pair_id", "stratum", "recoverable_pilot8", "priority84"]],
                   on="pair_id", how="left")
    print("\nadequate-pool rate by pre-registered stratum:")
    z = (dg.merge(st[["pair_id", "stratum"]], on="pair_id", how="left")
         .assign(ok=lambda x: x.caliper.isin([c["name"] for c in CALIPERS]))
         .groupby("stratum").ok.agg(["size", "sum", "mean"]))
    print(z.to_string())

print("\nnull pairs built: " + str(len(nullp))
      + "   distinct null CpGs: " + str(nullp.null_feature_id.nunique() if len(nullp) else 0))
write_atomic(nullp, OUT + "/b2_matched_null/b2_null_pairs.tsv")
write_atomic(dg, OUT + "/b2_matched_null/b2_positivity_diagnostics.tsv")
if len(nullp):
    nt = (nullp[["null_chrom", "null_pos", "null_feature_id"]].drop_duplicates()
          .sort_values(["null_chrom", "null_pos"]).reset_index(drop=True))
    nt["null_target_id"] = "n" + pd.Series(np.arange(len(nt))).astype(str).str.zfill(6)
    write_atomic(nt[["null_target_id", "null_feature_id", "null_chrom", "null_pos"]],
                 OUT + "/b2_matched_null/b2_null_targets_cpg.tsv")
    print("\nthese null CpGs are lifted in the SAME v2 pass as the discovery targets,")
    print("so any lift artefact is shared and therefore measurable rather than confounding.")
