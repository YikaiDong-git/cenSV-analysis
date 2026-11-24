#!/usr/bin/env python
"""Assemble the FULL-substrate per-sample genotype matrix: all 11,129 loci x 188 individuals."""
import pandas as pd, numpy as np, os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
SVD = f"{ROOT}/sv_extraction/decompose/{SUF}"
FRESH = f"{ROOT}/analyses/cdr_qtl/canonical_geno/decompose"
EN = f"{ROOT}/results/elastic_net/decompose/cosine_nj"
ADM = f"{ROOT}/analyses/sv_catalog/output/ref_ukbb/admit/admit_decompose_cv01.tsv"
OUT = f"{ROOT}/analyses/cdr_qtl_full_panel/geno"
os.makedirs(OUT, exist_ok=True)

CANON_SAMPLES = open(f"{FRESH}/chr8/chr8_sv_samples.txt").read().split()
sidx = {s: i for i, s in enumerate(CANON_SAMPLES)}
NS = len(CANON_SAMPLES)
assert NS == 188, NS
print(f"canonical sample order: {NS}")


def memset(s):
    return frozenset(str(s).split(";"))


def gcp(s):
    d = {p.split(":")[0]: int(p.split(":")[1]) for p in str(s).split("/") if ":" in p}
    return (d.get("0", 0), d.get("1", 0), d.get("2", 0))


adf = pd.read_csv(ADM, sep="\t", dtype=str)
adf["cv_r2"] = adf["cv_r2"].astype(float)
adf["chrom"] = "chr" + adf["primary_chrom"].astype(str)
best = adf.sort_values(["cv_r2", "sv_locus"], ascending=[False, True]).drop_duplicates("sv_locus")
loci = best[["sv_locus", "chrom", "sv_id", "cv_r2", "strategy", "arm", "status"]].copy()
loci["memset"] = loci.sv_id.map(memset)
print(f"loci (all statuses): {len(loci)}  across {loci.chrom.nunique()} chroms")
print(loci.status.value_counts().to_string())


def load_truth(ch):
    for s in ["midpoint1000", "midpoint500", "edge500", "edge1000", "edge2000",
              "fksize100k", "fksize200k", "fksize500k"]:
        p = f"{EN}/{s}_{SUF}/{ch}/sv_metrics.tsv"
        if os.path.exists(p):
            m = pd.read_csv(p, sep="\t", usecols=["svpos", "svend", "svtype", "svlen", "svfreq", "gt_counts", "sv_id"])
            m["memset"] = m.sv_id.map(memset)
            m["gc"] = m.gt_counts.map(gcp)
            return {r.memset: (r.gc, int(r.svpos),
                               (int(r.svend) if not pd.isna(r.svend) else int(r.svpos)),
                               str(r.svtype), int(r.svlen)) for r in m.itertuples()}
    return {}


def load_path_vectors(variants_tsv, npz_path, samples_txt):
    if not (os.path.exists(variants_tsv) and os.path.exists(npz_path) and os.path.exists(samples_txt)):
        return {}
    v = pd.read_csv(variants_tsv, sep="\t")
    G = np.load(npz_path)
    G = G[list(G.keys())[0]]
    if G.shape[0] != len(v) and G.shape[1] == len(v):
        G = G.T
    samp = np.array(open(samples_txt).read().split())
    col_for = {s: j for j, s in enumerate(samp)}
    take = np.array([col_for.get(s, -1) for s in CANON_SAMPLES])
    present = take >= 0
    out = {}
    ids = v["id"].astype(str).values
    for i in range(len(v)):
        raw = np.rint(G[i, :].astype(float)).clip(0, 2).astype(int)
        vec = np.full(NS, -1, dtype=np.int64)
        vec[present] = raw[take[present]]
        pv = vec[vec >= 0]
        gc = (int((pv == 0).sum()), int((pv == 1).sum()), int((pv == 2).sum()))
        out.setdefault(memset(ids[i]), []).append((vec, gc))
    return out


rows = []
MAT = np.full((len(loci), NS), -1, dtype=np.int64)
ri = 0
counts = {"path1": 0, "path2": 0, "drift_approx": 0, "unrecovered": 0}
loci_by_chr = {ch: g for ch, g in loci.groupby("chrom")}
for ch in sorted(loci_by_chr, key=lambda c: (len(c), c)):
    truth = load_truth(ch)
    p1 = load_path_vectors(f"{SVD}/{ch}/{ch}_sv_variants.tsv", f"{SVD}/{ch}/{ch}_sv_genotypes.npz", f"{SVD}/{ch}/{ch}_sv_samples.txt")
    p2 = load_path_vectors(f"{FRESH}/{ch}/{ch}_sv_variants.tsv", f"{FRESH}/{ch}/{ch}_sv_genotypes.npz", f"{FRESH}/{ch}/{ch}_sv_samples.txt")
    nch = 0
    for r in loci_by_chr[ch].itertuples():
        t = truth.get(r.memset)
        gc_true = t[0] if t else None
        pos = t[1] if t else -1
        end = t[2] if t else -1
        svt = t[3] if t else ""
        svl = t[4] if t else -1
        chosen = None; path = None; fid = None; delta = None
        for pname, pdict in (("path1", p1), ("path2", p2)):
            for vec, gc in pdict.get(r.memset, []):
                if gc_true is not None and gc == gc_true:
                    chosen = vec; path = pname; fid = "exact"; delta = 0
                    break
            if chosen is not None:
                break
        if chosen is None:
            bestd = 1e9; bestv = None; bestp = None
            for pname, pdict in (("path1", p1), ("path2", p2)):
                for vec, gc in pdict.get(r.memset, []):
                    d = sum(abs(a - b) for a, b in zip(gc, gc_true)) if gc_true else 1e9
                    if d < bestd:
                        bestd = d; bestv = vec; bestp = pname
            if bestv is not None:
                chosen = bestv; path = bestp; fid = "drift_approx"; delta = int(bestd)
                counts["drift_approx"] += 1
            else:
                counts["unrecovered"] += 1
        else:
            counts[path] += 1
        if chosen is not None:
            MAT[ri, :] = chosen
            pv = chosen[chosen >= 0]
            ac = int(pv[pv > 0].sum()); an = 2 * len(pv); af = ac / an if an else 0
            rows.append(dict(sv_locus=r.sv_locus, chrom=ch, pos=pos, end=end, svtype=svt, svlen=svl,
                             ac=ac, an=an, af=round(af, 6), maf=round(min(af, 1 - af), 6),
                             n_miss=int((chosen < 0).sum()), sv_id=r.sv_id, cv_r2=round(r.cv_r2, 4),
                             strategy=r.strategy, arm=r.arm, status=r.status,
                             path=path, fidelity=fid, gt_delta=delta))
            ri += 1; nch += 1
    print(f"  {ch:6s} assembled {nch:5d}/{len(loci_by_chr[ch]):5d}")

MAT = MAT[:ri, :]
vt = pd.DataFrame(rows)

keys = [MAT[i, :].tobytes() for i in range(MAT.shape[0])]
seen = {}
grp = np.empty(len(keys), dtype=np.int64)
for i, k in enumerate(keys):
    if k not in seen:
        seen[k] = len(seen)
    grp[i] = seen[k]
vt["geno_group"] = grp
n_uniq = len(seen)

np.savez_compressed(f"{OUT}/full_genotypes.npz", genotypes=MAT.astype(np.int8))
open(f"{OUT}/full_samples.txt", "w").write("\n".join(CANON_SAMPLES) + "\n")
tmp = f"{OUT}/full_variants.tsv.tmp.{os.getpid()}"
vt.to_csv(tmp, sep="\t", index=False)
os.replace(tmp, f"{OUT}/full_variants.tsv")

print(f"\n=== assembled FULL matrix: {MAT.shape} (loci x samples) ===")
print("fidelity:"); print(vt.fidelity.value_counts().to_string())
print("path:"); print(vt.path.value_counts().to_string())
print(f"unrecovered (no member match in either decompose): {counts['unrecovered']}")
print(f"\n=== framings ===")
print(f"  row-level loci        : {len(vt)}")
print(f"  unique genotypes      : {n_uniq}  (bit-identical 188-vectors collapsed)")
adm_v = vt[vt.status == "ADMIT"]
print(f"  ADMIT subset loci     : {len(adm_v)}   unique genotypes: {adm_v.geno_group.nunique()}")
print(f"  gain from full substrate: +{len(vt)-len(adm_v)} loci, +{n_uniq-adm_v.geno_group.nunique()} unique genotypes")
print("\nper-chrom (loci / unique genotypes), ADMIT vs ALL:")
tt = pd.DataFrame({"loci_admit": adm_v.groupby("chrom").size(),
                   "uniq_admit": adm_v.groupby("chrom").geno_group.nunique(),
                   "loci_all": vt.groupby("chrom").size(),
                   "uniq_all": vt.groupby("chrom").geno_group.nunique()}).fillna(0).astype(int)
tt = tt.sort_index(key=lambda s: [(len(x), x) for x in s])
print(tt.to_string())
print(f"TOTAL loci_admit={tt.loci_admit.sum()} uniq_admit={tt.uniq_admit.sum()} loci_all={tt.loci_all.sum()} uniq_all={tt.uniq_all.sum()}")

hl = vt[(vt.chrom == "chr8") & (vt.svtype == "DEL") & (vt.pos == 44300221) & (vt.svlen == 113)]
print("\nHEADLINE chr8 DEL-113 @44300221:")
print(hl[["sv_locus", "ac", "af", "cv_r2", "status", "path", "fidelity", "gt_delta", "geno_group"]].to_string(index=False)
      if len(hl) else "  NOT FOUND (check)")
print(f"\nwrote {OUT}/full_genotypes.npz + full_variants.tsv ({len(vt)} loci)")
