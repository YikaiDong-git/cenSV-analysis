#!/usr/bin/env python3
import pandas as pd, numpy as np, re, os, glob
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
ROOT = CENSV_ROOT + ""
ADMIT = f"{ROOT}/analyses/sv_catalog/output/ref_ukbb/admit/admit_decompose_cv01.tsv"
DBASE = f"{ROOT}/sv_extraction/decompose/mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
PERCH = f"{ROOT}/analyses/cdr_qtl/results/cdr_qtl_corrected_1e6/perchrom"
EN    = f"{ROOT}/results/elastic_net/decompose/cosine_nj"
OUTD  = f"{ROOT}/analyses/cdr_qtl/results/canonical_admit"
os.makedirs(OUTD, exist_ok=True)
HIT_P = 1e-4

def parse_ids(s):
    if not isinstance(s, str): return set()
    return set(x for x in re.split(r"[;,]", s.replace("_B", ";_B")) if x.startswith("svim_asm."))

adf = pd.read_csv(ADMIT, sep="\t", dtype=str)
adf = adf[adf["status"] == "ADMIT"].copy()
adf["cv_r2"] = adf["cv_r2"].astype(float)
adf["chr"] = adf["sv_locus"].str.extract(r"^(chr[0-9XY]+)_")
adf["ids"] = adf["sv_id"].map(parse_ids)
inv_by_chr, meta = {}, {}
for ch, g in adf.groupby("chr"):
    inv = {}
    for _, r in g.iterrows():
        meta[r["sv_locus"]] = dict(cv=r["cv_r2"], strat=r["strategy"], svidx=r["sv_index"], chr=ch)
        for i in r["ids"]:
            inv.setdefault(i, []).append(r["sv_locus"])
    inv_by_chr[ch] = inv

def _metrics_byid(strat, ch):
    k=(strat,ch)
    if k not in svfreq_cache:
        p=f"{EN}/{strat}_mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200/{ch}/sv_metrics.tsv"
        try:
            mt=pd.read_csv(p,sep="\t",usecols=["svfreq","svpos","sv_id"])
            svfreq_cache[k]=dict(zip(mt["sv_id"].astype(str), zip(mt["svfreq"].astype(float), mt["svpos"])))
        except Exception: svfreq_cache[k]={}
    return svfreq_cache[k]
svfreq_cache={}
_LOCRAW = dict(zip(adf["sv_locus"], adf["sv_id"]))
for loc,m in meta.items():
    d=_metrics_byid(m["strat"], m["chr"]).get(str(_LOCRAW.get(loc)), (np.nan,None))
    m["svfreq"]=d[0]; m["svpos"]=(int(d[1]) if d[1] is not None and not pd.isna(d[1]) else None)

_LOCID = dict(zip(adf["sv_locus"], adf["ids"]))
def meta_ids(loc):
    return _LOCID.get(loc, set())

def base_var(ch):
    vf = f"{DBASE}/{ch}/{ch}_sv_variants.tsv"
    if not os.path.exists(vf): return None
    v = pd.read_csv(vf, sep="\t", dtype=str)
    v["ids"] = v["id"].map(parse_ids)
    v["key"] = v["pos"].astype(str)+"|"+v["svtype"]+"|"+v["svlen"].astype(str)
    return dict(zip(v["key"], v["ids"]))

chroms = sorted({re.match(r"(chr[0-9XY]+)_", os.path.basename(f)).group(1)
                 for f in glob.glob(f"{PERCH}/*_perm.tsv")}, key=lambda c:(len(c), c))
recs = []
for ch in chroms:
    k2ids = base_var(ch)
    if k2ids is None: continue
    inv = inv_by_chr.get(ch, {})
    for f in glob.glob(f"{PERCH}/{ch}_*_perm.tsv"):
        pheno = re.search(rf"{ch}_(.+)_perm.tsv", os.path.basename(f)).group(1)
        d = pd.read_csv(f, sep="\t")
        d = d[d["p_emp_perm"] <= HIT_P]
        for _, h in d.iterrows():
            key = f"{int(h['pos'])}|{h['svtype']}|{int(h['svlen'])}"
            members = k2ids.get(key, set())
            pos, svt, svl = int(h["pos"]), str(h["svtype"]), int(h["svlen"])
            af = float(h["af"])
            cands = []
            for loc in {l for i in members for l in inv.get(i, [])}:
                m = meta[loc]
                sh = len(members & meta_ids(loc))
                ov = sh/len(members) if members else 0
                svf_c = m.get("svfreq", np.nan)
                afd = abs(af - svf_c) if not (svf_c is None or pd.isna(svf_c)) else 9.9
                cands.append((-afd, m["cv"], ov, sh, svf_c, loc))
            if not cands:
                best=None; cls,cv,strat,svf,bovl,bshare="non_admit_drop",np.nan,"",np.nan,0.0,0
            else:
                negafd, cv, bovl, bshare, svf, best = max(cands, key=lambda c:(c[0], c[1]))
                strat = meta[best]["strat"]; afd = -negafd
                cls = "canonical_confirmed" if (afd <= 0.05 and bovl >= 0.5) else "admit_divergent"
            recs.append(dict(chrom=ch, pos=int(h["pos"]), svtype=h["svtype"], svlen=int(h["svlen"]),
                             phenotype=pheno, t_obs=round(float(h["t_obs"]),3),
                             p_emp=float(h["p_emp_perm"]), af=round(af,4),
                             admit_locus=best, admit_cv_r2=(round(cv,3) if best else np.nan),
                             strategy=strat, overlap_frac=round(bovl,3), n_shared=bshare,
                             admit_svfreq=(round(svf,4) if not np.isnan(svf) else np.nan),
                             fidelity=cls))

R = pd.DataFrame(recs).sort_values(["phenotype","p_emp","chrom"]).reset_index(drop=True)
R.to_csv(f"{OUTD}/canonical_hit_map.tsv", sep="\t", index=False)

print(f"[hits<=1e-4] n={len(R)} across {R['chrom'].nunique()} chroms, {R['phenotype'].nunique()} phenotypes")
print("\n=== fidelity class counts ===")
print(R["fidelity"].value_counts().to_string())
print("\n=== canonical_confirmed hits (the defensible canonical map) ===")
cc = R[R["fidelity"]=="canonical_confirmed"].sort_values("p_emp")
cols = ["chrom","pos","svtype","svlen","phenotype","t_obs","p_emp","af","admit_locus","admit_cv_r2","overlap_frac","admit_svfreq"]
print(cc[cols].head(25).to_string(index=False))

n_tests = 6519 * 8
bonf = 0.05 / n_tests
floor = 1.0/1_000_001
print(f"\n=== genome-wide multiple-testing on canonical 6519 x 8 phenotypes ===")
print(f"  n_tests={n_tests}  Bonferroni 0.05 threshold = {bonf:.3e}")
print(f"  permutation floor (1e6 perms) = {floor:.3e}  -> floor {'ABOVE' if floor>bonf else 'BELOW'} Bonferroni")
print(f"  => floor-saturated hits (p={floor:.1e}) are {'NOT resolvable' if floor>bonf else 'resolvable'} at genome-wide sig; need >=1e7 perms")
res = R[(R["fidelity"]=="canonical_confirmed") & (R["p_emp"]<=bonf)]
print(f"  canonical_confirmed hits clearing Bonferroni outright (p<={bonf:.2e}): {len(res)}")
fl  = R[(R["fidelity"]=="canonical_confirmed") & (R["p_emp"]<=floor*2)]
print(f"  canonical_confirmed at perm floor (need more perms to resolve): {len(fl)}")
print(f"    -> {sorted(set(fl['chrom']+':'+fl['svtype']+fl['svlen'].astype(str)))}")
print(f"\nwrote {OUTD}/canonical_hit_map.tsv")
