#!/usr/bin/env python
import numpy as np, subprocess, os, csv
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P     = CENSV_WORKSPACE + ""
SSIMP = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
AGG   = f"{SSIMP}/validation/ssimp_impg/work/phewide_agg"
BUND  = f"{SSIMP}/output/bundles"
GD    = f"{SSIMP}/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
BT    = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
TABIX = f"{BT}/tabix"
CT    = f"{P}/Project1_Centromere/analyses/disease_replication"
OUT   = f"{CT}/phewide_222_handoff"; os.makedirs(OUT, exist_ok=True)

pairs = [r for r in csv.DictReader(open(f"{AGG}/phewide_pip_ge_0.01.tsv"), delimiter="\t")
         if float(r["impute_pip"]) >= 0.5 and float(r["matched_pip"]) >= 0.5]
assert len(pairs) == 222, f"expected 222 got {len(pairs)}"
svset = sorted({r["sv_locus"] for r in pairs}); assert len(svset) == 33
region_by_sv = {r["sv_locus"]: r["region_id"] for r in pairs}

robust = {r["sv_locus"] for r in csv.DictReader(open(f"{AGG}/candidate_audit_v2_owntag.tsv"), delimiter="\t")}
h50 = {r["sv_locus"] for r in csv.DictReader(open(f"{CT}/aou_handoff/sv_metadata.tsv"), delimiter="\t")}
existing_sign = {}
for r in csv.DictReader(open(f"{CT}/aou_handoff/discovery_sign_table.tsv"), delimiter="\t"):
    existing_sign[(r["sv_locus"], r["rep_trait"])] = r["discovery_sign"]

cache = {}
def bundle(region):
    chrom = region.split("_")[0].replace("chr", "")
    bp = f"{BUND}/chr{chrom}/{region}_bundle.npz"
    if bp not in cache: cache[bp] = np.load(bp, allow_pickle=True)
    return cache[bp]

wrows = []; mafrows = []; svmeta = {}; W = {}
for sv in svset:
    region = region_by_sv[sv]; chrom = region.split("_")[0].replace("chr", "")
    d = bundle(region)
    svn = np.array([str(x) for x in d["sv_locus"]]); ti = int(np.where(svn == sv)[0][0])
    w = d["B"][:, ti].astype(float); sup = np.where(w != 0)[0]
    pos = d["snp_pos"].astype(int)
    ref = np.array([str(x).upper() for x in d["snp_ref"]]); alt = np.array([str(x).upper() for x in d["snp_alt"]])
    maf = d["snp_maf_panel"].astype(float); wl = []
    for j in sup:
        cc = str(d["snp_chrom"][j]); cc = cc if cc.startswith("chr") else f"chr{cc}"
        wrows.append((sv, cc, int(pos[j]), ref[j], alt[j], alt[j], float(w[j])))
        mafrows.append((sv, cc, int(pos[j]), ref[j], alt[j], f"{float(maf[j]):.6f}"))
        wl.append((int(pos[j]), ref[j], alt[j], float(w[j])))
    W[sv] = (chrom, wl)
    svmeta[sv] = dict(region=region, chrom=chrom, cv_r2=float(d["sv_cv_r2"][ti]), n_weight_snps=len(sup))

zcache = {}
def zmap(study, trait, chrom, lo, hi):
    key = (study, trait, chrom, lo, hi)
    if key in zcache: return zcache[key]
    out = subprocess.run([TABIX, f"{GD}/{study}/{trait}.tsv.bgz", f"{chrom}:{lo}-{hi}"],
                         capture_output=True, text=True, check=True).stdout
    dz = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        try:
            pp = f[2].split("_"); rr, aa = pp[-2].upper(), pp[-1].upper()
            se = float(f[6]); bb = float(f[5])
            if se > 0: dz[(int(f[1]), rr, aa)] = bb / se
        except Exception: pass
    zcache[key] = dz; return dz

def bz(sv, study, trait):
    chrom, wl = W[sv]
    if not wl: return 0.0, 0
    lo = min(p for (p, r, a, b) in wl) - 1000; hi = max(p for (p, r, a, b) in wl) + 1000
    dz = zmap(study, trait, chrom, lo, hi); s = 0.0; n = 0
    for (p, r, a, b) in wl:
        zz = dz.get((p, r, a))
        if zz is None:
            z2 = dz.get((p, a, r)); zz = (-z2) if z2 is not None else None
        if zz is not None and zz == zz:
            s += b * zz; n += 1
    return s, n
def sgn(x): return "+" if x > 0 else ("-" if x < 0 else "0")

ntr = {}
for r in pairs: ntr[r["sv_locus"]] = ntr.get(r["sv_locus"], 0) + 1

prows = []; nchk = nmis = 0
for r in pairs:
    sv = r["sv_locus"]; study = r["study"]; trait = r["pheno"]
    b, n = bz(sv, study, trait); cvr2 = svmeta[sv]["cv_r2"]
    tier = ("MEGA_PLEIOTROPIC_REGIONAL" if ntr[sv] >= 20
            else "SNP_EQUIVALENT" if cvr2 >= 0.7 else "SV_CANDIDATE")
    degen = int(float(r["region_n_sv_ge_robust_imp"]) > 3 or float(r["region_n_sv_ge_robust_mat"]) > 3)
    es = existing_sign.get((sv, trait))
    if es is not None:
        nchk += 1; nmis += (es != sgn(b))
    prows.append([sv, study, trait, f"{cvr2:.3f}", r["impute_pip"], r["matched_pip"],
                  r["converged_imp"], r["converged_mat"], r["region_n_sv_ge_robust_imp"], r["region_n_sv_ge_robust_mat"],
                  r["sv_category"], r["sv_strategy"], r["sv_arm"], f"{b:.4f}", sgn(b), n,
                  tier, degen, str(sv in robust), str(sv in h50), ntr[sv]])
print(f"VALIDATION vs existing discovery_sign_table (overlap): checked={nchk} mismatches={nmis}")

def wtab(path, header, rows):
    with open(path + ".tmp", "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows: fh.write("\t".join(map(str, r)) + "\n")
    os.replace(path + ".tmp", path)

TRANK = {"SV_CANDIDATE": 0, "SNP_EQUIVALENT": 1, "MEGA_PLEIOTROPIC_REGIONAL": 2}
prows.sort(key=lambda r: (TRANK[r[16]], -float(r[3]), r[0], r[2]))
wtab(f"{OUT}/pairs_222.tsv",
     ["sv_locus","study","trait","cv_r2","impute_pip","matched_pip","converged_imp","converged_mat",
      "region_n_sv_ge_robust_imp","region_n_sv_ge_robust_mat","sv_category","sv_strategy","sv_arm",
      "discovery_bz","discovery_sign","n_obs_in_gwas","tier","region_degenerate","in_robust24","in_handoff50","n_traits_hit"], prows)
wrows.sort(key=lambda r: (r[0], r[1], r[2]))
wtab(f"{OUT}/sv_weights_long_222.tsv", ["sv_locus","chrom","pos","ref","alt","a1_effect","beta_std"], wrows)
mafrows.sort(key=lambda r: (r[0], r[1], r[2]))
wtab(f"{OUT}/train_maf_222.tsv", ["sv_locus","chrom","pos","ref","alt","train_maf_panel"], mafrows)
mrows = []
for sv in svset:
    m = svmeta[sv]; cvr2 = m["cv_r2"]
    mtier = ("MEGA_PLEIOTROPIC_REGIONAL" if ntr[sv] >= 20
             else "SNP_EQUIVALENT" if cvr2 >= 0.7 else "SV_CANDIDATE")
    mrows.append([sv, m["region"], m["chrom"], f"{cvr2:.3f}", m["n_weight_snps"], ntr[sv],
                  mtier, str(sv in robust), str(sv in h50)])
mrows.sort(key=lambda r: (TRANK[r[6]], -float(r[3])))
wtab(f"{OUT}/sv_metadata_222.tsv", ["sv_locus","region","chrom","cv_r2","n_weight_snps","n_traits_hit","tier","in_robust24","in_handoff50"], mrows)
tr = {}
for r in pairs: tr[(r["study"], r["pheno"])] = tr.get((r["study"], r["pheno"]), 0) + 1
wtab(f"{OUT}/trait_list_222.tsv", ["study","trait","n_SVs"],
     sorted([[s, t, n] for (s, t), n in tr.items()], key=lambda x: (-x[2], x[0], x[1])))

from collections import Counter
tc = Counter(r[16] for r in prows); svtc = Counter(m[6] for m in mrows)
print(f"SVs={len(svset)} pairs={len(prows)} weight_rows={len(wrows)} unique_traits={len(tr)}")
print("tier pairs:", dict(tc))
print("tier SVs :", dict(svtc))
print(f"region_degenerate pairs={sum(1 for r in prows if r[17]==1)}  in_handoff50 pairs={sum(1 for r in prows if r[19]=='True')}")
print("DONE_222_BUILD")
