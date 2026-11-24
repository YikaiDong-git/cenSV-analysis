#!/usr/bin/env python
"""build_phewide_all_handoff.py -- FULL-phewide AoU weight handoff."""
import numpy as np, subprocess, os, csv
from collections import defaultdict
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P = CENSV_WORKSPACE + ""
SSIMP = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
AGG = f"{SSIMP}/validation/ssimp_impg/work/phewide_agg"
BUND = f"{SSIMP}/output/bundles"
GD = f"{SSIMP}/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
TABIX = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin/tabix"
CT = f"{P}/Project1_Centromere/analyses/disease_replication"
OUT = f"{CT}/phewide_all_handoff"
os.makedirs(OUT, exist_ok=True)
FIT_PANEL = "all_1kGP_n188"
MAXROWS = 45000

rows = list(csv.DictReader(open(f"{AGG}/phewide_pip_ge_0.01.tsv"), delimiter="\t"))
sv_pipmeta = {}
maxpip = defaultdict(float)
npairs05 = defaultdict(int)
pairs05 = []
for r in rows:
    sv = r["sv_locus"]
    if sv not in sv_pipmeta:
        sv_pipmeta[sv] = dict(region=r["region_id"], category=r["sv_category"],
                              strategy=r["sv_strategy"], arm=r["sv_arm"])
    mx = max(float(r["impute_pip"]), float(r["matched_pip"]))
    if mx > maxpip[sv]:
        maxpip[sv] = mx
    if float(r["impute_pip"]) >= 0.5 or float(r["matched_pip"]) >= 0.5:
        npairs05[sv] += 1
        pairs05.append(r)
svset = sorted(sv_pipmeta)
region_svs = defaultdict(list)
for sv in svset:
    region_svs[sv_pipmeta[sv]["region"]].append(sv)
print(f"SV set: {len(svset)} sv_locus in {len(region_svs)} regions; pip>=0.5 pairs: {len(pairs05)}", flush=True)

W = {}
wrows_by_chrom = defaultdict(list)
maf_by_chrom = defaultdict(list)
svmeta_out = {}
for region, svs in region_svs.items():
    chrom = region.split("_")[0].replace("chr", "")
    d = np.load(f"{BUND}/chr{chrom}/{region}_bundle.npz", allow_pickle=True)
    svn = np.array([str(x) for x in d["sv_locus"]])
    pos = d["snp_pos"].astype(int)
    ref = np.array([str(x).upper() for x in d["snp_ref"]])
    alt = np.array([str(x).upper() for x in d["snp_alt"]])
    maf = d["snp_maf_panel"].astype(float)
    schr = d["snp_chrom"]
    cvr2arr = d["sv_cv_r2"]
    for sv in svs:
        idx = np.where(svn == sv)[0]
        if len(idx) == 0:
            print(f"WARN sv_locus not in bundle: {sv}", flush=True)
            continue
        ti = int(idx[0])
        w = d["B"][:, ti].astype(float)
        sup = np.where(w != 0)[0]
        wl = []
        for j in sup:
            cc = str(schr[j])
            cc = cc if cc.startswith("chr") else f"chr{cc}"
            wrows_by_chrom[cc].append((sv, cc, int(pos[j]), ref[j], alt[j], alt[j], f"{float(w[j]):.10g}"))
            maf_by_chrom[cc].append((sv, cc, int(pos[j]), ref[j], alt[j], f"{float(maf[j]):.6f}"))
            wl.append((int(pos[j]), ref[j], alt[j], float(w[j])))
        W[sv] = (chrom, wl)
        svmeta_out[sv] = dict(region=region, chrom=chrom, cv_r2=float(cvr2arr[ti]), n_weight_snps=len(sup))

zcache = {}
missing_bgz = set()


def zmap(study, trait, chrom, lo, hi):
    key = (study, trait, chrom, lo, hi)
    if key in zcache:
        return zcache[key]
    bgz = f"{GD}/{study}/{trait}.tsv.bgz"
    if not os.path.exists(bgz):
        missing_bgz.add((study, trait))
        zcache[key] = {}
        return {}
    try:
        out = subprocess.run([TABIX, bgz, f"{chrom}:{lo}-{hi}"], capture_output=True, text=True, check=True).stdout
    except Exception:
        zcache[key] = {}
        return {}
    dz = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        try:
            pp = f[2].split("_")
            rr, aa = pp[-2].upper(), pp[-1].upper()
            se, bb = float(f[6]), float(f[5])
            if se > 0:
                dz[(int(f[1]), rr, aa)] = bb / se
        except Exception:
            pass
    zcache[key] = dz
    return dz


def bz(sv, study, trait):
    chrom, wl = W.get(sv, (None, []))
    if not wl:
        return None, 0
    lo = min(p for (p, r, a, b) in wl) - 1000
    hi = max(p for (p, r, a, b) in wl) + 1000
    dz = zmap(study, trait, chrom, lo, hi)
    s, n = 0.0, 0
    for (p, r, a, b) in wl:
        zz = dz.get((p, r, a))
        if zz is None:
            z2 = dz.get((p, a, r))
            zz = (-z2) if z2 is not None else None
        if zz is not None and zz == zz:
            s += b * zz
            n += 1
    return (s if n > 0 else None), n


def sgn(x):
    return "+" if (x is not None and x > 0) else ("-" if (x is not None and x < 0) else "NA")


prows = []
for r in pairs05:
    sv, study, trait = r["sv_locus"], r["study"], r["pheno"]
    b, n = bz(sv, study, trait)
    prows.append([sv, study, trait, f'{svmeta_out[sv]["cv_r2"]:.3f}', r["impute_pip"], r["matched_pip"],
                  (f"{b:.4f}" if b is not None else "NA"), sgn(b), n])

def write_split(prefix, header, by_chrom):
    files = []
    for cc in sorted(by_chrom, key=lambda c: (len(c), c)):
        rr = sorted(by_chrom[cc], key=lambda x: (x[0], x[2]))
        nparts = (len(rr) + MAXROWS - 1) // MAXROWS
        for k in range(nparts):
            chunk = rr[k * MAXROWS:(k + 1) * MAXROWS]
            suff = "" if nparts == 1 else f"_part{k + 1}"
            fn = f"{OUT}/{prefix}_{cc}{suff}.tsv"
            with open(fn, "w") as fh:
                fh.write("\t".join(header) + "\n")
                for row in chunk:
                    fh.write("\t".join(map(str, row)) + "\n")
            files.append((os.path.basename(fn), len(chunk), os.path.getsize(fn)))
    return files


wf = write_split("sv_weights_long", ["sv_locus", "chrom", "pos", "ref", "alt", "a1_effect", "beta_std"], wrows_by_chrom)
mf = write_split("train_maf", ["sv_locus", "chrom", "pos", "ref", "alt", "train_maf_panel"], maf_by_chrom)

with open(f"{OUT}/sv_metadata.tsv", "w") as fh:
    fh.write("sv_locus\tregion\tchrom\tcv_r2\tn_weight_snps\tsv_category\tsv_strategy\tsv_arm\tn_pairs_pip05\tmax_pip\tfit_panel\n")
    for sv in svset:
        m = svmeta_out.get(sv)
        if not m:
            continue
        p = sv_pipmeta[sv]
        fh.write(f"{sv}\t{m['region']}\t{m['chrom']}\t{m['cv_r2']:.3f}\t{m['n_weight_snps']}\t"
                 f"{p['category']}\t{p['strategy']}\t{p['arm']}\t{npairs05[sv]}\t{maxpip[sv]:.3f}\t{FIT_PANEL}\n")

prows.sort(key=lambda r: (-(float(r[4]) if r[4] not in ("", "NA") else 0.0), r[0], r[2]))
with open(f"{OUT}/sv_trait_pairs.tsv", "w") as fh:
    fh.write("sv_locus\tstudy\ttrait\tcv_r2\timpute_pip\tmatched_pip\tdiscovery_bz\tdiscovery_sign\tn_obs_in_gwas\n")
    for row in prows:
        fh.write("\t".join(map(str, row)) + "\n")

maxsz = max([s for _, _, s in wf + mf], default=0)
print(f"weight files: {len(wf)}  maf files: {len(mf)}  max_file_bytes: {maxsz}  under_4.9MB: {maxsz < 4.9e6}", flush=True)
print(f"sv_metadata rows: {len(svmeta_out)}  sv_trait_pairs: {len(prows)}  bz_computed: {sum(1 for r in prows if r[6] != 'NA')}", flush=True)
if missing_bgz:
    print(f"missing GWAS bgz for {len(missing_bgz)} (study,trait) -> discovery_sign=NA there; e.g. {sorted(missing_bgz)[:5]}", flush=True)
print("DONE_PHEWIDE_ALL_BUILD", flush=True)
