#!/usr/bin/env python
import subprocess, os, csv, math
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P    = CENSV_WORKSPACE + ""
BASE = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
GD   = f"{BASE}/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
BT   = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
TABIX = f"{BT}/tabix"
HAND = f"{P}/Project1_Centromere/analyses/disease_replication/aou_handoff"
AGG  = f"{BASE}/validation/ssimp_impg/work/phewide_agg"
OUT  = f"{HAND}/discovery_sign_table.tsv"

W = {}
with open(f"{HAND}/sv_weights_long.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        W.setdefault(r["sv_locus"], []).append(
            (int(r["pos"]), r["ref"].upper(), r["alt"].upper(), float(r["beta_std"])))

META = {}
with open(f"{HAND}/sv_metadata.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        META[r["sv_locus"]] = r

AUD = {}
with open(f"{AGG}/candidate_audit_v2_owntag.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        AUD[r["sv_locus"]] = (r["study"], r["rep_trait"], float(r["z_SV"]))

_cache = {}
def zmap_for(study, trait, chrom, lo, hi):
    key = (study, trait, chrom, lo, hi)
    if key in _cache:
        return _cache[key]
    out = subprocess.run([TABIX, f"{GD}/{study}/{trait}.tsv.bgz", f"{chrom}:{lo}-{hi}"],
                         capture_output=True, text=True, check=True).stdout
    dz = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        try:
            pp = f[2].split("_"); rr, aa = pp[-2].upper(), pp[-1].upper()
            se = float(f[6]); bb = float(f[5])
            if se > 0:
                dz[(int(f[1]), rr, aa)] = bb / se
        except Exception:
            pass
    _cache[key] = dz
    return dz

def bz_for(sv, study, trait):
    wl = W.get(sv, [])
    if not wl:
        return 0.0, 0
    chrom = META[sv]["region"].split("_")[0].replace("chr", "")
    lo = min(p for (p, r, a, b) in wl) - 1000
    hi = max(p for (p, r, a, b) in wl) + 1000
    dz = zmap_for(study, trait, chrom, lo, hi)
    bz = 0.0; nobs = 0
    for (p, r, a, b) in wl:
        zz = dz.get((p, r, a))
        if zz is None:
            zz2 = dz.get((p, a, r))
            zz = (-zz2) if zz2 is not None else None
        if zz is not None and math.isfinite(zz):
            bz += b * zz; nobs += 1
    return bz, nobs

def sgn(x):
    return "+" if x > 0 else ("-" if x < 0 else "0")

rows = []
for sv, meta in META.items():
    bz, nobs = bz_for(sv, meta["rep_study"], meta["rep_trait"])
    aud = AUD.get(sv)
    rows.append(dict(sv_locus=sv, rep_study=meta["rep_study"], rep_trait=meta["rep_trait"],
                     n_weight_snps=len(W.get(sv, [])), n_obs_in_gwas=nobs,
                     discovery_bz=f"{bz:.4f}", discovery_sign=sgn(bz),
                     in_robust24=meta["in_robust24"],
                     audit_trait=(aud[1] if aud else ""),
                     audit_zSV=(f"{aud[2]:.4f}" if aud else "")))

nchk = nmis = 0; mism = []
for sv, (astudy, atrait, azsv) in AUD.items():
    abz, anobs = bz_for(sv, astudy, atrait)
    nchk += 1
    if sgn(abz) != sgn(azsv):
        nmis += 1; mism.append((sv, atrait, azsv, round(abz, 4), anobs))
print(f"VALIDATION (b'z on audit trait) vs audit z_SV sign: checked={nchk} mismatches={nmis}")
for m in mism:
    print("   MISMATCH", m)

diff = [(r["sv_locus"], r["audit_trait"], r["rep_trait"])
        for r in rows if r["audit_trait"] and r["audit_trait"] != r["rep_trait"]]
if diff:
    print("NOTE handoff rep_trait != audit trait (deliverable sign is on the HANDOFF trait):")
    for d in diff:
        print("   ", d)

zero = [r["sv_locus"] for r in rows if r["n_obs_in_gwas"] == 0]
if zero:
    print(f"WARN {len(zero)} SVs had 0 weight SNPs in the discovery sumstat (sign=0): {zero}")

cols = ["sv_locus","rep_study","rep_trait","n_weight_snps","n_obs_in_gwas",
        "discovery_bz","discovery_sign","in_robust24","audit_trait","audit_zSV"]
tmp = OUT + ".tmp"
with open(tmp, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: (r["in_robust24"] != "True", r["sv_locus"])):
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
os.replace(tmp, OUT)
print(f"WROTE {OUT} ({len(rows)} SVs)")
print("DONE_DISCOVERY_SIGN")
