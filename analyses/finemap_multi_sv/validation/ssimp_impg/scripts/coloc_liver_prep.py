#!/usr/bin/env python
import subprocess
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]
P = CENSV_WORKSPACE + ""
BASE = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
V = f"{BASE}/validation/ssimp_impg"
GD = f"{BASE}/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
BT = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
WD = f"{V}/work/coloc_pros1"; EQ = f"{WD}/gtex_liver_PROS1.b38.tsv"

rows = {}
with open(EQ) as f:
    next(f)
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        try:
            pos = int(c[2]); ref = c[3].upper(); alt = c[4].upper()
            maf = float(c[7]); beta = float(c[9]); se = float(c[10])
        except Exception:
            continue
        if se > 0 and 0 < maf < 1 and pos not in rows:
            rows[pos] = (ref, alt, maf, beta, se)
print(f"[prep-liver] GTEx Liver PROS1 cis SNPs = {len(rows)}")

def parse_gwas(study, lo, hi):
    out = subprocess.run([f"{BT}/tabix", f"{GD}/{study}/I9_VTE.tsv.bgz", f"3:{lo}-{hi}"],
                         capture_output=True, text=True, check=True)
    d = {}
    for ln in out.stdout.splitlines():
        c = ln.split("\t")
        try:
            pos = int(c[1]); parts = c[2].split("_"); ref, alt = parts[-2].upper(), parts[-1].upper()
            maf = float(c[3]); n = int(float(c[4])); b = float(c[5]); se = float(c[6])
        except Exception:
            continue
        if se > 0 and 0 < maf < 1:
            d[pos] = (ref, alt, maf, n, b, se)
    return d

poss = list(rows.keys()); lo, hi = min(poss) - 1, max(poss) + 1
print(f"[prep-liver] b38 cis window 3:{lo}-{hi}")
for study in ["finngen", "ukbb_neale"]:
    g = parse_gwas(study, lo, hi); nmatch = nflip = 0
    outp = f"{WD}/coloc_input_gtexliver_{study}.tsv"
    with open(outp, "w") as o:
        o.write("snp\tpos\tbeta_eqtl\tse_eqtl\tmaf_eqtl\tbeta_gwas\tse_gwas\tmaf_gwas\tn_gwas\n")
        for pos, (eref, ealt, emaf, eb, ese) in rows.items():
            if pos not in g:
                continue
            ref, alt, maf, n, b, se = g[pos]
            if {eref, ealt} != {ref, alt}:
                continue
            if ealt == alt:
                be = eb
            elif ealt == ref:
                be = -eb; nflip += 1
            else:
                continue
            nmatch += 1
            o.write(f"chr3_{pos}_{ref}_{alt}\t{pos}\t{be:.6g}\t{ese:.6g}\t{emaf:.6g}\t{b:.6g}\t{se:.6g}\t{maf:.6g}\t{n}\n")
    print(f"[prep-liver] {study}: common SNPs matched={nmatch} (sign-flipped={nflip}) -> {outp}")
print("PREP_LIVER_DONE")
