#!/usr/bin/env python
import subprocess
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]
P = CENSV_WORKSPACE + ""
BASE = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
V = f"{BASE}/validation/ssimp_impg"
GD = f"{BASE}/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
BT = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
WD = f"{V}/work/coloc_pros1"; EQ = f"{WD}/proteinS_pqtl_PROS1.b38.tsv"
rows = {}
with open(EQ) as f:
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        try:
            pos = int(c[1]); ea = c[2].upper(); oa = c[3].upper()
            beta = float(c[4]); se = float(c[5]); eaf = float(c[6])
        except Exception:
            continue
        if se > 0 and 0 < eaf < 1 and pos not in rows:
            rows[pos] = (ea, oa, beta, se, min(eaf, 1 - eaf))
print(f"[prep-pqtl] Protein-S pQTL cis SNPs (b38) = {len(rows)}")

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
print(f"[prep-pqtl] b38 cis window 3:{lo}-{hi}")
for study in ["finngen", "ukbb_neale"]:
    g = parse_gwas(study, lo, hi); nmatch = nflip = 0
    outp = f"{WD}/coloc_input_pqtl_{study}.tsv"
    with open(outp, "w") as o:
        o.write("snp\tpos\tbeta_pqtl\tse_pqtl\tmaf_pqtl\tbeta_gwas\tse_gwas\tmaf_gwas\tn_gwas\n")
        for pos, (ea, oa, beta, se, maf_p) in rows.items():
            if pos not in g:
                continue
            ref, alt, maf, n, b, sg = g[pos]
            if {ea, oa} != {ref, alt}:
                continue
            if ea == alt:
                bp = beta
            elif ea == ref:
                bp = -beta; nflip += 1
            else:
                continue
            nmatch += 1
            o.write(f"chr3_{pos}_{ref}_{alt}\t{pos}\t{bp:.6g}\t{se:.6g}\t{maf_p:.6g}\t{b:.6g}\t{sg:.6g}\t{maf:.6g}\t{n}\n")
    print(f"[prep-pqtl] {study}: common SNPs matched={nmatch} (flipped={nflip}) -> {outp}")
print("PREP_PQTL_DONE")
