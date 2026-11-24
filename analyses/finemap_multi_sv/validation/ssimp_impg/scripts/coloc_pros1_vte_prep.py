#!/usr/bin/env python
import subprocess, os
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]
P = CENSV_WORKSPACE + ""
BASE = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
V = f"{BASE}/validation/ssimp_impg"
GD = f"{BASE}/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
BT = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
TOOL = f"{P}/_shared/tools/process_neale_lab_GWAS/v_master_83578d1_20260514"
LO = f"{TOOL}/.pixi/envs/default/bin/liftOver"; CHAIN = f"{TOOL}/hg19ToHg38.over.chain"
WD = f"{V}/work/coloc_pros1"; EQ = f"{WD}/eqtlgen_PROS1_hg19.tsv"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

rows = []
with open(EQ) as f:
    hdr = f.readline().rstrip("\n").split("\t"); ix = {h: i for i, h in enumerate(hdr)}
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        try:
            rows.append({"rs": c[ix["SNP"]], "chr": c[ix["SNPChr"]], "pos19": int(c[ix["SNPPos"]]),
                         "ae": c[ix["AssessedAllele"]].upper(), "oa": c[ix["OtherAllele"]].upper(),
                         "z": float(c[ix["Zscore"]]), "n": int(c[ix["NrSamples"]])})
        except Exception:
            pass
print(f"[prep] eQTLGen PROS1 cis SNPs (hg19) = {len(rows)}")

os.makedirs(WD, exist_ok=True)
bed, outb, unm = f"{WD}/eq_hg19.bed", f"{WD}/eq_hg38.bed", f"{WD}/eq_unmapped.bed"
with open(bed, "w") as o:
    for i, r in enumerate(rows):
        o.write(f"chr{r['chr']}\t{r['pos19']-1}\t{r['pos19']}\t{i}\n")
subprocess.run([LO, bed, CHAIN, outb, unm], check=True, capture_output=True, text=True)
b38 = {}
with open(outb) as f:
    for ln in f:
        c = ln.split("\t"); b38[int(c[3])] = int(c[2])
print(f"[prep] lifted {len(b38)}/{len(rows)} to b38 "
      f"(unmapped={len(rows)-len(b38)})")

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

poss = list(b38.values()); lo, hi = min(poss) - 1, max(poss) + 1
print(f"[prep] b38 cis window 3:{lo}-{hi}")
for study in ["finngen", "ukbb_neale"]:
    g = parse_gwas(study, lo, hi)
    outp = f"{WD}/coloc_input_{study}.tsv"
    nmatch = nflip = npal = 0
    with open(outp, "w") as o:
        o.write("snp\tpos\ta_eff\ta_oth\tz_eqtl\tn_eqtl\tbeta_gwas\tse_gwas\tmaf\tn_gwas\n")
        for i, r in enumerate(rows):
            if i not in b38:
                continue
            p38 = b38[i]
            if p38 not in g:
                continue
            ref, alt, maf, n, b, se = g[p38]; ae, oa = r["ae"], r["oa"]
            if {ae, oa} != {ref, alt}:
                continue
            if ae == COMP.get(oa, ""):
                npal += 1
            if ae == alt:
                z = r["z"]
            elif ae == ref:
                z = -r["z"]; nflip += 1
            else:
                continue
            nmatch += 1
            o.write(f"chr3_{p38}_{ref}_{alt}\t{p38}\t{alt}\t{ref}\t{z:.6f}\t{r['n']}\t{b:.6g}\t{se:.6g}\t{maf:.6g}\t{n}\n")
    print(f"[prep] {study}: common SNPs matched={nmatch} (sign-flipped={nflip}, palindromic={npal}) -> {outp}")
print("PREP_DONE")
