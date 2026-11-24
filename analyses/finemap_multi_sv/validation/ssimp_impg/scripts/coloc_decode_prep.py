#!/usr/bin/env python
import numpy as np, subprocess, os
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]
TAB = chr(9); NL = chr(10)
P = CENSV_WORKSPACE + ""
BASE = P + "/Project1_Centromere/analyses/finemap_multi_sv"
V = BASE + "/validation/ssimp_impg"
GD = BASE + "/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
BT = P + "/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
BCF = BT + "/bcftools"
VCF = P + "/Project1_Centromere/resources/1000G/1kGP_high_coverage_Illumina.chr3.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PANEL = P + "/_shared/datasets/1KGP/v_canonical_phase3_20130502/EUR/integrated_call_samples_v3.20130502.ALL.panel"
WD = V + "/work/coloc_pros1"
SLICE = WD + "/decode_proteinS_PROS1.b38.slice.tsv"
TMP = "/tmp/cond_coloc_decode"; os.makedirs(TMP, exist_ok=True)
EUR_POPS = set(["FIN", "CEU", "GBR", "TSI", "IBS"])

rows = {}; nmax = 0
with open(SLICE) as fh:
    for ln in fh:
        c = ln.rstrip(NL).split(TAB)
        if len(c) < 12:
            continue
        try:
            pos = int(c[1]); ea = c[4].upper(); oa = c[5].upper()
            beta = float(c[6]); se = float(c[9]); n = int(float(c[10])); maf = float(c[11])
        except Exception:
            continue
        if se > 0 and 0.0 < maf < 1.0 and pos not in rows:
            rows[pos] = (ea, oa, beta, se, min(maf, 1.0 - maf))
            if n > nmax: nmax = n
print("[decode-prep] deCODE Protein-S cis SNPs (b38) = %d ; Nmax = %d" % (len(rows), nmax))
open(WD + "/decode_npqtl.txt", "w").write(str(nmax) + NL)

def parse_gwas(study, lo, hi):
    out = subprocess.run([BT + "/tabix", GD + "/" + study + "/I9_VTE.tsv.bgz", "3:%d-%d" % (lo, hi)],
                         capture_output=True, text=True, check=True)
    d = {}
    for ln in out.stdout.splitlines():
        c = ln.split(TAB)
        try:
            pos = int(c[1]); parts = c[2].split("_"); ref = parts[-2].upper(); alt = parts[-1].upper()
            maf = float(c[3]); n = int(float(c[4])); b = float(c[5]); se = float(c[6])
        except Exception:
            continue
        if se > 0 and 0.0 < maf < 1.0:
            d[pos] = (ref, alt, maf, n, b, se)
    return d

poss = sorted(rows.keys()); lo = poss[0] - 1; hi = poss[-1] + 1
print("[decode-prep] b38 cis window 3:%d-%d" % (lo, hi))
finn_inp = None
for study in ["finngen", "ukbb_neale"]:
    g = parse_gwas(study, lo, hi); nmatch = 0; nflip = 0
    outp = WD + "/coloc_input_decode_" + study + ".tsv"
    with open(outp, "w") as o:
        o.write(TAB.join(["snp", "pos", "beta_pqtl", "se_pqtl", "maf_pqtl", "beta_gwas", "se_gwas", "maf_gwas", "n_gwas"]) + NL)
        for pos in poss:
            if pos not in g: continue
            ea, oa, beta, se, maf_p = rows[pos]
            ref, alt, maf, n, b, sg = g[pos]
            if set([ea, oa]) != set([ref, alt]): continue
            if ea == alt: bp = beta
            elif ea == ref: bp = -beta; nflip += 1
            else: continue
            nmatch += 1
            o.write(TAB.join(["chr3_%d_%s_%s" % (pos, ref, alt), str(pos), "%.6g" % bp, "%.6g" % se,
                              "%.6g" % maf_p, "%.6g" % b, "%.6g" % sg, "%.6g" % maf, str(n)]) + NL)
    print("[decode-prep] %s: matched=%d (flipped=%d) -> %s" % (study, nmatch, nflip, outp))
    if study == "finngen": finn_inp = outp

snps = []; pos = []; ref = []; alt = []; bp = []; sp = []; mp = []; bg = []; sg = []; mg = []; ng = None
with open(finn_inp) as fh:
    next(fh)
    for ln in fh:
        c = ln.rstrip(NL).split(TAB); parts = c[0].split("_")
        snps.append(c[0]); pos.append(int(c[1])); ref.append(parts[-2].upper()); alt.append(parts[-1].upper())
        bp.append(float(c[2])); sp.append(float(c[3])); mp.append(float(c[4]))
        bg.append(float(c[5])); sg.append(float(c[6])); mg.append(float(c[7])); ng = int(float(c[8]))
n0 = len(snps); print("[decode-prep] FinnGen coloc-input common SNPs = %d" % n0)
if n0 < 5:
    print("DECODE_PREP_DONE"); raise SystemExit(0)
vcf_s = set(subprocess.run([BCF, "query", "-l", VCF], capture_output=True, text=True, check=True).stdout.split())
eur = []
with open(PANEL) as fh:
    next(fh)
    for ln in fh:
        f = ln.split()
        if len(f) >= 3 and f[1] in EUR_POPS and f[0] in vcf_s: eur.append(f[0])
open(TMP + "/eur.txt", "w").write(NL.join(eur) + NL)
with open(TMP + "/regions.txt", "w") as fh:
    for p in sorted(set(pos)): fh.write("chr3" + TAB + str(p) + TAB + str(p) + NL)
fmt = "%POS" + TAB + "%REF" + TAB + "%ALT[" + TAB + "%GT]" + NL
cmd = BCF + " view -S " + TMP + "/eur.txt -R " + TMP + "/regions.txt -v snps -m2 -M2 " + VCF + " | " + BCF + " query -f '" + fmt + "'"
out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
posidx = {}
for i, p in enumerate(pos): posidx[p] = i
G = np.full((n0, len(eur)), np.nan); got = set()
for ln in out.splitlines():
    f = ln.split(TAB)
    if len(f) < 4: continue
    p = int(f[0]); rr = f[1].upper(); aa = f[2].upper()
    if p not in posidx: continue
    k = posidx[p]
    if set([rr, aa]) != set([ref[k], alt[k]]): continue
    flip = (rr == alt[k] and aa == ref[k]); gts = f[3:]; dos = np.empty(len(gts))
    for j, gg in enumerate(gts):
        a0 = gg[0:1]; a1 = gg[2:3] if len(gg) >= 3 else "."
        if a0 == "." or a1 == ".": dos[j] = np.nan
        else:
            cnt = (1 if a0 == "1" else 0) + (1 if a1 == "1" else 0)
            dos[j] = (2 - cnt) if flip else cnt
    m = np.isnan(dos)
    if m.any(): dos[m] = np.nanmean(dos)
    G[k, :] = dos; got.add(k)
keep = sorted(got); print("[decode-prep] FinnGen SNPs with EUR genotypes = %d/%d" % (len(keep), n0))
G = G[keep, :]; sd = np.nanstd(G, axis=1); poly = np.where(sd > 0)[0]; keep2 = [keep[i] for i in poly]; G = G[poly, :]
R = np.corrcoef(G); R = np.nan_to_num(R, nan=0.0); np.fill_diagonal(R, 1.0); idx = list(keep2)
with open(WD + "/cond_coloc_decode_finngen_data.tsv", "w") as o:
    o.write(TAB.join(["snp", "pos", "beta_pqtl", "varbeta_pqtl", "maf", "beta_gwas", "varbeta_gwas", "maf_gwas", "n_gwas"]) + NL)
    for i in idx:
        o.write(TAB.join([snps[i], str(pos[i]), "%.6g" % bp[i], "%.6g" % (sp[i] * sp[i]), "%.6g" % mp[i],
                          "%.6g" % bg[i], "%.6g" % (sg[i] * sg[i]), "%.6g" % mg[i], str(ng)]) + NL)
np.savetxt(WD + "/cond_coloc_decode_finngen_LD.tsv", R, fmt="%.5f", delimiter=TAB)
print("[decode-prep] FinnGen final SNPs = %d ; wrote cond_coloc_decode_finngen_data.tsv + _LD.tsv" % len(idx))
print("DECODE_PREP_DONE")
