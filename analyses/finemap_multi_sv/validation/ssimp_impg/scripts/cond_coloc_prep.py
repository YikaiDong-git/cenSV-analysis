#!/usr/bin/env python
import numpy as np, subprocess, os
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]
P = CENSV_WORKSPACE + ""
BASE = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
V = f"{BASE}/validation/ssimp_impg"
BT = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
BCF = f"{BT}/bcftools"
VCF = f"{P}/Project1_Centromere/resources/1000G/1kGP_high_coverage_Illumina.chr3.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PANEL = f"{P}/_shared/datasets/1KGP/v_canonical_phase3_20130502/EUR/integrated_call_samples_v3.20130502.ALL.panel"
WD = f"{V}/work/coloc_pros1"; TMP = "/tmp/cond_coloc"; os.makedirs(TMP, exist_ok=True)
EUR_POPS = {"FIN", "CEU", "GBR", "TSI", "IBS"}
STUDY = "finngen"
INP = f"{WD}/coloc_input_pqtl_{STUDY}.tsv"

snps, pos, ref, alt, bp, sp, mp, bg, sg, mg = [], [], [], [], [], [], [], [], [], []
ng = None
with open(INP) as f:
    next(f)
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        parts = c[0].split("_")
        snps.append(c[0]); pos.append(int(c[1])); ref.append(parts[-2].upper()); alt.append(parts[-1].upper())
        bp.append(float(c[2])); sp.append(float(c[3])); mp.append(float(c[4]))
        bg.append(float(c[5])); sg.append(float(c[6])); mg.append(float(c[7])); ng = int(float(c[8]))
n0 = len(snps)
print(f"[cond-prep] coloc-input common SNPs = {n0}")

vcf_s = set(subprocess.run([BCF, "query", "-l", VCF], capture_output=True, text=True, check=True).stdout.split())
eur = []
with open(PANEL) as fh:
    next(fh)
    for ln in fh:
        f = ln.split()
        if len(f) >= 3 and f[1] in EUR_POPS and f[0] in vcf_s:
            eur.append(f[0])
open(f"{TMP}/eur.txt", "w").write("\n".join(eur) + "\n")
print(f"[cond-prep] EUR samples = {len(eur)}")

with open(f"{TMP}/regions.txt", "w") as fh:
    for p in sorted(set(pos)):
        fh.write(f"chr3\t{p}\t{p}\n")

cmd = (f"{BCF} view -S {TMP}/eur.txt -R {TMP}/regions.txt -v snps -m2 -M2 {VCF} "
       f"| {BCF} query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'")
out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
posidx = {p: i for i, p in enumerate(pos)}
G = np.full((n0, len(eur)), np.nan)
got = set()
for ln in out.splitlines():
    f = ln.split("\t"); p = int(f[0]); rr = f[1].upper(); aa = f[2].upper()
    if p not in posidx:
        continue
    k = posidx[p]
    if {rr, aa} != {ref[k], alt[k]}:
        continue
    flip = (rr == alt[k] and aa == ref[k])
    gts = f[3:]; dos = np.empty(len(gts))
    for j, g in enumerate(gts):
        a0 = g[0]; a1 = g[2] if len(g) >= 3 else "."
        if a0 == "." or a1 == ".":
            dos[j] = np.nan
        else:
            cnt = (a0 == "1") + (a1 == "1")
            dos[j] = (2 - cnt) if flip else cnt
    m = np.isnan(dos)
    if m.any():
        dos[m] = np.nanmean(dos)
    G[k, :] = dos; got.add(k)
keep = sorted(got)
print(f"[cond-prep] SNPs with EUR genotypes = {len(keep)}/{n0}")

G = G[keep, :]
sd = np.nanstd(G, axis=1)
poly = np.where(sd > 0)[0]
keep2 = [keep[i] for i in poly]
G = G[poly, :]
R = np.corrcoef(G)
R = np.nan_to_num(R, nan=0.0); np.fill_diagonal(R, 1.0)
idx = list(keep2)
with open(f"{WD}/cond_coloc_{STUDY}_data.tsv", "w") as o:
    o.write("snp\tpos\tbeta_pqtl\tvarbeta_pqtl\tmaf\tbeta_gwas\tvarbeta_gwas\tmaf_gwas\tn_gwas\n")
    for i in idx:
        o.write(f"{snps[i]}\t{pos[i]}\t{bp[i]:.6g}\t{sp[i]**2:.6g}\t{mp[i]:.6g}\t"
                f"{bg[i]:.6g}\t{sg[i]**2:.6g}\t{mg[i]:.6g}\t{ng}\n")
np.savetxt(f"{WD}/cond_coloc_{STUDY}_LD.tsv", R, fmt="%.5f", delimiter="\t")
print(f"[cond-prep] final SNPs (polymorphic, with LD) = {len(idx)} ; "
      f"wrote cond_coloc_{STUDY}_data.tsv + _LD.tsv ({R.shape[0]}x{R.shape[1]})")
print("COND_PREP_DONE")
