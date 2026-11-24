#!/usr/bin/env python
"""ld_locus_specificity.py — Does each imputed-SV predictor RE-TAG a known GWAS locus, or tag a distinct signal?"""
import numpy as np, subprocess, os
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P     = CENSV_WORKSPACE + ""
SSIMP = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
BUND  = f"{SSIMP}/output/bundles"
GD    = f"{P}/Project1_Centromere/analyses/sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
BT    = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
BCF   = f"{BT}/bcftools"; TABIX = f"{BT}/tabix"
RES   = f"{P}/Project1_Centromere/resources/1000G"
PANEL = f"{P}/_shared/datasets/1KGP/v_canonical_phase3_20130502/EUR/integrated_call_samples_v3.20130502.ALL.panel"
EUR_POPS = {"FIN", "CEU", "GBR", "TSI", "IBS"}
GWS    = 5.4513
MARGIN = 1_000_000
TMP    = "/tmp/ld_locus_spec"; os.makedirs(TMP, exist_ok=True)

SVS = [
  ("chr3_sv20_DEL_non_best6_alt_edge1000",  "chr3_super1_85533081_94536472",   "3",  "neale",   "21002_irnt", "-", "body_weight(neale) NEW REPLICATED (AoU b=-0.087 q=9e-4); prior LD check: NOT CADM2"),
  ("chr4_sv349_DEL_non_best8_fksize200k",   "chr4_super1_47409879_53011266",   "4",  "neale",   "50_irnt",    "-", "height(neale)   REPLICATED  (AoU b=-0.047 q=0.015)"),
  ("chr11_sv263_DEL_non_best8_edge500",     "chr11_super1_49844498_54927998",  "11", "finngen", "BMI_IRN",    "-", "BMI(finngen)    REPLICATED  (AoU b=-0.060 q=0.020)"),
  ("chr4_sv398_INS_non_best8_fksize100k",   "chr4_super1_47409879_53011266",   "4",  "neale",   "30110_irnt", "-", "PDW/platelet(neale) NOMINAL survivor (AoU cond-robust) - novelty UNTESTED"),
  ("chr12_sv721_DEL_non_best6_alt_edge1000","chr12_super1_32012658_38833367",  "12", "finngen", "HEIGHT_IRN", "+", "height(finngen) NON-REPLIC (AoU EUR NS p=0.12; all p=0.02 flip; q=0.35)"),
]

def eur_samples(chrom):
    vcf = f"{RES}/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    vs = set(subprocess.run([BCF, "query", "-l", vcf], capture_output=True, text=True, check=True).stdout.split())
    eur = []
    with open(PANEL) as fh:
        next(fh)
        for ln in fh:
            f = ln.split()
            if len(f) >= 3 and f[0] in vs and f[1] in EUR_POPS: eur.append(f[0])
    return vcf, eur

def read_geno_pos(vcf, eurfile, positions, chrom, tag):
    tf = f"{TMP}/pos_{chrom}_{tag}.txt"
    with open(tf, "w") as fh:
        for p in sorted(set(int(x) for x in positions)): fh.write(f"chr{chrom}\t{p}\t{p}\n")
    cmd = (f"{BCF} view -S {eurfile} -R {tf} -v snps -m2 -M2 {vcf} "
           f"| {BCF} query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    g = {}
    for ln in out.splitlines():
        f = ln.split("\t"); p = int(f[0]); r = f[1].upper(); a = f[2].upper(); gts = f[3:]
        dos = np.array([(gg[0] == "1") + (gg[2] == "1") if (len(gg) >= 3 and gg[0] != "." and gg[2] != ".") else np.nan for gg in gts], float)
        m = np.isnan(dos)
        if m.all(): continue
        dos[m] = np.nanmean(dos); g[p] = (r, a, dos)
    return g

for svid, region, chrom, study, trait, dsign, label in SVS:
    print(f"\n{'='*90}\n{svid}\n  {label} | region {region} | discovery_sign {dsign}\n{'='*90}")
    d = np.load(f"{BUND}/chr{chrom}/{region}_bundle.npz", allow_pickle=True)
    svn = np.array([str(x) for x in d["sv_locus"]]); ti = int(np.where(svn == svid)[0][0])
    w = d["B"][:, ti].astype(float); sup = np.where(w != 0)[0]
    spos = d["snp_pos"].astype(int)[sup]
    sref = np.array([str(x).upper() for x in d["snp_ref"]])[sup]
    salt = np.array([str(x).upper() for x in d["snp_alt"]])[sup]
    sw = w[sup]
    pmin, pmax = int(spos.min()), int(spos.max())
    print(f"  predictor: {len(sup)} weight SNPs, span chr{chrom}:{pmin}-{pmax} ({(pmax-pmin)/1e6:.2f} Mb), cv_r2={float(d['sv_cv_r2'][ti]):.3f}")

    vcf, eur = eur_samples(chrom)
    eurfile = f"{TMP}/eur_{chrom}.txt"; open(eurfile, "w").write("\n".join(eur) + "\n")

    gp = read_geno_pos(vcf, eurfile, spos, chrom, "pred")
    pred = np.zeros(len(eur)); nused = 0
    for k in range(len(spos)):
        p = int(spos[k]); rr = sref[k]; aa = salt[k]
        if p in gp:
            vr, va, dos = gp[p]
            if {vr, va} == {rr, aa}:
                g_alt = (2 - dos) if (vr == aa and va == rr) else dos
                pred += sw[k] * g_alt; nused += 1
    pred = pred - pred.mean()
    print(f"  predictor dosage built from {nused}/{len(sup)} SNPs found biallelic in 1kGP-EUR (n_eur={len(eur)})")

    reg0 = int(region.split("_")[-2]); reg1 = int(region.split("_")[-1])
    lo = max(1, reg0 - MARGIN); hi = reg1 + MARGIN
    out = subprocess.run([TABIX, f"{GD}/{study}/{trait}.tsv.bgz", f"{chrom}:{lo}-{hi}"],
                         capture_output=True, text=True, check=True).stdout
    gws = []
    for ln in out.splitlines():
        f = ln.split("\t")
        try:
            pp = f[2].split("_"); rr, aa = pp[-2].upper(), pp[-1].upper()
            bb = float(f[5]); se = float(f[6])
            if se > 0:
                z = bb / se
                if abs(z) >= GWS: gws.append((int(f[1]), rr, aa, z))
        except Exception:
            pass
    print(f"  GWS trait SNPs (|z|>{GWS:.2f}, p<5e-8) in [{lo}-{hi}] ({(hi-lo)/1e6:.1f} Mb): {len(gws)}")
    if not gws:
        print("  >>> NO genome-wide-significant trait SNP anywhere in the window ->")
        print("  >>> the predictor CANNOT be re-tagging a known local GWS locus. VERDICT: DISTINCT (no known local signal).")
        continue
    gg = read_geno_pos(vcf, eurfile, [p for p, _, _, _ in gws], chrom, "gws")
    res = []
    for (p, rr, aa, z) in gws:
        if p in gg:
            vr, va, dos = gg[p]
            dd = dos - dos.mean()
            if dd.std() > 0 and pred.std() > 0:
                r = float(np.corrcoef(pred, dd)[0, 1]); res.append((r * r, p, rr, aa, z, r))
    res.sort(reverse=True)
    print(f"  GWS SNPs with usable 1kGP-EUR biallelic genotype: {len(res)}")
    if nused == 0 or not res:
        print("  >>> predictor unbuilt or no GWS SNP had a usable genotype -> cannot compute LD. (check inputs)")
        continue
    wset = set(int(x) for x in spos)
    print(f"  top GWS trait SNPs by r2 with predictor dosage (W = also a predictor weight SNP):")
    for r2, p, rr, aa, z, r in res[:8]:
        tags = ('IN-span' if pmin <= p <= pmax else 'flank') + (' W' if p in wset else '')
        print(f"    chr{chrom}:{p} {rr}/{aa}  trait_z={z:+6.1f}  r2={r2:.3f}  r={r:+.2f}  ({tags})")
    mx = res[0][0]
    res_nw = [t for t in res if t[1] not in wset]
    mx_nw = res_nw[0][0] if res_nw else 0.0
    n_selfsel = sum(1 for t in res if t[1] in wset)
    print(f"  >>> MAX r2 (ANY GWS SNP)        = {mx:.3f} at chr{chrom}:{res[0][1]} (z={res[0][4]:+.1f})")
    print(f"  >>> MAX r2 (NON-weight GWS SNP) = {mx_nw:.3f}" + (f" at chr{chrom}:{res_nw[0][1]} (z={res_nw[0][4]:+.1f})" if res_nw else " (none)"))
    print(f"  >>> {n_selfsel} of the region's GWS trait SNPs were SELECTED into this SV predictor (self-inclusion = strongest re-tagging)")
    verdict = ("RE-TAGS a known GWS locus (NOT a novel signal)" if mx_nw >= 0.5 else
               "DISTINCT: LD-independent of all local known GWS signals" if mx_nw < 0.1 else
               "PARTIAL / AMBIGUOUS independent LD (but check self-inclusion)")
    print(f"  >>> VERDICT (by non-weight LD): {verdict}")
print("\nDONE_LOCUS_SPEC")
