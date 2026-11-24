#!/usr/bin/env python
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
GWS = 5.4513; MARGIN = 1_000_000; MAXLEADS = 12; PRUNE = 0.1
TMP = "/tmp/export_leads"; os.makedirs(TMP, exist_ok=True)
OUT = f"{P}/Project1_Centromere/analyses/disease_replication/aou_handoff/leads_for_aou_conditional.tsv"

SVS = [
  ("chr3_sv20_DEL_non_best6_alt_edge1000",  "chr3_super1_85533081_94536472",  "3",  "neale",   "21002_irnt", "body_weight"),
  ("chr4_sv349_DEL_non_best8_fksize200k",   "chr4_super1_47409879_53011266",  "4",  "neale",   "50_irnt",    "height"),
  ("chr11_sv263_DEL_non_best8_edge500",     "chr11_super1_49844498_54927998", "11", "finngen", "BMI_IRN",    "BMI"),
  ("chr4_sv398_INS_non_best8_fksize100k",   "chr4_super1_47409879_53011266",  "4",  "neale",   "30110_irnt", "PDW_platelet"),
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
    cmd = (f"{BCF} view -S {eurfile} -R {tf} -v snps -m2 -M2 {vcf} | {BCF} query -f '%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    g = {}
    for ln in out.splitlines():
        f = ln.split("\t"); p = int(f[0]); rsid = f[1]; r = f[2].upper(); a = f[3].upper(); gts = f[4:]
        dos = np.array([(gg[0] == "1") + (gg[2] == "1") if (len(gg) >= 3 and gg[0] != "." and gg[2] != ".") else np.nan for gg in gts], float)
        m = np.isnan(dos)
        if m.all(): continue
        dos[m] = np.nanmean(dos); g[p] = (r, a, dos, rsid)
    return g

def zmap(study, trait, chrom, lo, hi):
    out = subprocess.run([TABIX, f"{GD}/{study}/{trait}.tsv.bgz", f"{chrom}:{lo}-{hi}"], capture_output=True, text=True, check=True).stdout
    dz = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        try:
            pp = f[2].split("_"); rr, aa = pp[-2].upper(), pp[-1].upper()
            se = float(f[6]); bb = float(f[5])
            if se > 0: dz[(int(f[1]), rr, aa)] = bb / se
        except Exception: pass
    return dz

rows = []
for svid, region, chrom, study, trait, htrait in SVS:
    d = np.load(f"{BUND}/chr{chrom}/{region}_bundle.npz", allow_pickle=True)
    svn = np.array([str(x) for x in d["sv_locus"]]); ti = int(np.where(svn == svid)[0][0])
    w = d["B"][:, ti].astype(float); sup = np.where(w != 0)[0]
    spos = d["snp_pos"].astype(int)[sup]
    sref = np.array([str(x).upper() for x in d["snp_ref"]])[sup]
    salt = np.array([str(x).upper() for x in d["snp_alt"]])[sup]; sw = w[sup]
    wset = set(int(x) for x in spos)
    vcf, eur = eur_samples(chrom); eurfile = f"{TMP}/eur_{chrom}.txt"; open(eurfile, "w").write("\n".join(eur) + "\n")
    gp = read_geno_pos(vcf, eurfile, spos, chrom, "pred")
    pred = np.zeros(len(eur))
    for k in range(len(spos)):
        p = int(spos[k]); rr = sref[k]; aa = salt[k]
        if p in gp:
            vr, va, dos, _ = gp[p]
            if {vr, va} == {rr, aa}:
                pred += sw[k] * ((2 - dos) if (vr == aa and va == rr) else dos)
    pred = pred - pred.mean()
    reg0 = int(region.split("_")[-2]); reg1 = int(region.split("_")[-1])
    dz = zmap(study, trait, chrom, max(1, reg0 - MARGIN), reg1 + MARGIN)
    gws = [(p, rr, aa, z) for (p, rr, aa), z in dz.items() if abs(z) >= GWS]
    gg = read_geno_pos(vcf, eurfile, [p for p, _, _, _ in gws], chrom, "gws")
    cand = []
    for (p, rr, aa, zl) in gws:
        if p not in gg: continue
        vr, va, dos, rsid = gg[p]
        if {vr, va} != {rr, aa}: continue
        L = (2 - dos) if (vr == aa and va == rr) else dos
        Lc = L - L.mean()
        if Lc.std() == 0 or pred.std() == 0: continue
        r = float(np.corrcoef(pred, Lc)[0, 1])
        cand.append((r * r, r, p, rr, aa, zl, rsid, p in wset, Lc))
    cand.sort(key=lambda c: -c[0])
    leads = []
    for c in cand:
        if all(float(np.corrcoef(c[8], s[8])[0, 1]) ** 2 < PRUNE for s in leads): leads.append(c)
        if len(leads) >= MAXLEADS: break
    print(f"{svid} [{htrait}, {study}]: {len(sup)} weight SNPs, {len(gws)} GWS in region -> {len(leads)} independent leads")
    for c in leads:
        r2, r, p, rr, aa, zl, rsid, inpred, _ = c
        rows.append([svid, htrait, study, trait, f"chr{chrom}", p, rsid, rr, aa, f"{zl:+.2f}", f"{r2:.3f}", str(inpred)])

hdr = ["sv_locus", "trait", "study", "trait_code", "chrom", "pos_hg38", "rsid", "ref", "alt_effect", "gwas_z", "r2_with_predictor", "is_predictor_weight_snp"]
with open(OUT, "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
print(f"\nWROTE {len(rows)} lead SNPs across {len(SVS)} SV-traits -> {OUT}")
print("DONE_EXPORT_LEADS")
