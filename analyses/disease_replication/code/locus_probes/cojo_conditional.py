#!/usr/bin/env python
import numpy as np, subprocess, os
from math import erf, sqrt as msqrt
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
GWS   = 5.4513; MARGIN = 1_000_000
TMP   = "/tmp/cojo_cond"; os.makedirs(TMP, exist_ok=True)

SVS = [
  ("chr3_sv20_DEL_non_best6_alt_edge1000",  "chr3_super1_85533081_94536472",  "3",  "neale",   "21002_irnt", -0.7065, "body weight"),
  ("chr4_sv349_DEL_non_best8_fksize200k",   "chr4_super1_47409879_53011266",  "4",  "neale",   "50_irnt",    -1.0746, "height"),
  ("chr11_sv263_DEL_non_best8_edge500",     "chr11_super1_49844498_54927998", "11", "finngen", "BMI_IRN",    -1.0308, "BMI"),
  ("chr4_sv398_INS_non_best8_fksize100k",   "chr4_super1_47409879_53011266",  "4",  "neale",   "30110_irnt", -0.4563, "PDW/platelet (nominal AoU survivor, novelty untested)"),
  ("chr12_sv721_DEL_non_best6_alt_edge1000","chr12_super1_32012658_38833367", "12", "finngen", "HEIGHT_IRN", +0.5115, "height (non-replicated ctrl)"),
]

def p_from_z(z):
    return 2.0 * (1.0 - 0.5 * (1.0 + erf(abs(z) / msqrt(2.0))))

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
    cmd = (f"{BCF} view -S {eurfile} -R {tf} -v snps -m2 -M2 {vcf} | {BCF} query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    g = {}
    for ln in out.splitlines():
        f = ln.split("\t"); p = int(f[0]); r = f[1].upper(); a = f[2].upper(); gts = f[3:]
        dos = np.array([(gg[0] == "1") + (gg[2] == "1") if (len(gg) >= 3 and gg[0] != "." and gg[2] != ".") else np.nan for gg in gts], float)
        m = np.isnan(dos)
        if m.all(): continue
        dos[m] = np.nanmean(dos); g[p] = (r, a, dos)
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

for svid, region, chrom, study, trait, disc_bz, tlabel in SVS:
    print(f"\n{'='*88}\n{svid}  [{tlabel}]\n{'='*88}")
    d = np.load(f"{BUND}/chr{chrom}/{region}_bundle.npz", allow_pickle=True)
    svn = np.array([str(x) for x in d["sv_locus"]]); ti = int(np.where(svn == svid)[0][0])
    w = d["B"][:, ti].astype(float); sup = np.where(w != 0)[0]
    spos = d["snp_pos"].astype(int)[sup]
    sref = np.array([str(x).upper() for x in d["snp_ref"]])[sup]
    salt = np.array([str(x).upper() for x in d["snp_alt"]])[sup]; sw = w[sup]
    vcf, eur = eur_samples(chrom); eurfile = f"{TMP}/eur_{chrom}.txt"; open(eurfile, "w").write("\n".join(eur) + "\n")
    gp = read_geno_pos(vcf, eurfile, spos, chrom, "pred")

    reg0 = int(region.split("_")[-2]); reg1 = int(region.split("_")[-1])
    dz = zmap(study, trait, chrom, max(1, reg0 - MARGIN), reg1 + MARGIN)

    rows = []
    for k in range(len(spos)):
        p = int(spos[k]); rr = sref[k]; aa = salt[k]
        if p not in gp: continue
        vr, va, dos = gp[p]
        if {vr, va} != {rr, aa}: continue
        g_alt = (2 - dos) if (vr == aa and va == rr) else dos
        z = dz.get((p, rr, aa))
        if z is None:
            z2 = dz.get((p, aa, rr)); z = (-z2) if z2 is not None else None
        if z is None or not np.isfinite(z): continue
        if np.std(g_alt) == 0: continue
        rows.append((p, sw[k], z, g_alt))
    if len(rows) < 3:
        print(f"  too few usable weight SNPs ({len(rows)}) -> skip"); continue
    wu = np.array([r[1] for r in rows]); zu = np.array([r[2] for r in rows])
    G = np.vstack([r[3] for r in rows])
    R = np.corrcoef(G)
    wRw = float(wu @ R @ wu); wz = float(wu @ zu)
    Z_SV = wz / msqrt(wRw)
    D = wu @ G
    wset = set(int(r[0]) for r in rows)
    print(f"  usable weight SNPs={len(rows)}  w'z={wz:+.4f} (stored discovery_bz={disc_bz:+.4f})  Z_SV(TWAS)={Z_SV:+.2f}  p={p_from_z(Z_SV):.1e}")

    gws = [(p, rr, aa, z) for (p, rr, aa), z in dz.items() if abs(z) >= GWS]
    gg = read_geno_pos(vcf, eurfile, [p for p, _, _, _ in gws], chrom, "gws")
    cand = []
    Dc = D - D.mean()
    for (p, rr, aa, zl) in gws:
        if p not in gg: continue
        vr, va, dos = gg[p]
        if {vr, va} != {rr, aa}: continue
        L = (2 - dos) if (vr == aa and va == rr) else dos
        Lc = L - L.mean()
        if Lc.std() == 0: continue
        rr_ = float(np.corrcoef(Dc, Lc)[0, 1])
        cand.append((rr_ * rr_, rr_, p, zl, p in wset, Lc))
    cand.sort(key=lambda c: -c[0])
    def cond1(r2, r, zl):
        if 1 - r2 <= 1e-9: return None
        return (Z_SV - r * zl) / msqrt(1 - r2)
    top_any = cand[0]
    top_nw = next((c for c in cand if not c[4]), None)
    for tagname, lead in [("TOP-ANY", top_any), ("TOP-NONWEIGHT", top_nw)]:
        if lead is None: print(f"  {tagname}: none"); continue
        r2, r, p, zl, isw, _ = lead; zc = cond1(r2, r, zl)
        print(f"  1-SNP cond on {tagname} chr{chrom}:{p} (z_lead={zl:+.1f}, r={r:+.2f}, r2={r2:.2f}{', SELF-WT' if isw else ''}): "
              f"Z_SV|lead={zc:+.2f} p={p_from_z(zc):.1e}  ({abs(zc)/abs(Z_SV)*100:.0f}% retained)")

    by_sig = sorted(cand, key=lambda c: -abs(c[3]))
    S = []
    for c in by_sig:
        if all(float(np.corrcoef(c[5], s[5])[0, 1]) ** 2 < 0.1 for s in S):
            S.append(c)
        if len(S) >= 10: break
    if S:
        M = np.vstack([s[5] for s in S])
        R_SS = np.corrcoef(M)
        if np.ndim(R_SS) == 0: R_SS = np.array([[1.0]])
        Z_S = np.array([s[3] for s in S])
        r_vS = np.array([float(np.corrcoef(Dc, s[5])[0, 1]) for s in S])
        try:
            a = np.linalg.solve(R_SS, Z_S); b = np.linalg.solve(R_SS, r_vS)
            denom2 = 1.0 - float(r_vS @ b)
            zc = (Z_SV - float(r_vS @ a)) / msqrt(denom2) if denom2 > 1e-6 else float("nan")
            R2_expl = float(r_vS @ b)
            print(f"  MULTI-SNP cond on {len(S)} independent GWS leads (r2-pruned<0.1): "
                  f"Z_SV|S={zc:+.2f} p={p_from_z(zc):.1e}  ({abs(zc)/abs(Z_SV)*100:.0f}% retained; leads explain {R2_expl*100:.0f}% of predictor var)")
            verdict = ("EXPLAINED by the known-locus block (signal collapses -> NOT novel)" if abs(zc) < 2 else
                       "PARTIALLY attenuated, residual survives full block conditioning" if abs(zc) < abs(Z_SV) * 0.7 else
                       "SURVIVES full block conditioning (independent of local known GWS signals)")
            print(f"  >>> VERDICT (multi-SNP, significance-pruned): {verdict}")
        except np.linalg.LinAlgError:
            print(f"  MULTI-SNP: R_SS singular (m={len(S)}) -> skip")

    by_r = sorted(cand, key=lambda c: -c[0])
    Sp = []
    for c in by_r:
        if all(float(np.corrcoef(c[5], s[5])[0, 1]) ** 2 < 0.1 for s in Sp): Sp.append(c)
        if len(Sp) >= 10: break
    if Sp:
        Mp = np.vstack([s[5] for s in Sp]); Rp = np.corrcoef(Mp)
        if np.ndim(Rp) == 0: Rp = np.array([[1.0]])
        Zp = np.array([s[3] for s in Sp]); rp = np.array([float(np.corrcoef(Dc, s[5])[0, 1]) for s in Sp])
        try:
            ap = np.linalg.solve(Rp, Zp); bp = np.linalg.solve(Rp, rp); den2 = 1.0 - float(rp @ bp)
            zc = (Z_SV - float(rp @ ap)) / msqrt(den2) if den2 > 1e-6 else float("nan")
            print(f"  PREDICTOR-AWARE cond on {len(Sp)} SNPs (top-|r|, r2-pruned<0.1): "
                  f"Z_SV|S={zc:+.2f} p={p_from_z(zc):.1e}  ({abs(zc)/abs(Z_SV)*100:.0f}% retained; these explain {float(rp@bp)*100:.0f}% of predictor var)")
            v2 = ("EXPLAINED by what it tags (collapses -> overlaps known locus, not novel)" if abs(zc) < 2 else
                  "PARTIALLY attenuated" if abs(zc) < abs(Z_SV) * 0.7 else "still SURVIVES")
            print(f"  >>> VERDICT (predictor-aware): {v2}")
        except np.linalg.LinAlgError:
            print(f"  PREDICTOR-AWARE: R singular -> skip")

    by_r_nw = [c for c in sorted(cand, key=lambda c: -c[0]) if not c[4]]
    Sq = []
    for c in by_r_nw:
        if all(float(np.corrcoef(c[5], s[5])[0, 1]) ** 2 < 0.1 for s in Sq): Sq.append(c)
        if len(Sq) >= 10: break
    if Sq:
        Mq = np.vstack([s[5] for s in Sq]); Rq = np.corrcoef(Mq)
        if np.ndim(Rq) == 0: Rq = np.array([[1.0]])
        Zq = np.array([s[3] for s in Sq]); rq = np.array([float(np.corrcoef(Dc, s[5])[0, 1]) for s in Sq])
        try:
            aq = np.linalg.solve(Rq, Zq); bq = np.linalg.solve(Rq, rq); dq = 1.0 - float(rq @ bq)
            zc = (Z_SV - float(rq @ aq)) / msqrt(dq) if dq > 1e-6 else float("nan")
            atten = (1 - abs(zc) / abs(Z_SV)) * 100
            v3 = ("SURVIVES (SV-specific vs non-predictor regional SNPs)" if abs(zc) >= abs(Z_SV) * 0.7 else
                  "PARTIALLY attenuated" if abs(zc) >= 2 else "COLLAPSES (regional tag)")
            print(f"  [AoU-MATCHED excl-predictor-positions] {len(Sq)} non-weight SNPs: Z_SV|S={zc:+.2f} "
                  f"p={p_from_z(zc):.1e} ({atten:.0f}% attenuated) -> {v3}")
        except np.linalg.LinAlgError:
            print(f"  [AoU-MATCHED]: singular")
print("\nDONE_COJO")
