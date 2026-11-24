#!/usr/bin/env python
"""Cross-cohort (FinnGen <-> UKB/Neale) replication meta for the phewide UNRESOLVED survivors."""
import numpy as np, subprocess, os, math
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P = CENSV_WORKSPACE + ""
SSIMP = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
GD = f"{SSIMP}/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
BT = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
TABIX = f"{BT}/tabix"; BCFTOOLS = f"{BT}/bcftools"
RES = f"{P}/Project1_Centromere/resources/1000G"
PANEL = f"{P}/_shared/datasets/1KGP/v_canonical_phase3_20130502/EUR/integrated_call_samples_v3.20130502.ALL.panel"
EUR_POPS = {"FIN", "CEU", "GBR", "TSI", "IBS"}
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(f"{HERE}/../output/xcohort"); os.makedirs(OUT, exist_ok=True)
TMP = "/tmp/xcohort_meta"; os.makedirs(TMP, exist_ok=True)

CANDS = [
    ("chr12_super1_32012658_38833367", "chr12_sv721_DEL_non_best6_alt_edge1000", "height",
        [("finngen", "HEIGHT_IRN", "FIN"), ("neale", "50_irnt", "EUR")]),
    ("chr11_super1_49844498_54927998", "chr11_sv263_DEL_non_best8_edge500", "BMI",
        [("finngen", "BMI_IRN", "FIN"), ("neale", "21001_irnt", "EUR")]),
    ("chr11_super2_55315181_58690022", "chr11_sv271_DEL_ori_best8_edge500", "albumin",
        [("neale", "30600_irnt", "EUR")]),
]


def pop_samples(vcf):
    vs = set(subprocess.run([BCFTOOLS, "query", "-l", vcf], capture_output=True, text=True, check=True).stdout.split())
    fin, eur = [], []
    with open(PANEL) as fh:
        next(fh)
        for ln in fh:
            f = ln.split()
            if len(f) < 3 or f[0] not in vs:
                continue
            if f[1] == "FIN": fin.append(f[0])
            if f[1] in EUR_POPS: eur.append(f[0])
    return fin, eur


def load_z(chrom, start, end, study, trait, pos, ref, alt):
    out = subprocess.run([TABIX, f"{GD}/{study}/{trait}.tsv.bgz", f"{chrom}:{start}-{end}"],
                         capture_output=True, text=True, check=True).stdout
    dz = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        try:
            pp = f[2].split("_"); r, a = pp[-2].upper(), pp[-1].upper(); se = float(f[6]); bb = float(f[5])
            if se > 0: dz[(int(f[1]), r, a)] = bb / se
        except Exception:
            pass
    z = np.full(len(pos), np.nan)
    for i in range(len(pos)):
        k = (int(pos[i]), ref[i], alt[i]); ks = (int(pos[i]), alt[i], ref[i])
        if k in dz: z[i] = dz[k]
        elif ks in dz: z[i] = -dz[ks]
    return z


def build_R(sample_file, vcf, vcf_chrom, S_pos, S_ref, S_alt):
    posS = {int(p): k for k, p in enumerate(S_pos)}
    with open(f"{TMP}/regions.txt", "w") as fh:
        for p in sorted(set(int(x) for x in S_pos)):
            fh.write(f"{vcf_chrom}\t{p}\t{p}\n")
    cmd = (f"{BCFTOOLS} view -S {sample_file} -R {TMP}/regions.txt -v snps -m2 -M2 {vcf} "
           f"| {BCFTOOLS} query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    n = sum(1 for _ in open(sample_file))
    G = np.full((len(S_pos), n), np.nan); matched = 0
    for ln in out.splitlines():
        f = ln.split("\t"); p = int(f[0]); rr = f[1].upper(); aa = f[2].upper()
        if p not in posS: continue
        k = posS[p]; want_r, want_a = S_ref[k].upper(), S_alt[k].upper()
        if {rr, aa} != {want_r, want_a}: continue
        flip = (rr == want_a and aa == want_r)
        gts = f[3:]; dos = np.empty(len(gts))
        for j, g in enumerate(gts):
            a0 = g[0]; a1 = g[2] if len(g) >= 3 else "."
            if a0 == "." or a1 == ".": dos[j] = np.nan
            else:
                c = (a0 == "1") + (a1 == "1"); dos[j] = (2 - c) if flip else c
        m = np.isnan(dos)
        if m.any(): dos[m] = np.nanmean(dos)
        G[k, :] = dos; matched += 1
    sd = np.nanstd(G, axis=1); good = np.where(sd > 0)[0]
    R = np.eye(len(S_pos))
    if len(good):
        Rg = np.corrcoef(np.nan_to_num(G[good, :]))
        for ii, gi in enumerate(good):
            for jj, gj in enumerate(good):
                R[gi, gj] = Rg[ii, jj]
    R = np.nan_to_num(R, nan=0.0); np.fill_diagonal(R, 1.0)
    return R, matched


def zsv(b_s, z_s, R_ss):
    num = float(b_s @ z_s)
    D = float(np.sqrt(max(b_s @ (R_ss @ b_s), 1e-12)))
    return num, D, num / D


def p2(x):
    return math.erfc(abs(x) / math.sqrt(2.0))


def atomic_tsv(path, recs, cols):
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in recs:
            fh.write("\t".join(
                "" if r.get(c) is None else (f"{r[c]:.6g}" if isinstance(r.get(c), float) else str(r.get(c)))
                for c in cols) + "\n")
    os.replace(tmp, path)


rows, summary = [], []
for region, sv_name, label, cohorts in CANDS:
    parts = region.split("_"); chrom = parts[0].replace("chr", ""); start = int(parts[-2]); end = int(parts[-1])
    vcf = f"{RES}/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    vcf_chrom = f"chr{chrom}"
    bundle = f"{SSIMP}/output/bundles/chr{chrom}/{region}_bundle.npz"
    b = np.load(bundle, allow_pickle=True)
    R_bundle = b["R_snp_snp"].astype(float); B = b["B"].astype(float)
    pos = b["snp_pos"].astype(int)
    ref = np.array([str(x).upper() for x in b["snp_ref"]]); alt = np.array([str(x).upper() for x in b["snp_alt"]])
    sv = np.array([str(x) for x in b["sv_locus"]]); ti = int(np.where(sv == sv_name)[0][0])
    cvr2 = float(b["sv_cv_r2"][ti]) if "sv_cv_r2" in b else float("nan")
    print(f"\n######## {label}: {region} | {sv_name} | cv_r2={cvr2:.3f} ########")
    fin = eur = None; leg = {}
    for study, trait, anc in cohorts:
        z = load_z(chrom, start, end, study, trait, pos, ref, alt)
        obs = np.where(np.isfinite(z) & (B[:, ti] != 0))[0]
        if len(obs) == 0:
            print(f"  [{study}/{trait}] no observed support; skipping"); continue
        b_s = B[obs, ti]; z_s = z[obs]
        num, Db, zb = zsv(b_s, z_s, R_bundle[np.ix_(obs, obs)])
        if fin is None:
            fin, eur = pop_samples(vcf)
            open(f"{TMP}/fin.txt", "w").write("\n".join(fin) + "\n")
            open(f"{TMP}/eur.txt", "w").write("\n".join(eur) + "\n")
        sf = f"{TMP}/{'fin' if anc == 'FIN' else 'eur'}.txt"
        Rm, nmatch = build_R(sf, vcf, vcf_chrom, pos[obs], ref[obs], alt[obs])
        num_m, Dm, zm = zsv(b_s, z_s, Rm)
        npop = len(fin) if anc == "FIN" else len(eur)
        print(f"  [{study:7s}/{trait:11s} anc={anc} n_pop={npop:3d}] n_support={len(obs):4d}  "
              f"num=b'z={num:+.2f}  z_SV(bundle)={zb:+.3f}  z_SV(matched {nmatch}/{len(obs)})={zm:+.3f}")
        leg[study] = dict(sv=sv_name, region=region, label=label, cv_r2=cvr2, study=study, trait=trait, anc=anc,
                          n_support=int(len(obs)), n_pop=int(npop), num=num, D_bundle=Db, z_sv_bundle=zb,
                          D_matched=Dm, n_ld_matched=int(nmatch), z_sv_matched=zm)
        rows.append(leg[study])
    if len(leg) == 2:
        (s1, l1), (s2, l2) = list(leg.items())
        z1, z2 = l1["z_sv_matched"], l2["z_sv_matched"]
        sign_conc = bool(np.sign(z1) == np.sign(z2))
        Q = (z1 - z2) ** 2 / 2.0
        I2 = max(0.0, (Q - 1.0) / Q) if Q > 0 else 0.0
        het_p = p2((z1 - z2) / math.sqrt(2.0))
        z_meta = (z1 + z2) / math.sqrt(2.0); meta_p = p2(z_meta)
        verdict = ("REPLICATES(sign+homog)" if (sign_conc and het_p > 0.05)
                   else "SIGN-CONCORDANT(heterogeneous)" if sign_conc else "DISCORDANT(sign-flip)")
        print(f"  ==> XCOHORT [{s1} z_SV={z1:+.2f}(num {l1['num']:+.1f}) | {s2} z_SV={z2:+.2f}(num {l2['num']:+.1f})] "
              f"sign_conc={sign_conc} Q={Q:.2f} I2={I2:.2f} het_p={het_p:.3g} meta_z={z_meta:+.2f} meta_p={meta_p:.3g} => {verdict}")
        summary.append(dict(sv=sv_name, label=label, cv_r2=cvr2, cohort1=s1, z_sv_1=z1, num_1=l1["num"],
                            cohort2=s2, z_sv_2=z2, num_2=l2["num"], sign_concordant=sign_conc,
                            Q=Q, I2=I2, het_p=het_p, meta_z=z_meta, meta_p=meta_p, verdict=verdict))
    else:
        print(f"  ==> single-cohort only ({label}); independent cross-cohort needs an external GWAS")
        summary.append(dict(sv=sv_name, label=label, cv_r2=cvr2, verdict="SINGLE-COHORT(no FinnGen counterpart)"))

if rows:
    atomic_tsv(f"{OUT}/xcohort_per_cohort.tsv", rows,
               ["sv", "region", "label", "cv_r2", "study", "trait", "anc", "n_support", "n_pop",
                "num", "D_bundle", "z_sv_bundle", "D_matched", "n_ld_matched", "z_sv_matched"])
if summary:
    atomic_tsv(f"{OUT}/xcohort_summary.tsv", summary,
               ["sv", "label", "cv_r2", "cohort1", "z_sv_1", "num_1", "cohort2", "z_sv_2", "num_2",
                "sign_concordant", "Q", "I2", "het_p", "meta_z", "meta_p", "verdict"])
print(f"\nwrote {OUT}/xcohort_per_cohort.tsv + xcohort_summary.tsv")
print("READ: sign-concordant + het_p>0.05 => replicates across INDEPENDENT cohorts (numerator b'z is the independent signal).")
print("CAVEAT: shared elastic-net weights b => a predictor-baked artifact would also replicate; concordance != SV-specificity (cv_r2/WGS).")
print("DONE_XCOHORT_META")
