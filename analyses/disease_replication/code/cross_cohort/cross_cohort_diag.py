#!/usr/bin/env python
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
TMP = "/tmp/xcohort_diag"; os.makedirs(TMP, exist_ok=True)

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
            if len(f) < 3 or f[0] not in vs: continue
            if f[1] == "FIN": fin.append(f[0])
            if f[1] in EUR_POPS: eur.append(f[0])
    return fin, eur


def load_gwas(chrom, start, end, study, trait, pos, ref, alt):
    out = subprocess.run([TABIX, f"{GD}/{study}/{trait}.tsv.bgz", f"{chrom}:{start}-{end}"],
                         capture_output=True, text=True, check=True).stdout
    dz = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        try:
            pp = f[2].split("_"); r, a = pp[-2].upper(), pp[-1].upper()
            maf = float(f[3]); N = float(f[4]); bb = float(f[5]); se = float(f[6])
            if se > 0: dz[(int(f[1]), r, a)] = (bb / se, N, maf)
        except Exception:
            pass
    z = np.full(len(pos), np.nan); Nn = np.full(len(pos), np.nan); mf = np.full(len(pos), np.nan)
    for i in range(len(pos)):
        k = (int(pos[i]), ref[i], alt[i]); ks = (int(pos[i]), alt[i], ref[i])
        if k in dz: z[i], Nn[i], mf[i] = dz[k]
        elif ks in dz: v = dz[ks]; z[i], Nn[i], mf[i] = -v[0], v[1], v[2]
    return z, Nn, mf


def build_R_AF(sample_file, vcf, vcf_chrom, S_pos, S_ref, S_alt):
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
        k = posS[p]; wr, wa = S_ref[k].upper(), S_alt[k].upper()
        if {rr, aa} != {wr, wa}: continue
        flip = (rr == wa and aa == wr); gts = f[3:]; dos = np.empty(len(gts))
        for j, g in enumerate(gts):
            a0 = g[0]; a1 = g[2] if len(g) >= 3 else "."
            if a0 == "." or a1 == ".": dos[j] = np.nan
            else:
                c = (a0 == "1") + (a1 == "1"); dos[j] = (2 - c) if flip else c
        m = np.isnan(dos)
        if m.any(): dos[m] = np.nanmean(dos)
        G[k, :] = dos; matched += 1
    af = np.nanmean(G, axis=1) / 2.0
    sd = np.nanstd(G, axis=1); good = np.where(sd > 0)[0]; R = np.eye(len(S_pos))
    if len(good):
        Rg = np.corrcoef(np.nan_to_num(G[good, :]))
        for ii, gi in enumerate(good):
            for jj, gj in enumerate(good): R[gi, gj] = Rg[ii, jj]
    R = np.nan_to_num(R, nan=0.0); np.fill_diagonal(R, 1.0)
    return R, af, matched


def diag(b_s, R_ss):
    bRb = float(b_s @ (R_ss @ b_s)); ss2 = float(b_s @ b_s)
    vr = bRb / ss2 if ss2 > 0 else float("nan")
    ev = np.linalg.eigvalsh(R_ss); ev = ev[ev > 1e-10]
    er = float((ev.sum() ** 2) / np.sum(ev ** 2)) if len(ev) else 0.0
    return bRb, ss2, vr, math.sqrt(max(bRb, 1e-12)), er


def atomic_tsv(path, recs, cols):
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in recs:
            fh.write("\t".join("" if r.get(c) is None else
                               (f"{r[c]:.6g}" if isinstance(r.get(c), float) else str(r.get(c))) for c in cols) + "\n")
    os.replace(tmp, path)


rows = []
for region, sv_name, label, cohorts in CANDS:
    parts = region.split("_"); chrom = parts[0].replace("chr", ""); start = int(parts[-2]); end = int(parts[-1])
    vcf = f"{RES}/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"; vcf_chrom = f"chr{chrom}"
    b = np.load(f"{SSIMP}/output/bundles/chr{chrom}/{region}_bundle.npz", allow_pickle=True)
    R_bundle = b["R_snp_snp"].astype(float); B = b["B"].astype(float); pos = b["snp_pos"].astype(int)
    ref = np.array([str(x).upper() for x in b["snp_ref"]]); alt = np.array([str(x).upper() for x in b["snp_alt"]])
    svn = np.array([str(x) for x in b["sv_locus"]]); ti = int(np.where(svn == sv_name)[0][0]); cvr2 = float(b["sv_cv_r2"][ti])
    bw = B[:, ti]; full_sup = np.where(bw != 0)[0]
    _, _, vr_full, D_full, er_full = diag(bw[full_sup], R_bundle[np.ix_(full_sup, full_sup)])
    print(f"\n######## {label}: {sv_name} cv_r2={cvr2:.3f} ########")
    print(f"  [predictor intrinsic | full support n={len(full_sup)}] var_ratio={vr_full:.4f} D={D_full:.4f} "
          f"eff_rank={er_full:.1f}/{len(full_sup)}  {'<< DEGENERATE (<0.10)' if vr_full < 0.10 else 'OK (>=0.10)'}")
    order = full_sup[np.argsort(-np.abs(bw[full_sup]))]; tags = order[:min(20, len(order))]
    fin, eur = pop_samples(vcf)
    open(f"{TMP}/fin.txt", "w").write("\n".join(fin) + "\n"); open(f"{TMP}/eur.txt", "w").write("\n".join(eur) + "\n")
    _, fin_af, _ = build_R_AF(f"{TMP}/fin.txt", vcf, vcf_chrom, pos[tags], ref[tags], alt[tags])
    _, eur_af, _ = build_R_AF(f"{TMP}/eur.txt", vcf, vcf_chrom, pos[tags], ref[tags], alt[tags])
    daf = float(np.nanmean(np.abs(fin_af - eur_af)))
    print(f"  [top-{len(tags)} tags | 1kGP AF] mean|FIN-EUR dAF|={daf:.3f}  FIN~{np.nanmean(fin_af):.3f} EUR~{np.nanmean(eur_af):.3f}")
    for study, trait, anc in cohorts:
        z, Nn, mf = load_gwas(chrom, start, end, study, trait, pos, ref, alt)
        obs = np.where(np.isfinite(z) & (bw != 0))[0]
        if len(obs) == 0:
            print(f"  [{study}/{trait}] no observed support"); continue
        b_s = bw[obs]; z_s = z[obs]
        bRb, ss2, vr, D, er = diag(b_s, R_bundle[np.ix_(obs, obs)])
        num = float(b_s @ z_s); zsv = num / D
        sf = f"{TMP}/{'fin' if anc == 'FIN' else 'eur'}.txt"
        Rm, _, nm = build_R_AF(sf, vcf, vcf_chrom, pos[obs], ref[obs], alt[obs])
        _, _, vrm, Dm, erm = diag(b_s, Rm); zsvm = num / Dm
        terms = b_s * z_s; top = int(np.argmax(np.abs(terms))); topfrac = abs(terms[top]) / max(abs(num), 1e-9)
        zsv_lo = (num - terms[top]) / D; dztop = abs(zsv - zsv_lo) / max(abs(zsv), 1e-9)
        Nmed = float(np.nanmedian(Nn[obs]))
        flag = "DEGEN" if vr < 0.10 else "ok"
        print(f"  [{study:7s}/{trait:11s}] n_sup={len(obs):3d} N_med={Nmed:,.0f} "
              f"var_ratio(bundle)={vr:.4f}/{flag} var_ratio(matched)={vrm:.4f} eff_rank={er:.1f}/{len(obs)} "
              f"z_SV={zsv:+.2f}(m{zsvm:+.2f}) num={num:+.2f} top_tag_frac={topfrac:.2f} dz_drop_top={dztop:.2f}")
        rows.append(dict(sv=sv_name, label=label, cv_r2=cvr2, study=study, trait=trait, anc=anc,
                         n_support=len(obs), N_median=Nmed, vr_full=vr_full, var_ratio_bundle=vr,
                         var_ratio_matched=vrm, eff_rank=er, z_sv_bundle=zsv, z_sv_matched=zsvm, num=num,
                         top_tag_frac=topfrac, dz_drop_top=dztop, mean_dAF_FIN_EUR=daf))

if rows:
    atomic_tsv(f"{OUT}/xcohort_diag.tsv", rows,
               ["sv", "label", "cv_r2", "study", "trait", "anc", "n_support", "N_median", "vr_full",
                "var_ratio_bundle", "var_ratio_matched", "eff_rank", "z_sv_bundle", "z_sv_matched",
                "num", "top_tag_frac", "dz_drop_top", "mean_dAF_FIN_EUR"])
print(f"\nwrote {OUT}/xcohort_diag.tsv")
print("GATE: var_ratio<0.10 => z_SV is degeneracy-inflated (BOTH legs) => sign-concordance is a shared-b artifact, not biology.")
print("N: compare N_median FinnGen vs UKB -> which leg is actually better-powered (do NOT assume).")
print("CONC: top_tag_frac near 1 / dz_drop_top large => numerator driven by ONE tag (tag-driven, not distributed).")
print("AF: large mean|FIN-EUR dAF| => tags are population-differentiated (enrichment/stratification candidate).")
print("DONE_XCOHORT_DIAG")
