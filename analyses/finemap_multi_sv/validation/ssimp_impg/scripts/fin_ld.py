#!/usr/bin/env python
import numpy as np, subprocess, os, sys
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P = CENSV_WORKSPACE + ""
BASE = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
BUNDLE = f"{BASE}/output/bundles/chr3/chr3_super1_85533081_94536472_bundle.npz"
GD = f"{BASE}/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"
BT = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
TABIX = f"{BT}/tabix"; BCFTOOLS = f"{BT}/bcftools"
VCF = f"{P}/Project1_Centromere/resources/1000G/1kGP_high_coverage_Illumina.chr3.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PANEL = f"{P}/_shared/datasets/1KGP/v_canonical_phase3_20130502/EUR/integrated_call_samples_v3.20130502.ALL.panel"
CHROM, START, END = "3", 85533081, 94536472
VCF_CHROM = "chr3"
TARGET = "chr3_sv40_DEL_non_best6_alt_edge500"
TMP = "/tmp/fin_ld_redo"; os.makedirs(TMP, exist_ok=True)

genes = {"PROS1": (93873038, 93980003), "ARL13B": (93980128, 94055678), "DHFR2": (93981180, 94063389),
         "STX19": (94012043, 94028658), "NSUN3": (94062980, 94131832), "HTR1F": (87792548, 87993839)}
def gene_of(p):
    for g, (s, e) in genes.items():
        if s <= p <= e: return g
    return "prox(87)" if p < 88500000 else ("distPROS1nbhd" if p >= 93000000 else "interior")

b = np.load(BUNDLE, allow_pickle=True)
R_bundle = b["R_snp_snp"].astype(float); B = b["B"].astype(float)
pos = b["snp_pos"].astype(int)
ref = np.array([str(x).upper() for x in b["snp_ref"]]); alt = np.array([str(x).upper() for x in b["snp_alt"]])
sv = np.array([str(x) for x in b["sv_locus"]])
ti = int(np.where(sv == TARGET)[0][0])

def load_z(study, trait):
    out = subprocess.run([TABIX, f"{GD}/{study}/{trait}.tsv.bgz", f"{CHROM}:{START}-{END}"],
                         capture_output=True, text=True, check=True)
    dz = {}
    for ln in out.stdout.splitlines():
        f = ln.split("\t")
        try:
            parts = f[2].split("_"); r, a = parts[-2].upper(), parts[-1].upper()
            se = float(f[6]); bb = float(f[5])
            if se > 0: dz[(int(f[1]), r, a)] = bb / se
        except Exception: pass
    z = np.full(len(pos), np.nan)
    for i in range(len(pos)):
        k = (int(pos[i]), ref[i], alt[i]); ks = (int(pos[i]), alt[i], ref[i])
        if k in dz: z[i] = dz[k]
        elif ks in dz: z[i] = -dz[ks]
    return z

z = load_z("finngen", "I9_VTE")
obs = np.where(np.isfinite(z))[0]
zo = z[obs].copy(); poso = pos[obs]; refo = ref[obs]; alto = alt[obs]; bo = B[obs, :][:, ti]

def lead_in(lo, hi):
    m = (poso >= lo) & (poso < hi)
    if not m.any(): return None
    idx = np.where(m)[0]
    return idx[int(np.argmax(np.abs(zo[idx])))]
LEADS = {"interior(88.5-93.0)": lead_in(88500000, 93000000),
         "PROS1nbhd(93.8-94.4)": lead_in(93800000, 94400000),
         "prox(87.5-88.0)": lead_in(87500000, 88000000)}

tag_idx = np.where(bo != 0)[0]
lead_idx = [j for j in LEADS.values() if j is not None]
S_idx = sorted(set(tag_idx.tolist()) | set(lead_idx))
S_pos = poso[S_idx]; S_ref = refo[S_idx]; S_alt = alto[S_idx]
posS = {int(p): k for k, p in enumerate(S_pos)}
print("support obs=%d/%d ; sv40 tags=%d ; leads=%d ; |S|=%d" %
      (len(obs), len(pos), len(tag_idx), len(lead_idx), len(S_idx)))
print("LD-FREE invariants (same for every panel):")
contrib = bo * zo
isprox = np.array([gene_of(p) in ("prox(87)", "HTR1F") for p in poso])
isdist = np.array([gene_of(p) in ("distPROS1nbhd","PROS1","ARL13B","DHFR2","STX19","NSUN3") for p in poso])
print("  beta'z total=%.3f  proximal=%.3f (%.0f%%)  distal/PROS1nbhd=%.3f (%.0f%%)" %
      (contrib.sum(), contrib[isprox].sum(), 100*contrib[isprox].sum()/contrib.sum(),
       contrib[isdist].sum(), 100*contrib[isdist].sum()/contrib.sum()))
for name, j in LEADS.items():
    if j is not None: print("  lead %-22s pos=%d z=%+.2f gene=%s" % (name, poso[j], zo[j], gene_of(poso[j])))

vcf_samples = set(subprocess.run([BCFTOOLS, "query", "-l", VCF], capture_output=True, text=True, check=True).stdout.split())
fin, eur = [], []
EUR_POPS = {"FIN", "CEU", "GBR", "TSI", "IBS"}
with open(PANEL) as fh:
    next(fh)
    for ln in fh:
        f = ln.split()
        if len(f) < 3: continue
        s, popc = f[0], f[1]
        if s not in vcf_samples: continue
        if popc == "FIN": fin.append(s)
        if popc in EUR_POPS: eur.append(s)
for nm, lst in [("fin", fin), ("eur", eur)]:
    open(f"{TMP}/{nm}.txt", "w").write("\n".join(lst) + "\n")
print("\nFIN n=%d ; EUR n=%d (both intersected with VCF's %d samples)" % (len(fin), len(eur), len(vcf_samples)))

with open(f"{TMP}/regions.txt", "w") as fh:
    for p in sorted(set(int(x) for x in S_pos)):
        fh.write(f"{VCF_CHROM}\t{p}\t{p}\n")

def build_R(sample_file, label):
    cmd = (f"{BCFTOOLS} view -S {sample_file} -R {TMP}/regions.txt -v snps -m2 -M2 {VCF} "
           f"| {BCFTOOLS} query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    n = sum(1 for _ in open(sample_file))
    G = np.full((len(S_idx), n), np.nan)
    matched = 0
    for ln in out.splitlines():
        f = ln.split("\t")
        p = int(f[0]); rr = f[1].upper(); aa = f[2].upper()
        if p not in posS: continue
        k = posS[p]; want_r, want_a = S_ref[k].upper(), S_alt[k].upper()
        if {rr, aa} != {want_r, want_a}: continue
        flip = (rr == want_a and aa == want_r)
        gts = f[3:]
        dos = np.empty(len(gts))
        for j, g in enumerate(gts):
            a0 = g[0]; a1 = g[2] if len(g) >= 3 else "."
            if a0 == "." or a1 == ".":
                dos[j] = np.nan
            else:
                c = (a0 == "1") + (a1 == "1")
                dos[j] = (2 - c) if flip else c
        m = np.isnan(dos)
        if m.any(): dos[m] = np.nanmean(dos)
        G[k, :] = dos; matched += 1
    sd = np.nanstd(G, axis=1)
    mono = np.where(sd == 0)[0]
    Gc = G.copy()
    R = np.eye(len(S_idx))
    good = np.where(sd > 0)[0]
    if len(good):
        Rg = np.corrcoef(Gc[good, :])
        for ii, gi in enumerate(good):
            for jj, gj in enumerate(good):
                R[gi, gj] = Rg[ii, jj]
    print("  [%s] matched %d/%d SNPs ; monomorphic in panel: %d" % (label, matched, len(S_idx), len(mono)))
    return R, matched, mono

boS = bo[S_idx]; zoS = zo[S_idx]
tag_local = [S_idx.index(i) for i in tag_idx]
lead_local = {name: (S_idx.index(j) if j is not None else None) for name, j in LEADS.items()}

def analyse(R, panel):
    R = np.nan_to_num(R, nan=0.0); np.fill_diagonal(R, 1.0)
    Doo = float(np.sqrt(max(boS @ (R @ boS), 1e-12)))
    zsv = float((boS @ zoS) / Doo)
    Rb = R @ boS
    print("\n=== panel: %s ===" % panel)
    print("  D (normalizer) = %.4f ; uncond z_SV = %.3f" % (Doo, zsv))
    def cond_single(jl):
        rj = R[:, jl]; den = np.sqrt(np.clip(1 - rj**2, 1e-6, None))
        zc = (zoS - rj * zoS[jl]) / den; zc[jl] = 0.0; return zc
    def cond_joint(idxs):
        C = np.array(idxs); Rcc = R[np.ix_(C, C)]; Ric = R[:, C]
        inv = np.linalg.pinv(Rcc); zc = zoS - Ric @ (inv @ zoS[C])
        var = np.clip(1 - np.einsum("ic,cd,id->i", Ric, inv, Ric), 1e-6, None)
        zc = zc / np.sqrt(var); zc[C] = 0.0; return zc
    res = {"D": Doo, "zsv": zsv}
    for name, jl in lead_local.items():
        if jl is None: continue
        zc = cond_single(jl)
        cz = float((boS @ zc) / Doo); r = float(Rb[jl] / Doo)
        print("  cond %-22s : z_SV %.3f -> %6.3f | r(sv40,lead)=%+.3f" % (name, zsv, cz, r))
        res[f"cond::{name}"] = cz; res[f"r::{name}"] = r
    alll = [jl for jl in lead_local.values() if jl is not None]
    zc = cond_joint(alll); cz = float((boS @ zc) / Doo)
    print("  cond JOINT(all leads)      : z_SV %.3f -> %6.3f" % (zsv, cz)); res["cond::JOINT"] = cz
    jp = lead_local["PROS1nbhd(93.8-94.4)"]
    if jp is not None:
        zc = cond_joint([jl for nm, jl in lead_local.items() if jl is not None and nm != "PROS1nbhd(93.8-94.4)"])
        cz = float((boS @ zc) / Doo)
        print("  cond all-EXCEPT-PROS1      : z_SV %.3f -> %6.3f" % (zsv, cz)); res["cond::exPROS1"] = cz
    return res

S_full = [int(obs[i]) for i in S_idx]
R_b_S = R_bundle[np.ix_(S_full, S_full)]
res_bundle = analyse(R_b_S, "EUR-525 (bundle, baseline)")
R_eur, m_eur, mono_eur = build_R(f"{TMP}/eur.txt", "EUR-1kGP")
res_eur = analyse(R_eur, "EUR-1kGP (n=%d)" % len(eur))
R_fin, m_fin, mono_fin = build_R(f"{TMP}/fin.txt", "FIN-1kGP")
res_fin = analyse(R_fin, "FIN-1kGP (n=%d)" % len(fin))

print("\n================ SUMMARY (key load-bearing numbers) ================")
hdr = ["metric", "EUR525bundle", "EUR1kGP(%d)" % len(eur), "FIN1kGP(%d)" % len(fin)]
print("%-26s %14s %14s %14s" % tuple(hdr))
keys = ["D", "zsv", "cond::PROS1nbhd(93.8-94.4)", "r::PROS1nbhd(93.8-94.4)",
        "cond::JOINT", "cond::exPROS1", "cond::interior(88.5-93.0)", "cond::prox(87.5-88.0)"]
for k in keys:
    row = [k] + ["%.3f" % d.get(k, float("nan")) for d in (res_bundle, res_eur, res_fin)]
    print("%-26s %14s %14s %14s" % tuple(row))
print("\nINTERNAL CONTROL: EUR-1kGP should ~match EUR525bundle (validates VCF->R pipeline).")
print("ANCESTRY TEST   : FIN vs EUR-1kGP isolates the Finnish-LD effect on the conditioning.")
print("DONE_FIN_LD_REDO")
