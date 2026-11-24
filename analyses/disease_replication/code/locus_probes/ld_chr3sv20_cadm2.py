import numpy as np, subprocess, os
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]
P    = CENSV_WORKSPACE + ""
SSIMP= f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
BUND = f"{SSIMP}/output/bundles/chr3/chr3_super1_85533081_94536472_bundle.npz"
BT   = f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"
BCF  = f"{BT}/bcftools"
RES  = f"{P}/Project1_Centromere/resources/1000G"
VCF  = f"{RES}/1kGP_high_coverage_Illumina.chr3.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PANEL= f"{P}/_shared/datasets/1KGP/v_canonical_phase3_20130502/EUR/integrated_call_samples_v3.20130502.ALL.panel"
EUR_POPS = {"FIN","CEU","GBR","TSI","IBS"}
TMP = "/tmp/chr3_cadm2_ld"; os.makedirs(TMP, exist_ok=True)
CADM2 = (85008349, 86101271)

d = np.load(BUND, allow_pickle=True)
svn = np.array([str(x) for x in d["sv_locus"]]); ti = int(np.where(svn=="chr3_sv20_DEL_non_best6_alt_edge1000")[0][0])
w = d["B"][:,ti].astype(float); sup = np.where(w!=0)[0]
spos = d["snp_pos"].astype(int)[sup]; sref = np.array([str(x).upper() for x in d["snp_ref"]])[sup]
salt = np.array([str(x).upper() for x in d["snp_alt"]])[sup]; sw = w[sup]
proximal = spos < 90900000
print(f"chr3_sv20: {len(sup)} weight SNPs; proximal(<cen)={proximal.sum()} distal={(~proximal).sum()}")

vs = set(subprocess.run([BCF,"query","-l",VCF],capture_output=True,text=True,check=True).stdout.split())
eur = []
with open(PANEL) as fh:
    next(fh)
    for ln in fh:
        f = ln.split()
        if len(f)>=3 and f[0] in vs and f[1] in EUR_POPS: eur.append(f[0])
open(f"{TMP}/eur.txt","w").write("\n".join(eur)+"\n")
print(f"EUR n={len(eur)}")

def read_geno(arg, region):
    cmd = (f"{BCF} view -S {TMP}/eur.txt {'-r' if region else '-R'} {arg} -v snps -m2 -M2 {VCF} "
           f"| {BCF} query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    g = {}
    for ln in out.splitlines():
        f = ln.split("\t"); p=int(f[0]); r=f[1].upper(); a=f[2].upper(); gts=f[3:]
        dos = np.array([ (gg[0]=="1")+(gg[2]=="1") if (len(gg)>=3 and gg[0]!="." and gg[2]!=".") else np.nan for gg in gts], float)
        m = np.isnan(dos)
        if m.all(): continue
        dos[m] = np.nanmean(dos); g[p] = (r,a,dos)
    return g

with open(f"{TMP}/pred_pos.txt","w") as fh:
    for p in sorted(set(int(x) for x in spos)): fh.write(f"chr3\t{p}\t{p}\n")
gp = read_geno(f"{TMP}/pred_pos.txt", region=False)

def build_pred(mask):
    v = np.zeros(len(eur)); n=0
    for k in range(len(spos)):
        if not mask[k]: continue
        p=int(spos[k]); rr=sref[k]; aa=salt[k]
        if p in gp:
            vr,va,dos = gp[p]
            if {vr,va}=={rr,aa}:
                g_alt = (2-dos) if (vr==aa and va==rr) else dos
                v += sw[k]*g_alt; n+=1
    return v-v.mean(), n

g_full,nf = build_pred(np.ones(len(spos),bool))
g_prox,npx= build_pred(proximal)
print(f"predictor built: full={nf}/{len(spos)} snps, proximal={npx}")

cadm2 = read_geno(f"chr3:{CADM2[0]}-{CADM2[1]}", region=True)
print(f"CADM2-window ({CADM2[0]}-{CADM2[1]}) common SNPs read: {len(cadm2)}")
def topld(gv, label):
    res=[]
    for p,(r,a,dos) in cadm2.items():
        dd=dos-dos.mean()
        if dd.std()>0 and gv.std()>0:
            rr=float(np.corrcoef(gv,dd)[0,1]); res.append((rr*rr,p,r,a,rr))
    res.sort(reverse=True)
    print(f"\n[{label}] top LD (r2) vs CADM2-window SNPs:")
    for r2,p,r,a,rr in res[:8]:
        print(f"   chr3:{p} {r}/{a}  r2={r2:.3f} r={rr:+.2f}  ({'inside' if CADM2[0]<=p<=CADM2[1] else 'flank'})")
    print(f"   >>> MAX r2 = {res[0][0]:.3f} at chr3:{res[0][1]}" if res else "   none")
topld(g_full, "FULL predictor")
topld(g_prox, "PROXIMAL-only (87.5Mb cluster)")
print("DONE_LD")
