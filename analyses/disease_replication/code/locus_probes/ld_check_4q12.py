#!/usr/bin/env python
import numpy as np, subprocess, os
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]
P=CENSV_WORKSPACE + ""; SS=f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
BUND=f"{SS}/output/bundles/chr4/chr4_super1_47409879_53011266_bundle.npz"
BT=f"{P}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin"; BCF=f"{BT}/bcftools"
VCF=f"{P}/Project1_Centromere/resources/1000G/1kGP_high_coverage_Illumina.chr4.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PANEL=f"{P}/_shared/datasets/1KGP/v_canonical_phase3_20130502/EUR/integrated_call_samples_v3.20130502.ALL.panel"
EUR={"FIN","CEU","GBR","TSI","IBS"}; TMP="/tmp/ld4q12"; os.makedirs(TMP,exist_ok=True)
vs=set(subprocess.run([BCF,"query","-l",VCF],capture_output=True,text=True,check=True).stdout.split())
eur=[f.split()[0] for f in open(PANEL).read().splitlines()[1:] if len(f.split())>=3 and f.split()[0] in vs and f.split()[1] in EUR]
open(f"{TMP}/eur.txt","w").write("\n".join(eur)+"\n")
def geno(positions,tag):
    tf=f"{TMP}/p_{tag}.txt"; open(tf,"w").write("".join(f"chr4\t{p}\t{p}\n" for p in sorted(set(positions))))
    out=subprocess.run(f"{BCF} view -S {TMP}/eur.txt -R {tf} -v snps -m2 -M2 {VCF} | {BCF} query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'",shell=True,capture_output=True,text=True,check=True).stdout
    g={}
    for ln in out.splitlines():
        f=ln.split("\t");p=int(f[0]);gts=f[3:]
        dos=np.array([(x[0]=="1")+(x[2]=="1") if(len(x)>=3 and x[0]!="."and x[2]!=".")else np.nan for x in gts],float)
        m=np.isnan(dos)
        if m.all():continue
        dos[m]=np.nanmean(dos);g[p]=(f[1].upper(),f[2].upper(),dos)
    return g
gg=geno([52023436,48180521],"leads")
a=gg[52023436][2]-gg[52023436][2].mean(); b=gg[48180521][2]-gg[48180521][2].mean()
r=float(np.corrcoef(a,b)[0,1]); print(f"LEAD-SNP LD  chr4:52023436(height) vs chr4:48180521(PDW): r={r:+.3f} r2={r*r:.3f}")
d=np.load(BUND,allow_pickle=True); svn=np.array([str(x) for x in d["sv_locus"]])
def pred(svid):
    ti=int(np.where(svn==svid)[0][0]); w=d["B"][:,ti].astype(float); sup=np.where(w!=0)[0]
    pos=d["snp_pos"].astype(int)[sup]; ref=np.array([str(x).upper() for x in d["snp_ref"]])[sup]; alt=np.array([str(x).upper() for x in d["snp_alt"]])[sup]
    gp=geno(list(pos),svid.split("_")[1]); v=np.zeros(len(eur))
    for k in range(len(pos)):
        pp=int(pos[k])
        if pp in gp:
            vr,va,dos=gp[pp]
            if {vr,va}=={ref[k],alt[k]}: v+=w[sup][k]*((2-dos) if(vr==alt[k] and va==ref[k]) else dos)
    return v-v.mean()
p1=pred("chr4_sv349_DEL_non_best8_fksize200k"); p2=pred("chr4_sv398_INS_non_best8_fksize100k")
rr=float(np.corrcoef(p1,p2)[0,1]); print(f"PREDICTOR LD  sv349(height) vs sv398(PDW): r={rr:+.3f} r2={rr*rr:.3f}")
print("VERDICT:", "DISTINCT sub-loci (r2<0.1)" if rr*rr<0.1 and r*r<0.1 else "LD-LINKED (not distinct)")
print("DONE_LD4Q12")
