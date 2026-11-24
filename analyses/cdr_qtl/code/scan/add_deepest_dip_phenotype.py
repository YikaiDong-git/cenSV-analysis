#!/usr/bin/env python
"""Add a DEEPEST-DIP position phenotype alongside the existing median-of-dips cdr_pos_offset."""
import os, re, glob, numpy as np, pandas as pd
from collections import defaultdict
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
D=CENSV_ROOT + "/analyses/cdr_calls"
OUT_DIR=CENSV_ROOT + "/analyses/cdr_qtl/output"
GLOBE_AA=f"{D}/inputs/globe_aa"; CORR=f"{D}/output/calls"; COHORT=f"{D}/output/cohort_corrected_final.tsv"
ROW_RE=re.compile(r"^(.+?)_(chr[0-9XY]+)_((.+?)#([12])#[^:]+):")

def build_globe_map():
    m=defaultdict(list)
    for bed in sorted(glob.glob(f"{GLOBE_AA}/*_AS-HOR_stv_row.ort.bed")):
        for line in open(bed):
            p=line.rstrip("\n").split("\t")
            if len(p)<3: continue
            mm=ROW_RE.match(p[0])
            if not mm: continue
            try: s,e=int(p[1]),int(p[2])
            except ValueError: continue
            m[(mm.group(1),f"hap{mm.group(5)}",mm.group(2))].append((mm.group(3),s,e))
    return m

def load_bed(path,ncol=3):
    out=defaultdict(list)
    if not os.path.exists(path): return out
    for line in open(path):
        p=line.rstrip("\n").split("\t")
        if len(p)<ncol: continue
        try: s,e=int(p[1]),int(p[2])
        except ValueError: continue
        out[p[0]].append((s,e,p[3:]))
    return out

gmap=build_globe_map()
cdrs_by={}; bins_by={}
for tagdir in sorted(glob.glob(f"{CORR}/*/")):
    for ci,lst in load_bed(f"{tagdir}/cdrs.bed").items(): cdrs_by.setdefault(ci,[]).extend(lst)
    for ci,lst in load_bed(f"{tagdir}/binned_freq.bed",4).items(): bins_by.setdefault(ci,[]).extend(lst)
print(f"[globe] {len(gmap)} groups; {len(cdrs_by)} contigs w/ CDRs, {len(bins_by)} w/ windows")

def cdr_min5mc(ci, cs, ce):
    vals=[]
    for (ws,we,extra) in bins_by.get(ci,[]):
        if we>cs and ws<ce and extra:
            try: vals.append(float(extra[0]))
            except ValueError: pass
    return min(vals) if vals else np.nan

rows=[]
for (sample,asm_hap,chrom),intervals in sorted(gmap.items()):
    if not any(ci in bins_by for (ci,s,e) in intervals): continue
    aa_s=min(s for (_,s,e) in intervals); aa_e=max(e for (_,s,e) in intervals)
    aa_len=aa_e-aa_s
    if aa_len<=0: continue
    aa_mid=(aa_s+aa_e)/2.0; aa_half=aa_len/2.0
    cdr_list=[]
    for (ci,s,e) in intervals:
        for (cs,ce,_x) in cdrs_by.get(ci,[]):
            if ce>s and cs<e:
                a,b=max(cs,s),min(ce,e); mid=(a+b)/2.0
                cdr_list.append((mid, cdr_min5mc(ci,a,b)))
    n=len(cdr_list)
    if n==0:
        continue
    mids=np.array([m for m,_ in cdr_list]); rel=(mids-aa_mid)/aa_half
    pos_med=float(np.median(rel))
    depths=np.array([d for _,d in cdr_list], dtype=float)
    if np.all(np.isnan(depths)):
        pos_deep=pos_med; dmin=np.nan; dmid=np.nan
    else:
        j=int(np.nanargmin(depths)); pos_deep=float(rel[j]); dmin=float(depths[j]); dmid=float(mids[j])
    wp=[]
    for (ci,s,e) in intervals:
        for (ws,we,extra) in bins_by.get(ci,[]):
            if we>s and ws<e and extra:
                try: wp.append(float(extra[0]))
                except ValueError: pass
    rows.append(dict(sample=sample,hap_local=asm_hap,chm13_chrom=chrom,
                     active_start=aa_s,active_end=aa_e,active_len=aa_len,n_cdrs=n,
                     cdr_pos_offset=pos_med, cdr_pos_offset_deepest=pos_deep,
                     deepest_min5mC=dmin, deepest_cdr_mid=dmid,
                     mean_5mC=float(np.mean(wp)) if wp else np.nan))
df=pd.DataFrame(rows)
out=f"{OUT_DIR}/per_hap_chrom_phenotypes_deepest.tsv"
df.to_csv(out,sep="\t",index=False,float_format="%.6g")
print(f"[written] {out}: {len(df)} rows, {df['sample'].nunique()} samples")
print(f"\nmedian-vs-deepest offset corr (all): r={df[['cdr_pos_offset','cdr_pos_offset_deepest']].corr().iloc[0,1]:.3f}")
c8=df[df.chm13_chrom=='chr8']
print(f"chr8: n={len(c8)}  median-offset mean={c8.cdr_pos_offset.mean():.3f}  deepest-offset mean={c8.cdr_pos_offset_deepest.mean():.3f}  corr={c8[['cdr_pos_offset','cdr_pos_offset_deepest']].corr().iloc[0,1]:.3f}")
print(f"chr8 deepest_min5mC: median={c8.deepest_min5mC.median():.1f}%  (should be LOW = deep valley)")
