import pandas as pd, os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
PC=CENSV_ROOT + "/analyses/cdr_qtl/results/cdr_qtl_corrected/perchrom"
phenos=["mean_5mC","n_cdrs","cdrs_per_mb","cdr_fraction","median_cdr_width","cdr_pos_offset","cdr_pos_var","cdr_pos_iqr"]
hdr=("phenotype","n_sv","n_indiv","min_p","n_floor","top_hit")
print("%-18s %6s %7s %11s %7s  %s" % hdr)
floor_cells=[]
for ph in phenos:
    f=os.path.join(PC,"chr10_%s_perm.tsv"%ph)
    if not os.path.exists(f):
        print("%-18s  MISSING"%ph); continue
    d=pd.read_csv(f, sep="\t")
    nfloor=int((d.count_perm_ge_obs==0).sum())
    t=d.sort_values("p_emp_perm").iloc[0]
    desc="%d/%s/%d/af%.3f/t%+.2f/p%.1e/cnt%d" % (int(t.pos),t.svtype,int(t.svlen),t.af,t.t_obs,t.p_emp_perm,int(t.count_perm_ge_obs))
    print("%-18s %6d %7d %11.2e %7d  %s" % (ph,len(d),int(t.n_indiv),t.p_emp_perm,nfloor,desc))
    if nfloor>0: floor_cells.append(ph)
print()
print("CHR10 FLOOR CELLS (need 1e6 escalation):", floor_cells if floor_cells else "NONE")
