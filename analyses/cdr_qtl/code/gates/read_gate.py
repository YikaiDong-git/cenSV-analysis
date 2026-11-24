import json, os, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
WU=CENSV_ROOT + "/analyses/cdr_qtl"
r=json.load(open(f"{WU}/results/figures/05_confound/chr10_gate.json"))
for lab,v in r.items():
    print("\n"+lab+":  OVERALL="+str(v["overall"]))
    print("  geom_ok=%s  (len~dos p=%.2f, width~dos p=%.2f, pheno~dos t=%+.1f -> +len t=%+.1f, atten=%.0f%%)" % (
        v["geom_ok"], v["T1_len_p"], v["T2_width_p"], v["T3_t"], v["T4_t"], v["atten_pct"]))
    print("  abs_ok=%s  (ABS t=%s p=%s)" % (v["abs_ok"], v["ABS_t"], v["ABS_p"]))
    print("  ancestry_gate=%s  dropped=%s" % (v["gate"], v["dropped"]))
    print("  meta=%s" % (v["meta"],))
    print("  LOSO=%s" % (v["loso"],))
print("\n--- chr10 in sv_density figure data tables? ---")
for t in ["04_sv_density/cross_chrom_table.tsv", "04_sv_density/cdr_position_by_superpop.tsv"]:
    p=f"{WU}/results/figures/{t}"
    if os.path.exists(p):
        d=pd.read_csv(p, sep="\t")
        ccol=next((c for c in d.columns if "chrom" in c.lower()), d.columns[0])
        has=d[ccol].astype(str).str.contains("chr10").any()
        print(f"  {os.path.basename(t)}: rows={len(d)}, chr10_present={has}")
