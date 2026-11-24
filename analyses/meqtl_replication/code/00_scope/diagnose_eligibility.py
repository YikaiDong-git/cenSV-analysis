#!/usr/bin/env python3
import os
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
AUD = PROJ + "/analyses/cdr_calls"
BED = PROJ + "/subproject_methylation/raw/bedmethyl_filt"
FAI = PROJ + "/subproject_methylation/raw/asm_fai"
BAM = PROJ + "/alignments/svim_asm"

H = pd.read_csv(SUB + "/output/mqtl/bam_hap_map.tsv", sep="\t")
H["resolved"] = H.bam_tag.notna() & (H.bam_tag.astype(str).str.len() > 0) \
    & (H.bam_tag.astype(str) != "nan")
print("bam_hap_map rows: " + str(len(H)))
print(H.groupby(["resolved", "agrees_with_manifest"], dropna=False).size().to_string())

R = H[H.resolved].copy()
R["has_bed"] = [os.path.exists(BED + "/" + t + ".active.bed") for t in R.meth_tag]
R["has_fai"] = [os.path.exists(FAI + "/" + t.rsplit("_", 1)[0] + "_" + str(fl) + ".fa.gz.fai")
                for t, fl in zip(R.meth_tag, R.fai_label)]
R["has_bam"] = [os.path.exists(BAM + "/" + str(b) + "_to_chm13.bam") for b in R.bam_tag]

chem = pd.read_csv(AUD + "/output/cohort_corrected_final.tsv", sep="\t")
chem = chem[chem.chemistry.isin(["R941", "R1041"])]
cs = chem.groupby("sample").chemistry.apply(
    lambda x: "R1041" if "R1041" in set(x) else "R941")
R["sample"] = R.meth_tag.str.rsplit("_", n=1).str[0]
R["has_chem"] = R["sample"].map(cs).notna()

print("\nfile availability by label agreement:")
print(R.groupby("agrees_with_manifest")[["has_bed", "has_fai", "has_bam", "has_chem"]]
      .agg(["sum", "size"]).to_string())

R["eligible"] = R.has_bed & R.has_fai & R.has_bam & R.has_chem
print("\neligible by agreement:")
print(R.groupby("agrees_with_manifest").eligible.agg(["sum", "size"]).to_string())

inv = R[~R.agrees_with_manifest.astype(bool)]
print("\n--- the 183 inverted rows, first blocking gate ---")
for col in ["has_bed", "has_fai", "has_bam", "has_chem"]:
    print("  " + col + " missing for " + str(int((~inv[col]).sum())) + " / " + str(len(inv)))
print("\nexample inverted rows:")
print(inv[["meth_tag", "bam_tag", "has_bed", "has_fai", "has_bam", "has_chem"]]
      .head(6).to_string(index=False))

if len(inv):
    s = inv["sample"].iloc[0]
    print("\nwhat exists for sample " + s + ":")
    for d, lab in ((BED, "bedmethyl_filt"), (FAI, "asm_fai"), (BAM, "svim_asm")):
        got = sorted(f for f in os.listdir(d) if f.startswith(s))
        print("  " + lab + ": " + (", ".join(got[:6]) if got else "(none)"))

elig = R[R.eligible].copy()
print("\n=== COHORT SIZE ===")
print("eligible haplotypes : " + str(len(elig)))
print("distinct individuals: " + str(elig["sample"].nunique()))
both = elig.groupby("sample").size()
print("individuals with BOTH haplotypes: " + str(int((both == 2).sum())))
print("individuals with ONE haplotype  : " + str(int((both == 1).sum())))
elig["chemistry"] = elig["sample"].map(cs)
print(elig.groupby("chemistry")["sample"].nunique().to_string())
