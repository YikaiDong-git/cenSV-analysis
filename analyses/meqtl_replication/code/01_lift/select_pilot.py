#!/usr/bin/env python3
"""Choose the 8 pilot haplotypes."""
import os
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
AUD = PROJ + "/analyses/cdr_calls"
BED = PROJ + "/subproject_methylation/raw/bedmethyl_filt"
FAI = PROJ + "/subproject_methylation/raw/asm_fai"
BAM = PROJ + "/alignments/svim_asm"

H = pd.read_csv(SUB + "/output/mqtl/bam_hap_map.tsv", sep="\t")
H = H[H.bam_tag.notna() & (H.bam_tag.astype(str) != "")]

chem = pd.read_csv(AUD + "/output/cohort_corrected_final.tsv", sep="\t")
chem = chem[chem.chemistry.isin(["R941", "R1041"])]
cs = chem.groupby("sample").chemistry.apply(
    lambda x: "R1041" if "R1041" in set(x) else "R941")

H["sample"] = H.meth_tag.str.rsplit("_", n=1).str[0]
H["chemistry"] = H["sample"].map(cs)
H["files_ok"] = [
    os.path.exists(BED + "/" + t + ".active.bed")
    and os.path.exists(FAI + "/" + t.rsplit("_", 1)[0] + "_" + str(fl) + ".fa.gz.fai")
    and os.path.exists(BAM + "/" + str(b) + "_to_chm13.bam")
    for t, b, fl in zip(H.meth_tag, H.bam_tag, H.fai_label)]

E = H[H.files_ok & H.chemistry.notna()].copy()
print("eligible haplotypes: " + str(len(E)))
print(E.groupby(["chemistry", "agrees_with_manifest"]).size().to_string())

picks = []
for c in ["R941", "R1041"]:
    for a in [True, False]:
        cell = E[(E.chemistry == c) & (E.agrees_with_manifest == a)].sort_values("meth_tag")
        got = cell.head(2)
        if len(got) < 2:
            print("  [cell thin] chemistry=" + c + " agrees=" + str(a)
                  + " has only " + str(len(cell)))
        picks.append(got)
P = pd.concat(picks).reset_index(drop=True)
print("\npilot selection:")
print(P[["meth_tag", "bam_tag", "chemistry", "agrees_with_manifest"]].to_string(index=False))

out = UNIT + "/output/a1_lift/pilot_tags.txt"
tmp = out + ".tmp." + str(os.getpid())
with open(tmp, "w") as f:
    for t in P.meth_tag:
        f.write(t + "\n")
os.replace(tmp, out)
print("\nwrote " + out + "  (" + str(len(P)) + " tags)")

allout = UNIT + "/output/a1_lift/all_tags.txt"
tmp = allout + ".tmp." + str(os.getpid())
with open(tmp, "w") as f:
    for t in sorted(E.meth_tag):
        f.write(t + "\n")
os.replace(tmp, allout)
print("wrote " + allout + "  (" + str(len(E)) + " tags, for the full array)")
