#!/usr/bin/env python
"""Resolve, per methylation haplotype, WHICH svim-asm BAM is actually the same haplotype."""
import os, gzip, glob
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SUB = f"{P1}/subproject_methylation"
WU = f"{P1}/analyses/cdr_qtl_full_panel"
AUD = f"{P1}/analyses/cdr_calls"
FAI = f"{SUB}/raw/asm_fai"
OUT = f"{WU}/output/mqtl/bam_hap_map.tsv"

MAN = pd.read_csv(f"{SUB}/manifest/methylation_manifest.tsv", sep="\t")
REMOTE = {(r["sample"], r.hap_local): str(r.hap_remote) for _, r in MAN.iterrows()}

tags = sorted(os.path.basename(p) for p in glob.glob(f"{AUD}/output/calls/*")
              if os.path.isdir(p))
print(f"methylation haplotypes: {len(tags)}")


def fai_names(p):
    return {l.split("\t")[0] for l in open(p)}


def fai_lens(p):
    return [int(l.split("\t")[1]) for l in open(p)]


def paf_lens(p):
    d = {}
    with gzip.open(p, "rt") as fh:
        for line in fh:
            f = line.split("\t", 3)
            if len(f) >= 2:
                d[f[0]] = int(f[1])
    return d


rows = []
for tag in tags:
    sample, hap = tag.rsplit("_", 1)
    bed = f"{AUD}/output/calls/{tag}/binned_freq.bed"
    if not os.path.exists(bed):
        continue
    bed_ctgs = {l.split("\t")[0] for l in open(bed)}

    fai_label, best = None, 0
    for lab in ("hap1", "hap2", "mat", "pat"):
        p = f"{FAI}/{sample}_{lab}.fa.gz.fai"
        if not os.path.exists(p):
            continue
        n = len(bed_ctgs & fai_names(p))
        if n > best:
            fai_label, best = lab, n
    if fai_label is None or best < 0.5 * len(bed_ctgs):
        rows.append(dict(meth_tag=tag, fai_label=fai_label or "", bam_tag="",
                         n_shared_lengths=0, agrees_with_manifest=""))
        continue

    want = set(fai_lens(f"{FAI}/{sample}_{fai_label}.fa.gz.fai"))
    bam_tag, bbest = "", 0
    for h in ("hap1", "hap2"):
        p = f"{P1}/alignments/dipcall_paf/{sample}_{h}_to_chm13.bam"
        if not os.path.exists(p):
            continue
        n = len(want & set(paf_lens(p).values()))
        if n > bbest:
            bam_tag, bbest = f"{sample}_{h}", n
    rows.append(dict(meth_tag=tag, fai_label=fai_label, bam_tag=bam_tag,
                     n_shared_lengths=bbest,
                     agrees_with_manifest=(bam_tag.rsplit("_", 1)[-1] == hap) if bam_tag else ""))

D = pd.DataFrame(rows)
D.to_csv(OUT, sep="\t", index=False)
res = D[D.bam_tag != ""]
print(f"resolved {len(res)}/{len(D)} haplotypes")
print(f"  BAM label AGREES with the methylation label : {int((res.agrees_with_manifest == True).sum())}")
print(f"  BAM label INVERTED vs the methylation label : {int((res.agrees_with_manifest == False).sum())}")
print(f"  unresolved (no .fai or no PAF)              : {int((D.bam_tag == '').sum())}")
print(f"\nfai label used, by count:\n{D.fai_label.value_counts().to_string()}")
inv = res[res.agrees_with_manifest == False]
if len(inv):
    print(f"\nfirst inverted examples:\n{inv.head(8).to_string(index=False)}")
print(f"\nwrote {OUT}")
