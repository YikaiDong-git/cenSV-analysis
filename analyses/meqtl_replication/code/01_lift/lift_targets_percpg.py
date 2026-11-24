#!/usr/bin/env python3
import os
import sys
import gzip
import numpy as np
import pandas as pd
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
AUD = PROJ + "/analyses/cdr_calls"
BAM_DIR = PROJ + "/alignments/svim_asm"
BED_DIR = PROJ + "/subproject_methylation/raw/bedmethyl_filt"
BEDM_DIR = PROJ + "/subproject_methylation/raw/bedmethyl_filt_m"
FAI_DIR = PROJ + "/subproject_methylation/raw/asm_fai"
OUT_DIR = UNIT + "/output/a1_lift/per_hap"

MIN_MAPQ = 20
CPG_OFFSETS = (0, -1, 1)
CTRL_OFFSETS = (5003, -5003)

CONSUME_REF = (0, 2, 3, 7, 8)
CONSUME_QRY_STORED = (0, 1, 4, 7, 8)

TAG = sys.argv[1]
os.makedirs(OUT_DIR, exist_ok=True)


def ref_to_query_stored(aln, p):
    rpos = aln.reference_start
    qpos = 0
    for op, ln in aln.cigartuples:
        if op in (0, 7, 8):
            if rpos <= p < rpos + ln:
                return qpos + (p - rpos), "LIFTED"
            rpos += ln
            qpos += ln
        elif op in (2, 3):
            if rpos <= p < rpos + ln:
                return None, "IN_DELETION"
            rpos += ln
        elif op in (1, 4):
            qpos += ln
        elif op in (5, 6):
            pass
        else:
            raise ValueError("unexpected CIGAR op " + str(op) + " in " + TAG)
    return None, "WALK_OVERRAN"


HMAP = pd.read_csv(SUB + "/output/mqtl/bam_hap_map.tsv", sep="\t")
row = HMAP[HMAP.meth_tag == TAG]
if not len(row):
    print("[UNRESOLVED] " + TAG + ": absent from bam_hap_map.tsv -- data gap, not error")
    sys.exit(0)
bam_tag = row.bam_tag.iloc[0]
if not isinstance(bam_tag, str) or not bam_tag:
    print("[UNRESOLVED] " + TAG + ": haplotype identity unresolved (no published .fai)")
    sys.exit(0)
print(TAG + ": methylation tag -> BAM tag " + bam_tag
      + ("  (label AGREES)" if bam_tag == TAG else "  (label INVERTED -- T1)"))

sample = TAG.rsplit("_", 1)[0]
fai_label = row.fai_label.iloc[0]

bam_path = BAM_DIR + "/" + bam_tag + "_to_chm13.bam"
bed_path = BED_DIR + "/" + TAG + ".active.bed"
fai_path = FAI_DIR + "/" + sample + "_" + str(fai_label) + ".fa.gz.fai"
for pth in (bam_path, bed_path, fai_path):
    if not os.path.exists(pth):
        raise FileNotFoundError(pth)

fai = pd.read_csv(fai_path, sep="\t", header=None,
                  names=["contig", "length", "offset", "linebases", "linewidth"])
len_counts = fai.length.value_counts()
len2contig = {int(l): c for l, c in zip(fai.length, fai.contig)}
ambiguous_lengths = set(int(l) for l in len_counts[len_counts > 1].index)
print("  fai " + os.path.basename(fai_path) + ": " + str(len(fai))
      + " contigs; lengths that are not unique: " + str(len(ambiguous_lengths)))

FAI_PANSN = set(c.rsplit("#", 1)[0] + "#" for c in fai.contig if "#" in c)
if len(FAI_PANSN) != 1:
    raise ValueError(TAG + ": .fai carries " + str(len(FAI_PANSN))
                     + " PanSN prefixes, expected exactly 1: " + str(sorted(FAI_PANSN)))
FAI_PREFIX = sorted(FAI_PANSN)[0]
print("  PanSN prefix from .fai: " + FAI_PREFIX)

tg = pd.read_csv(UNIT + "/output/targets/targets_cpg.tsv", sep="\t")
print("  targets " + str(len(tg)))

chem = pd.read_csv(AUD + "/output/cohort_corrected_final.tsv", sep="\t")
chem = chem[chem.chemistry.isin(["R941", "R1041"])]
sc = chem[chem["sample"] == sample].chemistry.unique()
chemistry = "R1041" if "R1041" in sc else ("R941" if "R941" in sc else "NA")

bam = pysam.AlignmentFile(bam_path, "rb")
refs = set(bam.references)

idx_total = sum(s.mapped + s.unmapped for s in bam.get_index_statistics())
print("  index check: " + str(idx_total) + " records across "
      + str(bam.nreferences) + " references (idxstats); "
      + "run code/01_lift/verify_bam_index.py for the full-stream comparison")
if idx_total == 0:
    raise ValueError(TAG + ": BAM index reports zero records -- index is stale or broken")

lifted = []
prefix_mismatch = 0
for t in tg.itertuples(index=False):
    rec = dict(meth_tag=TAG, sample=sample, chemistry=chemistry, bam_tag=bam_tag,
               target_id=t.target_id, chm13_chrom=t.chm13_chrom, chm13_pos=int(t.chm13_pos),
               lift_status="NO_ALIGNMENT", n_aln_covering=0, bam_query="", bed_contig="",
               contig_pos=-1, is_reverse=False, mapq=-1, is_supplementary=False,
               contig_len=-1)
    if t.chm13_chrom not in refs:
        rec["lift_status"] = "CHROM_ABSENT_FROM_BAM"
        lifted.append(rec)
        continue
    p = int(t.chm13_pos)
    cands = []
    for aln in bam.fetch(t.chm13_chrom, p, p + 1):
        if aln.is_unmapped or aln.cigartuples is None:
            continue
        if aln.mapping_quality < MIN_MAPQ:
            continue
        cands.append(aln)
    rec["n_aln_covering"] = len(cands)
    if not cands:
        lifted.append(rec)
        continue
    cands.sort(key=lambda a: (a.is_supplementary or a.is_secondary,
                              -a.mapping_quality,
                              -(a.reference_end - a.reference_start),
                              a.query_name,
                              a.reference_start))
    aln = cands[0]
    if "#" in aln.query_name and not aln.query_name.startswith(FAI_PREFIX):
        prefix_mismatch += 1
        rec["lift_status"] = "PANSN_PREFIX_MISMATCH"
        rec["bam_query"] = aln.query_name
        lifted.append(rec)
        continue
    qpos_stored, status = ref_to_query_stored(aln, p)
    rec["lift_status"] = status
    rec["bam_query"] = aln.query_name
    rec["is_reverse"] = bool(aln.is_reverse)
    rec["mapq"] = int(aln.mapping_quality)
    rec["is_supplementary"] = bool(aln.is_supplementary or aln.is_secondary)
    if status != "LIFTED":
        lifted.append(rec)
        continue
    lead_h = aln.cigartuples[0][1] if aln.cigartuples[0][0] == 5 else 0
    L_full = aln.infer_read_length()
    qpos_full = qpos_stored + lead_h
    assert 0 <= qpos_full < L_full, (TAG, t.target_id, qpos_full, L_full)
    rec["contig_len"] = int(L_full)
    rec["contig_pos"] = int(L_full - 1 - qpos_full) if aln.is_reverse else int(qpos_full)
    if L_full in ambiguous_lengths:
        rec["lift_status"] = "AMBIGUOUS_CONTIG_LENGTH"
    elif L_full not in len2contig:
        rec["lift_status"] = "NO_CONTIG_MAP"
    else:
        rec["bed_contig"] = len2contig[int(L_full)]
    lifted.append(rec)
bam.close()

n_with_aln = sum(1 for r in lifted if r["n_aln_covering"] > 0)
if n_with_aln and prefix_mismatch / float(n_with_aln) > 0.05:
    raise ValueError(TAG + ": " + str(prefix_mismatch) + "/" + str(n_with_aln)
                     + " covered targets have a BAM query outside PanSN prefix "
                     + FAI_PREFIX + " -- haplotype pairing is wrong (see T1/T1c)")

L = pd.DataFrame(lifted)
n_ok = int((L.lift_status == "LIFTED").sum())
print("  PASS1: " + str(n_ok) + "/" + str(len(L)) + " lifted; status breakdown:")
print(L.lift_status.value_counts().to_string())
print("  reverse-strand share of lifted: "
      + str(round(float(L[L.lift_status == "LIFTED"].is_reverse.mean()) * 100, 2)) + "%")

want = {}
for r in L[L.lift_status == "LIFTED"].itertuples(index=False):
    if not r.bed_contig:
        continue
    for off in CPG_OFFSETS:
        want.setdefault((r.bed_contig, r.contig_pos + off), []).append((r.target_id, "cpg", off))
    for off in CTRL_OFFSETS:
        want.setdefault((r.bed_contig, r.contig_pos + off), []).append((r.target_id, "ctrl", off))
want_contigs = set(c for c, _ in want)
print("  PASS2: " + str(len(want)) + " positions on " + str(len(want_contigs)) + " contigs")

found = {}
n_rows = 0
n_m = 0
opener = gzip.open if bed_path.endswith(".gz") else open
with opener(bed_path, "rt") as f:
    for line in f:
        i = line.find("\t")
        contig = line[:i]
        if contig not in want_contigs:
            continue
        n_rows += 1
        parts = line.rstrip("\n").split("\t")
        if parts[3] != "m":
            continue
        n_m += 1
        key = (contig, int(parts[1]))
        if key in want:
            found[key] = (int(parts[4]), float(parts[10]))
print("  scanned " + str(n_rows) + " rows on target contigs, " + str(n_m)
      + " were mod-type m; matched " + str(len(found)) + " positions")

cpg_off, cpg_cov, cpg_pct = [], [], []
ctrl_hit, ctrl_cov, ctrl_pct = [], [], []
for r in L.itertuples(index=False):
    o, cv, pc = np.nan, np.nan, np.nan
    if r.lift_status == "LIFTED" and r.bed_contig:
        for off in CPG_OFFSETS:
            k = (r.bed_contig, r.contig_pos + off)
            if k in found:
                o, (cv, pc) = off, found[k]
                break
    cpg_off.append(o)
    cpg_cov.append(cv)
    cpg_pct.append(pc)
    hits, cvs, pcs = 0, [], []
    if r.lift_status == "LIFTED" and r.bed_contig:
        for off in CTRL_OFFSETS:
            k = (r.bed_contig, r.contig_pos + off)
            if k in found:
                hits += 1
                cvs.append(found[k][0])
                pcs.append(found[k][1])
    ctrl_hit.append(hits)
    ctrl_cov.append(float(np.mean(cvs)) if cvs else np.nan)
    ctrl_pct.append(float(np.mean(pcs)) if pcs else np.nan)

L["cpg_hit_offset"] = cpg_off
L["cov"] = cpg_cov
L["pct_m"] = cpg_pct
L["n_ctrl_hit"] = ctrl_hit
L["ctrl_cov"] = ctrl_cov
L["ctrl_pct_m"] = ctrl_pct
L["has_filt_m_file"] = os.path.exists(BEDM_DIR + "/" + TAG + ".active.m.bed")

ok = L[L.lift_status == "LIFTED"]
hit = ok.cpg_hit_offset.notna()
ctrl_rate = float(ok.n_ctrl_hit.sum()) / max(1, len(ok) * len(CTRL_OFFSETS))
print("\n  === VALIDATION (this is the check that can fail) ===")
print("  CpG hit rate      : " + str(round(float(hit.mean()) * 100, 2)) + "%  ("
      + str(int(hit.sum())) + "/" + str(len(ok)) + ")")
print("  control hit rate  : " + str(round(ctrl_rate * 100, 2)) + "%   <- background, expect ~1.6%")
for rev in (False, True):
    m = ok.is_reverse == rev
    if int(m.sum()) == 0:
        continue
    sub = ok[m]
    h = sub.cpg_hit_offset.notna()
    offs = sub.cpg_hit_offset.dropna().astype(int).value_counts().sort_index().to_dict()
    print("  is_reverse=" + str(rev) + ": n=" + str(int(m.sum()))
          + "  hit=" + str(round(float(h.mean()) * 100, 2)) + "%  offsets=" + str(offs))

L = L.sort_values(["chm13_chrom", "chm13_pos", "target_id"]).reset_index(drop=True)
out = OUT_DIR + "/" + TAG + ".tsv.gz"
tmp = out + ".tmp." + str(os.getpid())
L.to_csv(tmp, sep="\t", index=False, float_format="%.6g", compression="gzip")
os.replace(tmp, out)
print("  wrote " + out + "  (" + str(len(L)) + " rows)")
