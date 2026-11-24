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
FAI_DIR = PROJ + "/subproject_methylation/raw/asm_fai"
OUT_DIR = UNIT + "/output/v2_lift/per_hap"

MIN_MAPQ = 20
VERIFY_OFFSETS = (-1, 1)
CTRL_OFFSETS = (5003, -5003)

os.makedirs(OUT_DIR, exist_ok=True)

TARGETS = pd.read_csv(UNIT + "/output/targets/targets_v2_lift.tsv", sep="\t")
assert TARGETS.target_class.notna().all(), "target_class has NaN -- NA-string collision"
assert set(TARGETS.target_class.unique()) == {"discovery", "matched_null"}
assert TARGETS.cpg_key.is_unique, "targets_v2_lift.tsv is not one row per CpG"
HMAP = pd.read_csv(SUB + "/output/mqtl/bam_hap_map.tsv", sep="\t")
CHEM = pd.read_csv(AUD + "/output/cohort_corrected_final.tsv", sep="\t")
CHEM = CHEM[CHEM.chemistry.isin(["R941", "R1041"])]


def ref_to_query_stored(aln, p, tag):
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
            raise ValueError("unexpected CIGAR op " + str(op) + " in " + tag)
    return None, "WALK_OVERRAN"


def lift_one(TAG):
    out_path = OUT_DIR + "/" + TAG + ".tsv.gz"
    if os.path.exists(out_path):
        print(TAG + ": already present, skipping")
        return "SKIP"

    row = HMAP[HMAP.meth_tag == TAG]
    if not len(row):
        print("[UNRESOLVED] " + TAG + ": absent from bam_hap_map.tsv")
        return "UNRESOLVED"
    bam_tag = row.bam_tag.iloc[0]
    if not isinstance(bam_tag, str) or not bam_tag:
        print("[UNRESOLVED] " + TAG + ": haplotype identity unresolved")
        return "UNRESOLVED"

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
    ambiguous = set(int(l) for l in len_counts[len_counts > 1].index)
    pansn = set(c.rsplit("#", 1)[0] + "#" for c in fai.contig if "#" in c)
    if len(pansn) != 1:
        raise ValueError(TAG + ": .fai carries " + str(len(pansn)) + " PanSN prefixes")
    prefix = sorted(pansn)[0]

    sc = CHEM[CHEM["sample"] == sample].chemistry.unique()
    chemistry = "R1041" if "R1041" in sc else ("R941" if "R941" in sc else "NA")

    bam = pysam.AlignmentFile(bam_path, "rb")
    refs = set(bam.references)
    idx_total = sum(s.mapped + s.unmapped for s in bam.get_index_statistics())
    if idx_total == 0:
        raise ValueError(TAG + ": BAM index reports zero records")

    lifted = []
    prefix_mismatch = 0
    for t in TARGETS.itertuples(index=False):
        rec = dict(meth_tag=TAG, sample=sample, chemistry=chemistry, bam_tag=bam_tag,
                   cpg_key=t.cpg_key, target_id=t.target_id, target_class=t.target_class,
                   chrom=t.chrom, plus_c_pos=int(t.plus_c_pos),
                   lift_status="NO_ALIGNMENT", n_aln_covering=0, bam_query="",
                   bed_contig="", contig_pos=-1, anchor_pos=-1, is_reverse=False,
                   mapq=-1, is_supplementary=False, contig_len=-1)
        if t.chrom not in refs:
            rec["lift_status"] = "CHROM_ABSENT_FROM_BAM"
            lifted.append(rec)
            continue
        p = int(t.plus_c_pos)
        cands = [a for a in bam.fetch(t.chrom, p, p + 1)
                 if not a.is_unmapped and a.cigartuples is not None
                 and a.mapping_quality >= MIN_MAPQ]
        rec["n_aln_covering"] = len(cands)
        if not cands:
            lifted.append(rec)
            continue
        cands.sort(key=lambda a: (a.is_supplementary or a.is_secondary,
                                  -a.mapping_quality,
                                  -(a.reference_end - a.reference_start),
                                  a.query_name, a.reference_start))
        aln = cands[0]
        if "#" in aln.query_name and not aln.query_name.startswith(prefix):
            prefix_mismatch += 1
            rec["lift_status"] = "PANSN_PREFIX_MISMATCH"
            rec["bam_query"] = aln.query_name
            lifted.append(rec)
            continue
        qpos_stored, status = ref_to_query_stored(aln, p, TAG)
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
        cpos = int(L_full - 1 - qpos_full) if aln.is_reverse else int(qpos_full)
        rec["contig_pos"] = cpos
        rec["anchor_pos"] = cpos - 1 if aln.is_reverse else cpos
        if L_full in ambiguous:
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
                         + " covered targets outside PanSN prefix " + prefix)

    L = pd.DataFrame(lifted)
    n_ok = int((L.lift_status == "LIFTED").sum())

    want = {}
    for r in L[L.lift_status == "LIFTED"].itertuples(index=False):
        if not r.bed_contig:
            continue
        want.setdefault((r.bed_contig, r.anchor_pos), []).append((r.cpg_key, "anchor", 0))
        for off in VERIFY_OFFSETS:
            want.setdefault((r.bed_contig, r.anchor_pos + off), []).append(
                (r.cpg_key, "verify", off))
        for off in CTRL_OFFSETS:
            want.setdefault((r.bed_contig, r.anchor_pos + off), []).append(
                (r.cpg_key, "ctrl", off))
    want_contigs = set(c for c, _ in want)

    found = {}
    n_rows = n_m = 0
    with open(bed_path, "rt") as f:
        for line in f:
            i = line.find("\t")
            if line[:i] not in want_contigs:
                continue
            n_rows += 1
            parts = line.rstrip("\n").split("\t")
            if parts[3] != "m":
                continue
            n_m += 1
            key = (parts[0], int(parts[1]))
            if key in want:
                found[key] = (int(parts[4]), int(parts[9]), float(parts[10]),
                              int(parts[11]), int(parts[12]), int(parts[13]))

    cols = {k: [] for k in ("score", "n_valid", "pct_mod", "n_mod", "n_canon", "n_other",
                            "hit_anchor", "hit_m1", "hit_p1", "hit_ctrl")}
    for r in L.itertuples(index=False):
        a = found.get((r.bed_contig, r.anchor_pos)) if r.bed_contig else None
        cols["hit_anchor"].append(int(a is not None))
        cols["score"].append(a[0] if a else -1)
        cols["n_valid"].append(a[1] if a else -1)
        cols["pct_mod"].append(a[2] if a else np.nan)
        cols["n_mod"].append(a[3] if a else -1)
        cols["n_canon"].append(a[4] if a else -1)
        cols["n_other"].append(a[5] if a else -1)
        if r.bed_contig:
            cols["hit_m1"].append(int((r.bed_contig, r.anchor_pos - 1) in found))
            cols["hit_p1"].append(int((r.bed_contig, r.anchor_pos + 1) in found))
            cols["hit_ctrl"].append(int(any((r.bed_contig, r.anchor_pos + o) in found
                                            for o in CTRL_OFFSETS)))
        else:
            cols["hit_m1"].append(0)
            cols["hit_p1"].append(0)
            cols["hit_ctrl"].append(0)
    for k, v in cols.items():
        L[k] = v

    lif = L[L.lift_status == "LIFTED"]
    ha = int(lif.hit_anchor.sum())
    print(TAG + ": lifted " + str(n_ok) + "/" + str(len(L))
          + "  anchor hits " + str(ha)
          + " (" + str(round(100.0 * ha / max(len(lif), 1), 1)) + "%)"
          + "  verify -1/+1 " + str(int(lif.hit_m1.sum())) + "/" + str(int(lif.hit_p1.sum()))
          + "  ctrl " + str(int(lif.hit_ctrl.sum()))
          + "  [" + chemistry + ", " + ("label INVERTED" if bam_tag != TAG else "label agrees") + "]")

    L = L.sort_values(["chrom", "plus_c_pos"]).reset_index(drop=True)
    tmp = out_path + ".tmp." + str(os.getpid())
    L.to_csv(tmp, sep="\t", index=False, compression="gzip", float_format="%.6g")
    os.replace(tmp, out_path)
    return "OK"


if __name__ == "__main__":
    tags = sys.argv[1:]
    assert tags, "usage: lift_targets.py <METH_TAG> [...]"
    print("targets " + str(len(TARGETS)) + " ("
          + str(dict(TARGETS.target_class.value_counts())) + ")  tags " + str(len(tags)))
    res = {}
    for tg in tags:
        res[tg] = lift_one(tg)
    print("\nsummary: " + str(pd.Series(list(res.values())).value_counts().to_dict()))
