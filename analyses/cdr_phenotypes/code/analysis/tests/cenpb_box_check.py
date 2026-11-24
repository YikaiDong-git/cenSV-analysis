#!/usr/bin/env python
"""CENP-B box test for the chr8 CDR lead deletion."""
import os
import re
import subprocess
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
UNIT = ROOT + "/analyses/cdr_phenotypes"
FA = ROOT + "/resources/chm13v2.0.fa.gz"
CENSAT = ROOT + "/resources/chm13v2.0_censat_v2.1.active.bed"
SAMTOOLS = ROOT + "/.pixi/envs/default/bin/samtools"
VARIANTS = ROOT + "/analyses/cdr_qtl_full_panel/geno/full_variants.tsv"
SCAN = UNIT + "/output/agg/rescan_all_frames.tsv.gz"
OUT = UNIT + "/output/cenpb"

LEAD = "chr8_sv30_DEL_mon_best8_midpoint1000"

FWD = re.compile(r"(?=([ACGT]TTCG[ACGT]{4}A[ACGT]{2}CGGG[ACGT]))")
REV = re.compile(r"(?=([ACGT]CCCG[ACGT]{2}T[ACGT]{4}CGAA[ACGT]))")
INFORM_F = [(1, "T"), (2, "T"), (3, "C"), (4, "G"), (9, "A"),
            (12, "C"), (13, "G"), (14, "G"), (15, "G")]
INFORM_R = [(1, "C"), (2, "C"), (3, "C"), (4, "G"), (7, "T"),
            (12, "C"), (13, "G"), (14, "A"), (15, "A")]


def fetch(region):
    txt = subprocess.run([SAMTOOLS, "faidx", FA, region], check=True,
                         capture_output=True, text=True).stdout
    return "".join(txt.split("\n")[1:]).upper()


def exact_hits(seq):
    f = [m.start() for m in FWD.finditer(seq)]
    r = [m.start() for m in REV.finditer(seq)]
    return f, r


def mismatch_hits(seq, max_mm):
    arr = np.frombuffer(seq.encode(), dtype="S1")
    n = len(seq)
    out = {}
    for tag, spec in (("+", INFORM_F), ("-", INFORM_R)):
        mm = np.zeros(n - 16, dtype=np.int8)
        for p, b in spec:
            mm += (arr[p:n - 16 + p] != b.encode()).astype(np.int8)
        idx = np.flatnonzero(mm <= max_mm)
        out[tag] = list(zip(idx.tolist(), mm[idx].tolist()))
    return out["+"], out["-"]


def selfcheck():
    box = "ATTCGTTGGAAACGGGA"
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    rc = "".join(comp[b] for b in reversed(box))
    assert FWD.search(box) and not FWD.search(rc), "forward motif regex is wrong"
    assert REV.search(rc) and not REV.search(box), "reverse motif regex is wrong"
    print("instrument check: motif regexes PASS (fwd fires on box only, rev on revcomp only)")


def main():
    os.makedirs(OUT, exist_ok=True)
    selfcheck()

    ncol = len(open(CENSAT).readline().rstrip("\n").split("\t"))
    names = ["chrom", "start", "end", "name", "score", "strand",
             "tstart", "tend", "rgb"][:ncol]
    cen = pd.read_csv(CENSAT, sep="\t", header=None, names=names)
    c8 = cen[cen.chrom == "chr8"]
    assert len(c8) >= 1, "no chr8 row in " + CENSAT
    a_start, a_end = int(c8.start.min()), int(c8.end.max())
    print("chr8 active array (CenSat v2.1): chr8:%d-%d  (%.3f Mb)"
          % (a_start, a_end, (a_end - a_start) / 1e6))
    for _, r in c8.iterrows():
        print("   row  %d-%d  %s" % (r.start, r.end, r.get("name", "")))

    seq = fetch("chr8:%d-%d" % (a_start + 1, a_end))
    assert len(seq) == a_end - a_start, (len(seq), a_end - a_start)
    n_amb = sum(1 for b in seq if b not in "ACGT")
    print("fetched %d bp hull, %d non-ACGT" % (len(seq), n_amb))

    in_row = np.zeros(len(seq), dtype=bool)
    rows0 = []
    for _, rr in c8.iterrows():
        s0, e0 = int(rr.start) - a_start, int(rr.end) - a_start
        in_row[s0:e0] = True
        rows0.append((s0, e0))
    hull_bp, row_bp = len(seq), int(in_row.sum())
    print("hull %d bp, annotated rows %d bp, inter-row gap %d bp (%.2f%% of hull)"
          % (hull_bp, row_bp, hull_bp - row_bp, 100.0 * (hull_bp - row_bp) / hull_bp))

    f, r = exact_hits(seq)
    all_boxes = sorted([(p, "+") for p in f] + [(p, "-") for p in r])
    boxes = [(p, st) for p, st in all_boxes if in_row[p] and in_row[min(p + 16, len(seq) - 1)]]
    nf = sum(1 for _, st in boxes if st == "+")
    nr = len(boxes) - nf
    print("")
    print("EXACT CENP-B boxes inside annotated chr8 HOR rows: %d  (+%d / -%d)"
          % (len(boxes), nf, nr))
    print("   (%d of %d hull boxes fell in the inter-row gap and are excluded)"
          % (len(all_boxes) - len(boxes), len(all_boxes)))
    print("   density: 1 per %.0f bp   ( = 1 per %.2f alpha-satellite monomers of 171 bp)"
          % (row_bp / max(len(boxes), 1), row_bp / max(len(boxes), 1) / 171.0))

    mf0, mr0 = mismatch_hits(seq, 1)
    mf = [(p, m) for p, m in mf0 if in_row[p]]
    mr = [(p, m) for p, m in mr0 if in_row[p]]
    print("<=1 mismatch over the 9 informative positions: %d  (+%d / -%d)"
          % (len(mf) + len(mr), len(mf), len(mr)))

    covered = np.zeros(len(seq), dtype=bool)
    for p, _ in boxes:
        covered[p:p + 17] = True
    print("   %.2f%% of annotated-row bases sit inside an exact box"
          % (100.0 * (covered & in_row).sum() / row_bp))

    vt = pd.read_csv(VARIANTS, sep="\t", usecols=["sv_locus", "chrom", "pos", "end",
                                                  "svtype", "svlen", "af"])
    lead = vt[vt.sv_locus == LEAD]
    assert len(lead) == 1, "lead not unique: %d" % len(lead)
    lead = lead.iloc[0]
    print("")
    print("LEAD  %s" % LEAD)
    print("   chr8:%d-%d  %s  %d bp  AF %.3f"
          % (lead.pos, lead.end, lead.svtype, lead.svlen, lead.af))

    del_s0 = int(lead.pos) - a_start
    del_e0 = int(lead.end) - a_start
    assert 0 <= del_s0 < del_e0 <= len(seq), (del_s0, del_e0, len(seq))
    deleted = seq[del_s0:del_e0]
    print("   deleted sequence: %d bp  (= %.2f alpha-satellite monomers)"
          % (len(deleted), len(deleted) / 171.0))

    inside = [(p, s) for p, s in boxes if p >= del_s0 and p + 17 <= del_e0]
    touch = [(p, s) for p, s in boxes if p < del_e0 and p + 17 > del_s0]
    print("")
    print("   EXACT boxes fully inside the deleted interval: %d" % len(inside))
    for p, s in inside:
        print("      chr8:%d  strand %s  %s" % (a_start + p + 1, s, seq[p:p + 17]))
    print("   EXACT boxes overlapping the deleted interval (incl. partial): %d" % len(touch))
    for p, s in touch:
        if (p, s) not in inside:
            print("      chr8:%d  strand %s  %s  (partial)" % (a_start + p + 1, s, seq[p:p + 17]))

    if boxes:
        d = min(abs((a_start + p) - int(lead.pos)) for p, _ in boxes)
        print("   distance from the deletion start to the nearest exact box: %d bp" % d)

    mm_in = [(p, m, "+") for p, m in mf if p >= del_s0 and p + 17 <= del_e0] + \
            [(p, m, "-") for p, m in mr if p >= del_s0 and p + 17 <= del_e0]
    print("   <=1-mismatch boxes fully inside: %d" % len(mm_in))
    for p, m, s in sorted(mm_in):
        print("      chr8:%d  strand %s  %d mm  %s" % (a_start + p + 1, s, m, seq[p:p + 17]))

    L = int(lead.svlen)
    counts = np.bincount([p for p, _ in boxes], minlength=len(seq))
    csum = np.concatenate([[0], np.cumsum(counts)])
    row_csum = np.concatenate([[0], np.cumsum(in_row.astype(np.int64))])
    starts = np.arange(0, len(seq) - L)
    fully_in_row = (row_csum[starts + L] - row_csum[starts]) == L
    starts = starts[fully_in_row]
    hi = np.minimum(starts + L - 17 + 1, len(seq))
    n_in = csum[hi] - csum[starts]
    p_any = float((n_in >= 1).mean())
    obs = len(inside)
    print("")
    print("CALIBRATION over every %d bp window lying FULLY inside an annotated row"
          " (n = %d windows):" % (L, len(starts)))
    print("   P(window fully contains >=1 exact box) = %.4f" % p_any)
    print("   mean boxes per window = %.3f" % n_in.mean())
    if obs > 0:
        print("   one-sided empirical P(random window contains >= %d) = %.4f"
              % (obs, float((n_in >= obs).mean())))
    else:
        print("   observed 0 boxes; empirical P(random window contains 0) = %.4f"
              % float((n_in == 0).mean()))
    if p_any > 0.5:
        print("   -> containing a box is the COMMON case; presence alone is NOT evidence.")

    sc = pd.read_csv(SCAN, sep="\t")
    sc = sc[(sc.chrom == "chr8") & (sc.frame == "splice") & (sc.cohort == "all")]
    strength = sc.groupby("sv_locus").abs_t_obs.max().rename("max_abs_t")
    v8 = vt[(vt.chrom == "chr8") & (vt.svtype == "DEL") &
            (vt.pos >= a_start) & (vt.end <= a_end)].copy()
    v8 = v8.join(strength, on="sv_locus")
    v8 = v8[v8.max_abs_t.notna()].copy()
    print("")
    print("chr8 deletions inside the active array with a scan result: %d" % len(v8))

    nb, nbmm = [], []
    mf_pos = set(p for p, _ in mf)
    mr_pos = set(p for p, _ in mr)
    for _, rr in v8.iterrows():
        s0, e0 = int(rr.pos) - a_start, int(rr.end) - a_start
        nb.append(int(csum[min(e0 - 17 + 1, len(seq))] - csum[max(s0, 0)]) if e0 - 17 + 1 > s0 else 0)
        nbmm.append(sum(1 for p in mf_pos if s0 <= p and p + 17 <= e0) +
                    sum(1 for p in mr_pos if s0 <= p and p + 17 <= e0))
    v8["n_box"] = nb
    v8["n_box_mm1"] = nbmm
    v8["removes_box"] = v8.n_box > 0

    strong = v8[v8.max_abs_t > 4]
    weak = v8[v8.max_abs_t <= 4]
    for tag, g in (("|t|>4 ", strong), ("|t|<=4", weak)):
        if len(g):
            print("   %s: %3d deletions, %d remove >=1 box (%.1f%%), median svlen %.0f bp"
                  % (tag, len(g), int(g.removes_box.sum()),
                     100 * g.removes_box.mean(), g.svlen.median()))
        else:
            print("   %s: 0 deletions" % tag)

    if len(strong) and len(weak):
        from scipy.stats import fisher_exact, mannwhitneyu
        tab = [[int(strong.removes_box.sum()), int((~strong.removes_box).sum())],
               [int(weak.removes_box.sum()), int((~weak.removes_box).sum())]]
        orr, pf = fisher_exact(tab)
        print("   Fisher exact, removes-a-box x strong/weak: OR = %.3f, P = %.4f" % (orr, pf))
        print("   NOTE longer deletions remove more boxes mechanically; box RATE per kb:")
        s_rate = strong.n_box / (strong.svlen / 1000.0)
        w_rate = weak.n_box / (weak.svlen / 1000.0)
        u, pu = mannwhitneyu(s_rate, w_rate, alternative="two-sided")
        print("      strong median %.2f boxes/kb, weak median %.2f boxes/kb, MWU P = %.4f"
              % (s_rate.median(), w_rate.median(), pu))

    v8.sort_values("max_abs_t", ascending=False).to_csv(
        OUT + "/chr8_del_cenpb.tsv", sep="\t", index=False)
    pd.DataFrame({"pos_chm13": [a_start + p + 1 for p, _ in boxes],
                  "strand": [s for _, s in boxes],
                  "seq": [seq[p:p + 17] for p, _ in boxes]}).to_csv(
        OUT + "/chr8_array_cenpb_boxes.tsv", sep="\t", index=False)
    print("")
    print("[written] " + OUT + "/chr8_del_cenpb.tsv")
    print("[written] " + OUT + "/chr8_array_cenpb_boxes.tsv")


main()
