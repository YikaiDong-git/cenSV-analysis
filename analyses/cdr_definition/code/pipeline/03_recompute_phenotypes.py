#!/usr/bin/env python
import os
import re
import glob
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

D = CENSV_ROOT + "/analyses/cdr_calls"
OUT = CENSV_ROOT + "/analyses/cdr_definition/output"
GLOBE_AA = D + "/inputs/globe_aa"
CORR = D + "/output/calls"
COHORT = D + "/output/cohort_corrected_final.tsv"

ROW_RE = re.compile(r"^(.+?)_(chr[0-9XY]+)_((.+?)#([12])#[^:]+):")


def _isnum(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def build_globe_map():
    m = defaultdict(lambda: defaultdict(list))
    for bed in sorted(glob.glob(GLOBE_AA + "/*_AS-HOR_stv_row.ort.bed")):
        with open(bed) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                mm = ROW_RE.match(p[0])
                if not mm:
                    continue
                try:
                    s, e = int(p[1]), int(p[2])
                except ValueError:
                    continue
                m[(mm.group(1), "hap" + mm.group(5), mm.group(2))][mm.group(3)].append((s, e))
    return m


def load_bed_by_contig(path, ncol_min=3):
    out = defaultdict(list)
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < ncol_min:
                continue
            try:
                out[p[0]].append((int(p[1]), int(p[2]), p[3:]))
            except ValueError:
                continue
    return out


def merge(iv):
    iv = sorted(iv)
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def cum_covered(x, blocks):
    t = 0
    for a, b in blocks:
        if x >= b:
            t += b - a
        elif x > a:
            t += x - a
            break
        else:
            break
    return t


def frames(m, blocks, aa_s, aa_e, S):
    half = (aa_e - aa_s) / 2.0
    r_span = (m - (aa_s + aa_e) / 2.0) / half if half > 0 else np.nan
    r_block = np.nan
    for a, b in blocks:
        if a <= m <= b:
            h = (b - a) / 2.0
            r_block = (m - (a + b) / 2.0) / h if h > 0 else np.nan
            break
    r_splice = (cum_covered(m, blocks) - S / 2.0) / (S / 2.0) if S > 0 else np.nan
    return r_span, r_block, r_splice


def main(dips_root, out_name):
    os.makedirs(OUT, exist_ok=True)
    gmap = build_globe_map()

    chem = {}
    if os.path.exists(COHORT):
        c = pd.read_csv(COHORT, sep="\t")
        scol = "sample" if "sample" in c.columns else c.columns[0]
        for col in ("basecaller", "chemistry", "is_r1041"):
            if col in c.columns:
                chem = dict(zip(c[scol], c[col]))
                break

    cdrs_by, bins_by = {}, {}
    tags = 0
    for tagdir in sorted(glob.glob(CORR + "/*")):
        if not os.path.isdir(tagdir):
            continue
        tags += 1
        dip_bed = os.path.join(dips_root, os.path.basename(tagdir), "cdrs.bed")
        for ci, lst in load_bed_by_contig(dip_bed).items():
            cdrs_by.setdefault(ci, []).extend(lst)
        for ci, lst in load_bed_by_contig(tagdir + "/binned_freq.bed", 4).items():
            bins_by.setdefault(ci, []).extend(lst)
    print("[data] %d corrected tags; %d contigs with dips, %d with windows"
          % (tags, len(cdrs_by), len(bins_by)))

    rows = []
    n_multicontig = n_nodata = n_singleblock = n_frames_agree = 0
    for (sample, hap, chrom), by_contig in sorted(gmap.items()):
        if len(by_contig) > 1:
            n_multicontig += 1
            continue
        contig, raw = next(iter(by_contig.items()))
        if contig not in bins_by:
            n_nodata += 1
            continue

        blocks = merge(raw)
        aa_s, aa_e = blocks[0][0], blocks[-1][1]
        span = aa_e - aa_s
        S = sum(b - a for a, b in blocks)
        if S <= 0:
            continue

        win = [(ws, we, ex) for (ws, we, ex) in bins_by.get(contig, [])
               if any(we > a and ws < b for a, b in blocks)]

        def min5mc(cs, ce):
            v = [float(ex[0]) for (ws, we, ex) in win
                 if we > cs and ws < ce and ex and _isnum(ex[0])]
            return min(v) if v else np.nan

        dips = []
        for (cs, ce, _x) in cdrs_by.get(contig, []):
            for a, b in blocks:
                if ce > a and cs < b:
                    lo, hi = max(cs, a), min(ce, b)
                    dips.append(((lo + hi) / 2.0, hi - lo, min5mc(lo, hi)))
                    break

        wp = [float(ex[0]) for (_, _, ex) in win if ex and _isnum(ex[0])]
        wc = [float(ex[1]) for (_, _, ex) in win if len(ex) > 1 and _isnum(ex[1])]
        mean5 = float(np.mean(wp)) if wp else np.nan
        if wp and len(wc) == len(wp) and np.sum(wc) > 0:
            mean5cw = float(np.sum(np.array(wp) * np.array(wc)) / np.sum(wc))
        else:
            mean5cw = mean5

        rec = dict(sample=sample, hap_local=hap, chm13_chrom=chrom, asm_contig=contig,
                   basecaller=chem.get(sample, "NA"),
                   active_start=aa_s, active_end=aa_e,
                   active_len_span=span, active_len_sum=S,
                   n_intervals=len(raw), n_blocks=len(blocks), gap_bp=span - S,
                   n_cdrs=len(dips), mean_5mC=mean5, mean_5mC_covwt=mean5cw,
                   n_windows=len(win))

        if dips:
            fr = np.array([frames(m, blocks, aa_s, aa_e, S) for m, _, _ in dips])
            widths = np.array([w for _, w, _ in dips], dtype=float)
            depths = np.array([d for _, _, d in dips], dtype=float)
            total = float(widths.sum())
            j = int(np.nanargmin(depths)) if not np.all(np.isnan(depths)) else None
            for i, tag in enumerate(("span", "block", "splice")):
                col = fr[:, i]
                rec["cdr_pos_offset_" + tag] = float(np.nanmedian(col))
                rec["cdr_pos_iqr_" + tag] = float(np.nanpercentile(col, 75) - np.nanpercentile(col, 25))
                rec["cdr_pos_var_" + tag] = float(np.nanvar(col))
                rec["cdr_pos_offset_deepest_" + tag] = (
                    float(col[j]) if j is not None else float(np.nanmedian(col)))
            rec.update(total_cdr_bp=total, median_cdr_width=float(np.median(widths)),
                       cdr_fraction_span=total / span, cdr_fraction_sum=total / S,
                       cdrs_per_mb_span=len(dips) / (span / 1e6),
                       cdrs_per_mb_sum=len(dips) / (S / 1e6),
                       deepest_min5mC=float(depths[j]) if j is not None else np.nan)
            if len(blocks) == 1:
                n_singleblock += 1
                if (np.allclose(fr[:, 0], fr[:, 1], equal_nan=True)
                        and np.allclose(fr[:, 0], fr[:, 2], equal_nan=True)):
                    n_frames_agree += 1
        else:
            rec.update(total_cdr_bp=0.0, median_cdr_width=np.nan,
                       cdr_fraction_span=0.0, cdr_fraction_sum=0.0,
                       cdrs_per_mb_span=0.0, cdrs_per_mb_sum=0.0, deepest_min5mC=np.nan)
            for tag in ("span", "block", "splice"):
                for k in ("cdr_pos_offset_", "cdr_pos_iqr_", "cdr_pos_var_",
                          "cdr_pos_offset_deepest_"):
                    rec[k + tag] = np.nan
        rows.append(rec)

    df = pd.DataFrame(rows)
    out = os.path.join(OUT, out_name)
    df.to_csv(out, sep="\t", index=False, float_format="%.6g")

    print("[written] %s: %d rows, %d samples" % (out, len(df), df["sample"].nunique()))
    print("[qc] dropped %d cross-contig groups (defect 2); %d without methylation data"
          % (n_multicontig, n_nodata))
    print("[qc] single-block groups with dips: %d; all three frames identical in %d"
          % (n_singleblock, n_frames_agree))

    assert n_singleblock == n_frames_agree, (
        "FRAME INVARIANCE VIOLATED: the three parameterisations must agree exactly where the "
        "array is a single contiguous block. %d of %d disagree."
        % (n_singleblock - n_frames_agree, n_singleblock))
    print("[gate] frame-invariance assertion PASSED")
    assert (df.active_len_sum <= df.active_len_span + 1).all(), "SUM exceeds SPAN somewhere"
    assert df.active_len_sum.max() < 12e6, (
        "an implausible array length survived: %.2f Mb" % (df.active_len_sum.max() / 1e6))
    print("[gate] length sanity assertions PASSED (max SUM = %.2f Mb)"
          % (df.active_len_sum.max() / 1e6))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dips-root", default=CORR,
                    help="directory of <tag>/cdrs.bed defining the call set")
    ap.add_argument("--out-name", default="per_hap_chrom_phenotypes_reproduced.tsv")
    a = ap.parse_args()
    print("[call set] %s" % a.dips_root)
    main(a.dips_root, a.out_name)
