#!/usr/bin/env python
import os, re, glob
import numpy as np
import pandas as pd
from collections import defaultdict
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

D = CENSV_ROOT + "/analyses/cdr_calls"
OUT_DIR = CENSV_ROOT + "/analyses/cdr_qtl/output"
GLOBE_AA = f"{D}/inputs/globe_aa"
CORR = f"{D}/output/calls"
COHORT = f"{D}/output/cohort_corrected_final.tsv"

ROW_RE = re.compile(r"^(.+?)_(chr[0-9XY]+)_((.+?)#([12])#[^:]+):")


def build_globe_map():
    m = defaultdict(list)
    for bed in sorted(glob.glob(f"{GLOBE_AA}/*_AS-HOR_stv_row.ort.bed")):
        with open(bed) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                mm = ROW_RE.match(p[0])
                if not mm:
                    continue
                sample, chrom, contig, hapnum = mm.group(1), mm.group(2), mm.group(3), mm.group(5)
                asm_hap = f"hap{hapnum}"
                try:
                    s, e = int(p[1]), int(p[2])
                except ValueError:
                    continue
                m[(sample, asm_hap, chrom)].append((contig, s, e))
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
            contig = p[0]
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            extra = p[3:]
            out[contig].append((s, e, extra))
    return out


def main():
    chem = {}
    coh = pd.read_csv(COHORT, sep="\t")
    for _, r in coh.iterrows():
        chem[r["sample"]] = r["chemistry"]

    gmap = build_globe_map()
    print(f"[globe] {len(gmap)} (sample,asm_hap,chrom) groups")

    cdrs_by_contig = {}
    bins_by_contig = {}
    n_tags = 0
    for tagdir in sorted(glob.glob(f"{CORR}/*/")):
        cp = f"{tagdir}/cdrs.bed"
        bp = f"{tagdir}/binned_freq.bed"
        if not os.path.exists(cp):
            continue
        n_tags += 1
        for ci, lst in load_bed_by_contig(cp).items():
            cdrs_by_contig.setdefault(ci, []).extend(lst)
        for ci, lst in load_bed_by_contig(bp, ncol_min=4).items():
            bins_by_contig.setdefault(ci, []).extend(lst)
    print(f"[data] pooled {n_tags} corrected tags; "
          f"{len(cdrs_by_contig)} contigs with CDRs, {len(bins_by_contig)} with windows")

    rows = []
    n_missing_corr = 0
    for (sample, asm_hap, chrom), intervals in sorted(gmap.items()):
        if True:
            cdrs = cdrs_by_contig
            bins = bins_by_contig
            bc = chem.get(sample, "NA")
            if not any(ci in bins_by_contig for (ci, s, e) in intervals):
                n_missing_corr += 1
                continue
            contigs = sorted(set(ci for (ci, s, e) in intervals))
            aa_s = min(s for (ci, s, e) in intervals)
            aa_e = max(e for (ci, s, e) in intervals)
            aa_len = aa_e - aa_s
            if aa_len <= 0:
                continue
            aa_mid = (aa_s + aa_e) / 2.0
            aa_half = aa_len / 2.0

            cdr_list = []
            for (ci, s, e) in intervals:
                for (cs, ce, _x) in cdrs.get(ci, []):
                    if ce > s and cs < e:
                        cdr_list.append((max(cs, s), min(ce, e)))
            win_pct = []
            win_cov = []
            for (ci, s, e) in intervals:
                for (ws, we, extra) in bins.get(ci, []):
                    if we > s and ws < e and len(extra) >= 1:
                        try:
                            pct = float(extra[0])
                        except ValueError:
                            continue
                        win_pct.append(pct)
                        if len(extra) >= 2:
                            try:
                                win_cov.append(float(extra[1]))
                            except ValueError:
                                win_cov.append(np.nan)
            n_win = len(win_pct)
            mean_5mc = float(np.mean(win_pct)) if n_win else np.nan
            if n_win and len(win_cov) == n_win and np.nansum(win_cov) > 0:
                wc = np.array(win_cov); wp = np.array(win_pct)
                ok = ~np.isnan(wc)
                mean_5mc_cw = float(np.sum(wp[ok] * wc[ok]) / np.sum(wc[ok])) if ok.any() else np.nan
            else:
                mean_5mc_cw = mean_5mc

            n = len(cdr_list)
            if n > 0:
                widths = np.array([e - s for s, e in cdr_list], dtype=float)
                mids = np.array([(s + e) / 2.0 for s, e in cdr_list], dtype=float)
                total_bp = float(widths.sum())
                rel = (mids - aa_mid) / aa_half
                pos_offset = float(np.median(rel))
                pos_var = float(np.var(rel))
                pos_iqr = float(np.percentile(rel, 75) - np.percentile(rel, 25))
                med_w = float(np.median(widths))
            else:
                total_bp = 0.0
                pos_offset = pos_var = pos_iqr = med_w = np.nan

            rows.append(dict(
                sample=sample, hap_local=asm_hap, basecaller=bc,
                chm13_chrom=chrom, asm_contig=",".join(contigs),
                active_start=aa_s, active_end=aa_e, active_len=aa_len,
                n_cdrs=n, total_cdr_bp=total_bp,
                cdrs_per_mb=n / (aa_len / 1e6),
                cdr_fraction=total_bp / aa_len,
                median_cdr_width=med_w,
                cdr_pos_offset=pos_offset, cdr_pos_var=pos_var, cdr_pos_iqr=pos_iqr,
                mean_5mC=mean_5mc, mean_5mC_covwt=mean_5mc_cw, n_windows=n_win,
            ))

    df = pd.DataFrame(rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    out = f"{OUT_DIR}/per_hap_chrom_phenotypes.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.6g")
    print(f"[written] {out}: {len(df)} (sample,hap,chrom) rows; "
          f"{df['sample'].nunique()} samples; {n_missing_corr} haps missing corrected calls")
    g = df.groupby("chm13_chrom").agg(
        n_haps=("sample", "count"),
        med_n_cdrs=("n_cdrs", "median"),
        med_cdrs_per_mb=("cdrs_per_mb", "median"),
        med_mean_5mC=("mean_5mC", "median"),
        med_pos_offset=("cdr_pos_offset", "median"),
    )
    print(g.to_string())
    print(f"\n[overall] median n_cdrs/hap-chrom = {df['n_cdrs'].median()}; "
          f"median mean_5mC = {df['mean_5mC'].median():.2f}%")


if __name__ == "__main__":
    main()
