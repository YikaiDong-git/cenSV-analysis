#!/usr/bin/env python3
"""Per-SV cross-cohort burden replication of the 2 AoU SV-causal fine-map candidates in FinnGen R12."""
import sys
import math
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

MS = Path(CENSV_ROOT + "/analyses/finemap_multi_sv")
sys.path.insert(0, str(MS / "code" / "02_prepare_jwbundle_multisv"))
from prepare_jwbundle_multisv import harmonize_gwas  # noqa: E402  (also bootstraps lib path)

GWS = 5e-8

SVREG = {
    "chr5_sv324": dict(
        bundle=MS / "output/bundles/chr5/chr5_super1_43983397_52838788_bundle.npz",
        locus="chr5_sv324_DEL_non_best8_edge500", chrom="chr5", rs=43983397, re=52838788,
        target_panel="vertigo"),
    "chr8_sv5": dict(
        bundle=MS / "output/bundles/chr8/chr8_super1_41863936_49169911_bundle.npz",
        locus="chr8_sv5_INS_non_best8_midpoint500", chrom="chr8", rs=41863936, re=49169911,
        target_panel="hearing"),
}


def two_sided_p(z):
    if not np.isfinite(z):
        return float("nan")
    return math.erfc(abs(z) / math.sqrt(2.0))


def load_sv(svkey):
    d = SVREG[svkey]
    b = np.load(d["bundle"], allow_pickle=False)
    svl = b["sv_locus"].astype(str)
    cols = np.where(svl == d["locus"])[0]
    if len(cols) != 1:
        raise ValueError(f"{svkey}: expected 1 matching sv_locus, got {len(cols)}")
    col = int(cols[0])
    R = b["R_snp_snp"]
    w = b["B"][:, col].astype(np.float64)
    snp_meta = {k: b[k] for k in ("snp_id", "snp_chrom", "snp_pos", "snp_ref", "snp_alt", "snp_maf_panel", "snp_swap")}
    return dict(d, col=col, R=R, w=w, snp_meta=snp_meta,
                cv_r2=float(b["sv_cv_r2"][col]), n_tags_total=int((np.abs(w) > 1e-12).sum()))


def stream_region_rows(tsv_path, regions):
    buf = {c: [] for c in regions}
    with open(tsv_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            snp = f[0]
            p = snp.split("_")
            if len(p) < 4:
                continue
            chrom = p[0]
            reg = regions.get(chrom)
            if reg is None:
                continue
            pos = int(p[1])
            if pos < reg[0] or pos > reg[1]:
                continue
            buf[chrom].append((snp, float(f[3]), float(f[4]), float(f[5]), pos))
    out = {}
    for c, rows in buf.items():
        out[c] = pd.DataFrame(rows, columns=["SNP", "BETA", "SE", "PVALUE", "POS"])
    return out


def per_sv_stats(sv, region_df):
    z_snp, keep = harmonize_gwas(sv["snp_meta"], region_df)
    w = sv["w"]
    tag = (np.abs(w) > 1e-12) & keep
    n_tag_match = int(tag.sum())
    if n_tag_match >= 1:
        idx = np.where(tag)[0]
        wt = w[idx]
        zt = z_snp[idx]
        Rt = sv["R"][np.ix_(idx, idx)]
        denom = float(wt @ Rt @ wt)
        zsv = float(wt @ zt) / math.sqrt(denom) if denom > 0 else float("nan")
        aw = np.abs(wt)
        wmean = float((aw * zt).sum() / aw.sum()) if aw.sum() > 0 else float("nan")
    else:
        zsv = float("nan"); wmean = float("nan")
    n_reg = len(region_df)
    if n_reg:
        pv = region_df["PVALUE"].to_numpy()
        min_p = float(np.nanmin(pv))
        n_gws = int((pv < GWS).sum())
        j = int(np.nanargmin(pv))
        top_snp = str(region_df["SNP"].iloc[j]); top_p = float(pv[j])
        top_z = float(region_df["BETA"].iloc[j] / region_df["SE"].iloc[j])
    else:
        min_p = float("nan"); n_gws = 0; top_snp = "NA"; top_p = float("nan"); top_z = float("nan")
    return dict(n_tag_match=n_tag_match, z_SV=zsv, p_SV=two_sided_p(zsv), wmean_z=wmean,
                region_n_snps=n_reg, region_min_p=min_p, region_n_gws=n_gws,
                region_top_snp=top_snp, region_top_p=top_p, region_top_z=top_z)


def main():
    grid_path, out_path = sys.argv[1], sys.argv[2]
    svs = {k: load_sv(k) for k in SVREG}
    regions = {sv["chrom"]: (sv["rs"], sv["re"]) for sv in svs.values()}

    grid = pd.read_csv(grid_path, sep="\t")
    rows = []
    for _, g in grid.iterrows():
        endpoint, panel, role, htsv = g["endpoint"], g["panel"], g["role"], g["harmonized_tsv"]
        if not Path(htsv).exists():
            sys.stderr.write(f"[skip] {endpoint}: harmonized TSV missing {htsv}\n")
            continue
        sys.stderr.write(f"[stream] {endpoint} <- {htsv}\n")
        reg_rows = stream_region_rows(htsv, regions)
        for svkey, sv in svs.items():
            st = per_sv_stats(sv, reg_rows[sv["chrom"]])
            is_target = (panel == sv["target_panel"])
            rows.append(dict(
                sv=svkey, sv_locus=sv["locus"], cv_r2=sv["cv_r2"], chrom=sv["chrom"],
                endpoint=endpoint, panel=panel, role=role,
                cell="TARGET" if is_target else ("NEG_CTRL_unrelated" if panel == "control" else "NEG_CTRL_crosstrait"),
                n_tags_total=sv["n_tags_total"], **st))
    out = pd.DataFrame(rows)
    cols = ["sv", "sv_locus", "cv_r2", "chrom", "endpoint", "panel", "role", "cell",
            "n_tags_total", "n_tag_match", "z_SV", "p_SV", "wmean_z",
            "region_n_snps", "region_min_p", "region_n_gws", "region_top_snp", "region_top_p", "region_top_z"]
    out = out[cols].sort_values(["sv", "cell", "p_SV"], na_position="last")
    out.to_csv(out_path, sep="\t", index=False, float_format="%.4g")
    sys.stderr.write(f"[done] wrote {out_path} ({len(out)} cells)\n")
    lead = out[out["endpoint"].isin(["H8_VERTIGO", "H8_HL_SEN_NAS"])]
    print(lead[["sv", "endpoint", "cell", "n_tag_match", "z_SV", "p_SV", "wmean_z", "region_min_p", "region_n_gws"]].to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
