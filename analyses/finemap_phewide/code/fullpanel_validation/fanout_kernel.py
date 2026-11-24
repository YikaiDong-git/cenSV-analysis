#!/usr/bin/env python3
"""Fan-out kernel (reticulate-loaded by fanout_run.R)."""
import sys
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

MULTI = CENSV_ROOT + "/analyses/finemap_multi_sv"
NEW   = CENSV_ROOT + "/analyses/finemap_phewide"
sys.path.insert(0, MULTI + "/code/01_build_multisv_bundle")
sys.path.insert(0, MULTI + "/code/lib")
sys.path.insert(0, NEW + "/code")
import build_multisv_bundle as bmb
import build_cells as bc
from raug import build_canonical_raug, build_z_aug, psd_safeguard
import subprocess

_S = {}


def _query_all(chrom, samples, ref_vcf, start, end):
    cmd = [str(bmb.BCFTOOLS), "query", "-r", f"{chrom}:{start}-{end}", "-s", ",".join(samples),
           "-f", "%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(ref_vcf)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    n = len(samples)
    gt = {"0|0": 0., "0/0": 0., "0|1": 1., "1|0": 1., "0/1": 1., "1/0": 1., "1|1": 2., "1/1": 2.}
    rows_gt, meta = [], []
    for ln in proc.stdout.splitlines():
        c = ln.split("\t")
        if len(c) < 5 + n:
            continue
        _id, vch, vpos, vref, valt = c[:5]
        if len(vref) != 1 or len(valt) != 1 or vref not in "ACGT" or valt not in "ACGT":
            continue
        dos = np.array([gt.get(g[:3], np.nan) for g in c[5:5 + n]], dtype=np.float64)
        m = np.isnan(dos)
        if m.all():
            continue
        if m.any():
            dos[m] = np.nanmean(dos)
        af = dos.sum() / (2.0 * n)
        maf = min(af, 1.0 - af)
        if maf < bmb.MAF_MIN:
            continue
        rows_gt.append(dos)
        meta.append({"snp_id": f"{vch.replace('chr','')}:{vpos}:{vref}:{valt}", "snp_chrom": chrom,
                     "snp_pos": int(vpos), "snp_ref": vref, "snp_alt": valt, "maf_panel": maf, "swap": 0})
    return np.array(rows_gt).T, pd.DataFrame(meta).reset_index(drop=True)


def load_region(chrom, rstart, rend):
    samples = bmb.read_pedigree_samples(bmb.KG_EUR_525_PEDIGREE)
    vcf = bmb.ONE_KG_VCF_TEMPLATE.format(chrom=chrom)
    G, meta = _query_all(chrom, samples, vcf, int(rstart), int(rend))
    _S.update(chrom=chrom, G=G, meta=meta, pos=meta["snp_pos"].values,
              ref=meta["snp_ref"].values, alt=meta["snp_alt"].values, zc={})
    return int(G.shape[1])


def _z(trait, bgz, chrom_gwas, lo, hi):
    if trait not in _S["zc"]:
        gw = bc.load_gwas(bgz, chrom_gwas, int(lo), int(hi))
        _S["zc"][trait] = bc.match_z(_S["pos"], _S["ref"], _S["alt"], gw)
    return _S["zc"][trait]


def build_cell(sv_index, sv_locus, trait, bgz, chrom_gwas, lo, hi, arm, pcap=12000):
    wp = bmb.weight_tsv_path(bmb.WEIGHTS_BASE, _S["chrom"], int(sv_index), sv_locus)
    wdf = bmb.load_weight_tsv(wp)
    cv = float(wdf["cv_r2"].iloc[0]) if "cv_r2" in wdf.columns else float("nan")
    Bfull = bmb.build_B_matrix([wdf], _S["meta"])[:, 0]
    supp_mask = Bfull != 0
    z_all = _z(trait, bgz, chrom_gwas, lo, hi)
    az = np.abs(np.nan_to_num(z_all, nan=0.0))
    keep = list(np.where(supp_mask)[0])
    for i in np.argsort(-az):
        if len(keep) >= pcap:
            break
        if not supp_mask[i]:
            keep.append(int(i))
    keep = np.array(sorted(set(keep)))
    Gk = _S["G"][:, keep]
    Rk = bmb.compute_R_snp_snp(Gk)
    Bk = Bfull[keep].reshape(-1, 1)
    zk = z_all[keep]
    finite = np.isfinite(zk)
    if arm == "matched":
        obs = np.where(finite)[0]
        if obs.size < 5:
            return None
        Roo, Boo, zoo = Rk[np.ix_(obs, obs)], Bk[obs, :], zk[obs]
    else:
        zoo = bc.impg_fill(zk, Rk, finite)
        Roo, Boo = Rk, Bk
    nz = int((Boo != 0).sum())
    quad = float(Boo[:, 0] @ (Roo @ Boo[:, 0]))
    if nz < 3 or quad <= 0:
        return None
    R_aug, D = build_canonical_raug(Roo, Boo)
    R_aug, _info = psd_safeguard(R_aug)
    z_aug = build_z_aug(zoo, Boo, D)
    n = int(R_aug.shape[0])
    return {"R_aug": R_aug, "z_aug": z_aug, "n": n, "p_snp": int(Roo.shape[0]),
            "z_SV": float(z_aug[-1]), "cv": cv,
            "supp_matched": int(supp_mask.sum()), "n_supp": int(len(wdf))}
