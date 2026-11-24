#!/usr/bin/env python3
"""finemap_group — BATCHED-SCORE engine: fine-map a whole (chrom, stratum) group in ONE window pass."""
from __future__ import annotations

import os
import time

import numpy as np
from scipy.stats import chi2

import finemap_core as core
from finemap_core import (
    ScoreNull, VREL, KEEP_DIR,
    load_keep, load_covar, load_pheno, _svpgen_index, _svpgen_prefixes, _psam_iids,
    read_sv_dosages, list_window_svs, extract_snp_pgen, iter_snp_blocks,
    masked_meanimpute, std_present, _finalize_fit, log,
)


def _load_pred_ids(chrom):
    return core.load_pred_ids(chrom)


def _res(f, status):
    return {"phecode": f["phecode"], "panel_mode": f["panel_mode"], "mcc": f["mcc"],
            "status": status, "out_dir": f["out_dir"]}


def _score_block_vectorized(Dm, R_c, W_c, WX_c, XtWXi_list, K, m):
    b = Dm.shape[0]
    U = Dm @ R_c
    Vraw = (Dm * Dm) @ W_c
    Astack = (Dm @ WX_c).reshape(b, m, K)
    z = np.empty((b, m)); beta = np.empty((b, m)); p = np.empty((b, m))
    for j in range(m):
        Aj = Astack[:, j, :]
        corr = np.einsum("bk,kl,bl->b", Aj, XtWXi_list[j], Aj)
        Vrawj = Vraw[:, j]
        Vj = Vrawj - corr
        good = np.isfinite(Vj) & (Vj > VREL * Vrawj) & (Vrawj > 0)
        Uj = U[:, j]
        Vsafe = np.where(good, Vj, 1.0)
        z[:, j] = np.where(good, Uj / np.sqrt(Vsafe), np.nan)
        beta[:, j] = np.where(good, Uj / Vsafe, np.nan)
        p[:, j] = np.where(good, chi2.sf(np.where(good, Uj * Uj / Vsafe, 0.0), 1), np.nan)
    return z, beta, p


def build_group(chrom, stratum, fits, *, top_assoc=60, top_tag=60, L_snp_cap=200,
                min_callrate=0.90, force=False):
    t_grp = time.time()
    for f in fits:
        f.setdefault("snp_select", "assoc_tag")
    iids, iidx = load_keep(stratum)
    X, _ = load_covar(stratum, iidx)
    N = len(iids)
    K = X.shape[1]
    sv_index, _ = _svpgen_index()
    sv_src = _psam_iids(_svpgen_prefixes()[0] + ".psam")

    phen = {}
    for f in fits:
        key = (f["phecode"], f["mcc"])
        if key in phen:
            continue
        y = load_pheno(f["mcc"], f["phecode"], iidx)
        null = ScoreNull(X, y)
        if not null.tested:
            phen[key] = {"skip": f"below floor (ncase={null.ncase},ncontrol={null.ncontrol})"}
        elif null.degenerate:
            phen[key] = {"skip": "degenerate null (|sum r|>IDENT_TOL)"}
        else:
            phen[key] = {"null": null, "y": y, "present": null.present,
                         "phash": null.present.tobytes()}
    live = [f for f in fits if "skip" not in phen[(f["phecode"], f["mcc"])]]
    log(f"[group {chrom}.{stratum}] {len(fits)} fits, {len(phen)} phecodes, "
        f"{len(live)} live ({len(fits)-len(live)} skipped-pheno)")
    if not live:
        return [_res(f, "SKIP:" + phen[(f["phecode"], f["mcc"])]["skip"]) for f in fits]

    fit_svids = {}
    union = []
    for f in live:
        svids = (list(dict.fromkeys(f["hit_sv_ids"])) if f["panel_mode"] == "hit_only"
                 else list_window_svs(chrom, sv_index))
        fit_svids[id(f)] = svids
        union.extend(svids)
    union = list(dict.fromkeys(union))
    ucol = {s: k for k, s in enumerate(union)}
    Gsv_all = read_sv_dosages(union, iids, iidx, sv_index=sv_index, src_iids=sv_src)

    scratch = os.environ.get("FINEMAP_SCRATCH", "/tmp/finemap_scratch_ext2")
    wpref = f"{scratch}/{stratum}/{chrom}.peri.keep"
    extract_snp_pgen(chrom, f"{KEEP_DIR}/{stratum}.keep.fid_iid.tsv", wpref, force=force)
    pred_ids = _load_pred_ids(chrom)

    classes = {}
    for key, pk in phen.items():
        if "skip" in pk:
            continue
        h = pk["phash"]
        if h not in classes:
            present = pk["present"]
            Gsv_imp, _ = masked_meanimpute(Gsv_all, present)
            classes[h] = {"present": present, "npres": int(present.sum()), "members": [],
                          "Gsv_std": std_present(Gsv_imp, present)}
            del Gsv_imp
        classes[h]["members"].append(key)
    for h, c in classes.items():
        members = c["members"]
        c["mindex"] = {key: j for j, key in enumerate(members)}
        c["R_c"] = np.column_stack([phen[k]["null"].r for k in members])
        c["W_c"] = np.column_stack([phen[k]["null"].w for k in members])
        c["WX_c"] = np.hstack([phen[k]["null"].WX for k in members])
        c["XtWXi_list"] = [phen[k]["null"].XtWXi for k in members]
    fits_by_class = {h: [] for h in classes}
    for f in live:
        h = phen[(f["phecode"], f["mcc"])]["phash"]
        fits_by_class[h].append(f)

    ids_all, pos_all = [], []
    phen_z = {k: [] for k in phen if "skip" not in phen[k]}
    phen_b = {k: [] for k in phen_z}
    phen_p = {k: [] for k in phen_z}
    cls_cr = {h: [] for h in classes}
    fit_mc = {id(f): [] for f in live}

    blk = max(256, min(2048, int(1.5e9 / (max(N, 1) * 8))))
    nblk = 0
    for ids_b, pos_b, D in iter_snp_blocks(wpref, iids, iidx, block=blk):
        ids_all.extend(ids_b); pos_all.extend(int(p) for p in pos_b)
        for h, c in classes.items():
            present = c["present"]; npres = c["npres"]; m = len(c["members"])
            Dm, cr = masked_meanimpute(D, present)
            z, b_, p_ = _score_block_vectorized(Dm, c["R_c"], c["W_c"], c["WX_c"], c["XtWXi_list"], K, m)
            for key in c["members"]:
                j = c["mindex"][key]
                phen_z[key].append(z[:, j]); phen_b[key].append(b_[:, j]); phen_p[key].append(p_[:, j])
            Ds = std_present(Dm, present)
            Csv = (Ds @ c["Gsv_std"].T) / npres
            for f in fits_by_class[h]:
                cols = [ucol[s] for s in fit_svids[id(f)]]
                fit_mc[id(f)].append(np.max(np.abs(Csv[:, cols]), axis=1))
            cls_cr[h].append(cr)
            del Dm, z, b_, p_, Ds, Csv
        nblk += 1
    ids_all = list(ids_all); pos_all = list(pos_all)
    for k in phen_z:
        phen_z[k] = np.concatenate(phen_z[k]); phen_b[k] = np.concatenate(phen_b[k]); phen_p[k] = np.concatenate(phen_p[k])
    for h in cls_cr:
        cls_cr[h] = np.concatenate(cls_cr[h])
    for f in live:
        fit_mc[id(f)] = np.concatenate(fit_mc[id(f)])
    log(f"[group {chrom}.{stratum}] streamed {len(ids_all)} SNPs x {nblk} blocks "
        f"({time.time()-t_grp:.0f}s); finalizing {len(live)} fits")

    for _c in classes.values():
        for _k in ("Gsv_std", "R_c", "W_c", "WX_c"):
            _c.pop(_k, None)

    ispred_all = np.array([sid in pred_ids for sid in ids_all], dtype=bool)
    results = []
    for f in fits:
        key = (f["phecode"], f["mcc"])
        if "skip" in phen[key]:
            results.append(_res(f, "SKIP:" + phen[key]["skip"])); continue
        t0 = time.time()
        try:
            os.makedirs(f["out_dir"], exist_ok=True)
            null = phen[key]["null"]; present = phen[key]["present"]
            zc = phen_z[key]; bc = phen_b[key]; pc = phen_p[key]
            crc = cls_cr[phen[key]["phash"]]; mc = fit_mc[id(f)]
            snp_recs = [(ids_all[k], pos_all[k], float(zc[k]), float(bc[k]), float(pc[k]),
                         float(mc[k]), bool(ispred_all[k]), float(crc[k])) for k in range(len(ids_all))]
            svids = fit_svids[id(f)]
            Gsv_fit = Gsv_all[[ucol[s] for s in svids]]
            zsv, _bsv, psv = null.score(Gsv_fit)
            _finalize_fit(
                f["out_dir"], chrom=chrom, phecode=f["phecode"], stratum=stratum, mcc=f["mcc"],
                panel_mode=f["panel_mode"], snp_select=f["snp_select"], sv_ids=svids,
                hit_sv_ids=f["hit_sv_ids"], Gsv=Gsv_fit, zsv=zsv, psv=psv, null=null, X=X,
                y=phen[key]["y"], present=present, iids=iids, iidx=iidx, wpref=wpref, pred_ids=pred_ids,
                snp_recs=snp_recs, nblk=nblk, min_callrate=min_callrate, L_snp_cap=L_snp_cap,
                top_assoc=top_assoc, top_tag=top_tag, t0=t0)
            results.append(_res(f, "built"))
        except Exception as e:   # noqa: BLE001 - keep the group going; record per-fit failure
            results.append(_res(f, "error:" + str(e)[:180]))
    if os.environ.get("FINEMAP_KEEP_EXTRACT") != "1":
        for ext in (".pgen", ".pvar", ".psam", ".build.lock", ".tmp.log", ".log"):
            try:
                os.remove(wpref + ext)
            except OSError:
                pass
    log(f"[group {chrom}.{stratum}] DONE {sum(1 for r in results if r['status']=='built')}/{len(fits)} "
        f"({time.time()-t_grp:.0f}s)")
    return results
