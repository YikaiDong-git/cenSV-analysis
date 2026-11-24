#!/usr/bin/env python3
import sys
import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
LIB = CENSV_ROOT + "/analyses/finemap_multi_sv/code/lib"
sys.path.insert(0, LIB)
from raug import build_canonical_raug, build_z_aug, psd_safeguard  # noqa

_S = {}


def load(bundle, zmat):
    b = np.load(bundle, allow_pickle=True)
    z = np.load(zmat, allow_pickle=True)
    _S["R"] = b["R_snp_snp"].astype(np.float64)
    _S["B"] = b["B"].astype(np.float64)
    _S["locus"] = [str(x) for x in b["sv_locus"]]
    _S["svidx"] = np.asarray(b["sv_index"])
    _S["cv"] = np.asarray(b["sv_cv_r2"], float)
    _S["z"] = np.asarray(z["z"], np.float64)
    _S["zf"] = np.asarray(z["zf"], np.float64)
    _S["nobs"] = np.asarray(z["nobs"], int)
    _S["traits"] = [str(x) for x in z["traits"]]
    _S["p"], _S["G"], _S["T"] = _S["R"].shape[0], _S["B"].shape[1], len(_S["traits"])
    B = _S["B"]; RB = _S["R"] @ B
    nz = (B != 0).sum(0); quad = np.einsum("ij,ij->j", B, RB)
    _S["keep_imp"] = (nz >= 3) & (quad > 0)
    return [_S["p"], _S["G"], _S["T"]]


def traits():   return _S["traits"]
def nobs():     return _S["nobs"].tolist()
def sv_info(j): return [_S["locus"][j], int(_S["svidx"][j]), float(_S["cv"][j])]

def hb_impute_raug(j):
    if not _S["keep_imp"][j]:
        return None
    Bj = _S["B"][:, [j]]
    R_aug, D = build_canonical_raug(_S["R"], Bj); R_aug, _ = psd_safeguard(R_aug)
    _S["cur"] = (Bj, D)
    return R_aug

def hb_impute_zaug(t):
    Bj, D = _S["cur"]
    return build_z_aug(_S["zf"][:, t], Bj, D)

def hb_matched(j, t):
    z = _S["z"][:, t]; obs = np.where(np.isfinite(z))[0]
    if obs.size < 5:
        return None
    Roo = _S["R"][np.ix_(obs, obs)]; Bj = _S["B"][obs][:, [j]]; zoo = z[obs]
    nz = int((Bj != 0).sum()); quad = float(Bj[:, 0] @ (Roo @ Bj[:, 0]))
    if nz < 3 or quad <= 0:
        return None
    R_aug, D = build_canonical_raug(Roo, Bj); R_aug, _ = psd_safeguard(R_aug)
    return [R_aug, build_z_aug(zoo, Bj, D)]

def aw_impute_raug():
    keep = np.where(_S["keep_imp"])[0]
    if keep.size == 0:
        return None
    Bk = _S["B"][:, keep]
    R_aug, D = build_canonical_raug(_S["R"], Bk); R_aug, _ = psd_safeguard(R_aug)
    _S["cur_aw"] = (Bk, D, keep)
    return R_aug

def aw_impute_zaug(t):
    Bk, D, keep = _S["cur_aw"]
    return build_z_aug(_S["zf"][:, t], Bk, D)

def aw_impute_keep():
    return _S["cur_aw"][2].tolist()

def aw_matched(t):
    z = _S["z"][:, t]; obs = np.where(np.isfinite(z))[0]
    if obs.size < 5:
        return None
    Roo = _S["R"][np.ix_(obs, obs)]; Boo = _S["B"][obs]; zoo = z[obs]
    nz = (Boo != 0).sum(0); quad = np.einsum("ij,ij->j", Boo, Roo @ Boo)
    svk = np.where((nz >= 3) & (quad > 0))[0]
    if svk.size == 0:
        return None
    R_aug, D = build_canonical_raug(Roo, Boo[:, svk]); R_aug, _ = psd_safeguard(R_aug)
    return [R_aug, build_z_aug(zoo, Boo[:, svk], D), svk.tolist()]
