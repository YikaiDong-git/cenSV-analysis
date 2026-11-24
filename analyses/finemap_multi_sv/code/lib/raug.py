#!/usr/bin/env python3
from __future__ import annotations

import numpy as np


__all__ = ["build_canonical_raug", "build_z_aug", "psd_safeguard"]


def build_canonical_raug(R: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p, p2 = R.shape
    if p != p2:
        raise ValueError(f"R must be square, got {R.shape}")
    p_b, G = B.shape
    if p_b != p:
        raise ValueError(f"B rows={p_b} must equal R dim={p}")

    RB = R @ B
    BRB = B.T @ RB
    quad = np.diag(BRB).copy()

    if (quad <= 0).any():
        bad = np.where(quad <= 0)[0]
        raise ValueError(
            f"Non-positive quadratic form beta_j' R beta_j for j={bad.tolist()} "
            f"(values={quad[bad].tolist()}). Likely a zero weight vector or "
            f"non-PSD R. Investigate before clipping."
        )

    D = np.sqrt(quad)

    T = RB / D[np.newaxis, :]
    G_gg = BRB / (D[:, np.newaxis] * D[np.newaxis, :])
    np.fill_diagonal(G_gg, 1.0)

    R_aug = np.empty((p + G, p + G), dtype=np.float64)
    R_aug[:p, :p] = R
    R_aug[:p, p:] = T
    R_aug[p:, :p] = T.T
    R_aug[p:, p:] = G_gg

    R_aug = 0.5 * (R_aug + R_aug.T)

    return R_aug, D


def build_z_aug(z_snp: np.ndarray, B: np.ndarray, D: np.ndarray) -> np.ndarray:
    p, G = B.shape
    if z_snp.shape != (p,):
        raise ValueError(f"z_snp shape {z_snp.shape} must be ({p},)")
    if D.shape != (G,):
        raise ValueError(f"D shape {D.shape} must be ({G},)")

    z_sv = (B.T @ z_snp) / D

    z_aug = np.empty(p + G, dtype=np.float64)
    z_aug[:p] = z_snp
    z_aug[p:] = z_sv
    return z_aug


def psd_safeguard(R_aug: np.ndarray, tol: float = -1e-6) -> tuple[np.ndarray, dict]:
    eig_min_pre = float(np.linalg.eigvalsh(R_aug)[0])

    info = {
        "eig_min_pre": eig_min_pre,
        "tikhonov_lambda": 0.0,
        "tikhonov_applied": False,
        "n_clipped": 0,
        "eig_min_post": eig_min_pre,
    }

    if eig_min_pre >= tol:
        return R_aug, info

    lam = abs(eig_min_pre) + 1e-6
    n = R_aug.shape[0]
    R_aug = R_aug + lam * np.eye(n)
    d = np.sqrt(np.diag(R_aug))
    R_aug = R_aug / np.outer(d, d)

    mask_high = R_aug > 1.0
    mask_low = R_aug < -1.0
    n_clipped = int(mask_high.sum() + mask_low.sum())
    R_aug = np.clip(R_aug, -0.999999, 0.999999)
    np.fill_diagonal(R_aug, 1.0)

    eig_min_post = float(np.linalg.eigvalsh(R_aug)[0])

    info.update({
        "tikhonov_lambda": lam,
        "tikhonov_applied": True,
        "n_clipped": n_clipped,
        "eig_min_post": eig_min_post,
    })
    return R_aug, info


def _self_test() -> None:
    rng = np.random.default_rng(seed=42)
    p = 50
    A = rng.standard_normal((p, p))
    R = A @ A.T
    d = np.sqrt(np.diag(R))
    R = R / np.outer(d, d)

    beta = rng.standard_normal(p) * 0.3
    B = beta.reshape(-1, 1)

    R_aug, D = build_canonical_raug(R, B)
    Dval = float(D[0])

    assert np.allclose(R_aug[:p, :p], R), "top-left != R"
    expected_topright = (R @ beta) / Dval
    assert np.allclose(R_aug[:p, -1], expected_topright), "top-right != (R beta)/D"
    assert np.isclose(R_aug[-1, -1], 1.0), f"bottom-right != 1 (got {R_aug[-1, -1]})"
    assert np.allclose(R_aug, R_aug.T), "not symmetric"

    z_snp = rng.standard_normal(p)
    z_aug = build_z_aug(z_snp, B, D)
    assert np.allclose(z_aug[:p], z_snp), "z_aug top != z_snp"
    expected_zsv = (beta @ z_snp) / Dval
    assert np.isclose(z_aug[-1], expected_zsv), "z_SV mismatch"

    R_aug_safe, info = psd_safeguard(R_aug)
    assert abs(info["eig_min_pre"]) < 1e-6, f"eig_min_pre should be ~0 (got {info['eig_min_pre']:.6e})"
    assert not info["tikhonov_applied"], "Tikhonov should NOT fire for rank-G numerical noise"

    G = 3
    B3 = rng.standard_normal((p, G)) * 0.3
    R_aug3, D3 = build_canonical_raug(R, B3)
    assert R_aug3.shape == (p + G, p + G)
    assert np.allclose(np.diag(R_aug3[p:, p:]), 1.0), "bottom-right diag != 1"
    assert np.allclose(R_aug3[:p, :p], R), "multi-SV top-left != R"
    assert (np.abs(R_aug3) <= 1.0 + 1e-9).all(), "|R_aug entries| > 1"
    assert np.allclose(R_aug3, R_aug3.T), "multi-SV not symmetric"
    eigs = np.linalg.eigvalsh(R_aug3)
    n_zero = int((np.abs(eigs) < 1e-6).sum())
    assert n_zero == G, f"expected G={G} ~zero eigenvalues, got {n_zero}"

    R_aug3_safe, info3 = psd_safeguard(R_aug3)
    assert not info3["tikhonov_applied"], "Tikhonov should NOT fire for rank-G multi-SV numerical noise"
    assert np.allclose(R_aug3_safe, R_aug3), "no-op safeguard should return identical matrix"

    R_bad = R_aug3.copy()
    R_bad[0, 1] = -0.999
    R_bad[1, 0] = -0.999
    R_bad_safe, info_bad = psd_safeguard(R_bad)
    if np.linalg.eigvalsh(R_bad)[0] < -1e-6:
        assert info_bad["tikhonov_applied"], "Tikhonov SHOULD fire for materially negative eig"
        assert info_bad["eig_min_post"] > 0, "post-ridge eig_min should be >0"

    print("OK: build_canonical_raug + build_z_aug + psd_safeguard self-test PASS")
    print(f"  p={p}, G={G}")
    print(f"  single-SV: eig_min_pre={info['eig_min_pre']:.6e}, no-ridge (canonical preserved)")
    print(f"  multi-SV:  eig_min_pre={info3['eig_min_pre']:.6e}, no-ridge (canonical preserved)")
    print(f"  multi-SV:  n_zero_eigs (expected G={G}) = {n_zero}")
    print(f"  forced bad: eig_min_pre={info_bad['eig_min_pre']:.4e} -> ridge applied: {info_bad['tikhonov_applied']}")


if __name__ == "__main__":
    _self_test()
