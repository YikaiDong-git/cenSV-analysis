#!/usr/bin/env python3
"""Cross-external-cohort Stouffer meta (FinnGen + UKB-Neale) of the 2 AoU SV-causal candidates."""
import sys
import math
import numpy as np
import pandas as pd

fg = pd.read_csv(sys.argv[1], sep="\t")
uk = pd.read_csv(sys.argv[2], sep="\t")

def neff(ncase, nctrl):
    return 4.0 / (1.0 / ncase + 1.0 / nctrl)

LEAD = {
    "chr5_sv324": dict(
        fg_endpoint="H8_VERTIGO", fg_neff=neff(21139, 474839),
        uk_endpoint="UKB_20002_1500_vertigo", uk_neff=None),
    "chr8_sv5": dict(
        fg_endpoint="H8_HL_SEN_NAS", fg_neff=neff(44745, 437331),
        uk_endpoint="UKB_2247_1_hearingdiff", uk_neff=None),
}
UKB_NEFF = {"UKB_2247_1_hearingdiff": 326600.0, "UKB_20002_1500_vertigo": None}


def get_z(df, sv, endpoint):
    r = df[(df["sv"] == sv) & (df["endpoint"] == endpoint)]
    if len(r) == 0 or not np.isfinite(r["z_SV"].iloc[0]) or int(r["n_tag_match"].iloc[0]) == 0:
        return None, 0
    return float(r["z_SV"].iloc[0]), int(r["n_tag_match"].iloc[0])


rows = []
for sv, L in LEAD.items():
    zf, nf_tag = get_z(fg, sv, L["fg_endpoint"])
    zu, nu_tag = get_z(uk, sv, L["uk_endpoint"])
    nf = L["fg_neff"]
    nu = UKB_NEFF.get(L["uk_endpoint"])
    contribs = []
    if zf is not None:
        contribs.append((zf, nf, "FinnGen:" + L["fg_endpoint"]))
    if zu is not None and nu is not None:
        contribs.append((zu, nu, "UKB:" + L["uk_endpoint"]))
    num = sum(math.sqrt(n) * z for z, n, _ in contribs)
    den = math.sqrt(sum(n for _, n, _ in contribs)) if contribs else float("nan")
    z_meta = num / den if den and np.isfinite(den) else float("nan")
    p_meta = math.erfc(abs(z_meta) / math.sqrt(2)) if np.isfinite(z_meta) else float("nan")
    zs = [z for z, _, _ in contribs]
    z_eq = sum(zs) / math.sqrt(len(zs)) if zs else float("nan")
    print(f"\n=== {sv} ({L['fg_endpoint']} + {L['uk_endpoint']}) ===")
    print(f"  FinnGen z = {zf}  (Neff={nf:.0f}, tags={nf_tag})")
    print(f"  UKB     z = {zu}  (Neff={nu}, tags={nu_tag})"
          + ("   [UNTESTABLE: 0 tags in UKB variant set]" if zu is None else ""))
    print(f"  --> Stouffer sqrtN-weighted Z_meta = {z_meta:.3f}  p = {p_meta:.3g}"
          f"   | equal-weight Z = {z_eq:.3f}   | n_cohorts = {len(contribs)}")
    rows.append(dict(sv=sv, fg_endpoint=L["fg_endpoint"], fg_z=zf, fg_neff=round(nf),
                     uk_endpoint=L["uk_endpoint"], uk_z=zu, uk_neff=nu,
                     n_cohorts=len(contribs), z_meta_sqrtN=round(z_meta, 4) if np.isfinite(z_meta) else None,
                     p_meta=p_meta, z_meta_equal=round(z_eq, 4) if np.isfinite(z_eq) else None))

pd.DataFrame(rows).to_csv(sys.argv[3], sep="\t", index=False, float_format="%.4g")
print(f"\n[wrote] {sys.argv[3]}")
print("\nNOTE: signs opposing (FinnGen + / UKB -) collapse chr8_sv5 to null under any positive weighting;")
print("chr5_sv324 external meta = FinnGen only (UKB untestable). AoU is a separate single-cohort anchor.")
