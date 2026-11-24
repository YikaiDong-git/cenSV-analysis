#!/usr/bin/env python3
"""Artefact leverage of the sign-concordance endpoint, per candidate stratum."""
import os
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
OUT = UNIT + "/output/null_design"
os.makedirs(OUT, exist_ok=True)

st = pd.read_csv(UNIT + "/output/targets/prereg_strata.tsv", sep="\t")
pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
st = st.merge(pairs[["pair_id", "source_target_distance_bp"]], on="pair_id", how="left")


def b_break(n, f_a, alpha=0.05):
    k = int(np.ceil(stats.binom.isf(alpha, n, 0.5))) + 1
    k = min(k, n)
    c_needed = k / n
    slope = 2.0 * f_a - 1.0
    if slope <= 1e-12:
        return np.nan, c_needed, k
    b = (c_needed - (1.0 - f_a)) / slope
    return b, c_needed, k


rows = []
strata = {
    "465 all (pre-registered)": st,
    "465 recoverable": st[st.recoverable_pilot8],
    "84 priority": st[st.priority84],
    "84 recoverable": st[st.priority84 & st.recoverable_pilot8],
    "465 recoverable, cis <5Mb": st[st.recoverable_pilot8
                                    & st.dist_class.isin(["<=100kb", "0.1-1Mb", "1-5Mb"])],
}
for name, sub in strata.items():
    if not len(sub):
        continue
    n = len(sub)
    f_a = float((sub.cohortA_sign < 0).mean())
    bb, c_need, k = b_break(n, f_a)
    rows.append(dict(
        stratum=name, n_pairs=n, n_clusters=int(sub.sv_id.nunique()),
        f_A=round(f_a, 5),
        leverage_2fA_minus_1=round(2 * f_a - 1, 4),
        c_at_b_0p50=round((1 - f_a) + (2 * f_a - 1) * 0.50, 4),
        c_at_b_0p70=round((1 - f_a) + (2 * f_a - 1) * 0.70, 4),
        k_needed=k, c_needed_p05=round(c_need, 4),
        b_break=round(bb, 4) if np.isfinite(bb) else np.nan,
        b_break_reachable=bool(np.isfinite(bb) and bb <= 1.0)))
lev = pd.DataFrame(rows)
print(lev.to_string(index=False))

print("\nreading: b_break is the marginal negative-beta rate that alone makes the endpoint")
print("significant at p<0.05 with ZERO pair-level correspondence. b_break > 1 means the")
print("artefact CANNOT fake a positive in that stratum by marginal shift alone.")

print("\ncluster penalty (DEFF = 1 + (m-1)*rho), m = mean pairs per cluster:")
for name, sub in strata.items():
    if not len(sub):
        continue
    g = sub.sv_id.nunique()
    m = len(sub) / g
    for rho in (0.2, 0.5):
        deff = 1 + (m - 1) * rho
        print("  " + name.ljust(28) + " G=" + str(g).rjust(3)
              + "  m=" + format(m, "5.1f") + "  rho=" + str(rho)
              + "  DEFF=" + format(deff, "6.2f")
              + "  n_eff=" + format(len(sub) / deff, "6.1f"))

tmp = OUT + "/b1_endpoint/b1_endpoint_leverage.tsv.tmp"
lev.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
os.replace(tmp, OUT + "/b1_endpoint/b1_endpoint_leverage.tsv")
print("\nwrote " + OUT + "/b1_endpoint/b1_endpoint_leverage.tsv")
