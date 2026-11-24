#!/usr/bin/env python3
"""Hardy-Weinberg check on the observed genotypes of every tested SV."""
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
TRA = f"{M}/results/meqtl_trans/chm13_observed_sv_all_family_pc10_mac8"
U = (CENSV_ROOT + "/analyses/"
     "meqtl_replication")
OUT = f"{U}/output/arms/refit_genotype_complete"

z = np.load(f"{TRA}/manifests/profile_dosages.npz", allow_pickle=True)
D, prof_ids = z["dosages"], list(z["profile_ids"])
fam = pd.read_csv(f"{TRA}/manifests/family_to_profile.tsv", sep="\t")
p2f = dict(zip(fam.profile_id, fam.physical_family_id))

R = pd.read_csv(f"{OUT}/r4_refit.tsv", sep="\t", low_memory=False)
R = R[R.status == "ok"]
arms = R.groupby("family_id").arm.agg(
    lambda s: "distal_or_trans" if set(s) & {"cis_distal", "trans"} else "local_only")

rows = []
for i, pid in enumerate(prof_ids):
    g = D[i]
    obs = g[np.isin(g, [0.0, 1.0, 2.0])]
    n = len(obs)
    n0, n1, n2 = int((obs == 0).sum()), int((obs == 1).sum()), int((obs == 2).sum())
    if n == 0:
        continue
    p = (n1 + 2 * n2) / (2 * n)
    e0, e1, e2 = n * (1 - p) ** 2, n * 2 * p * (1 - p), n * p ** 2
    chi = sum((o - e) ** 2 / e for o, e in ((n0, e0), (n1, e1), (n2, e2)) if e > 0)
    pv = float(stats.chi2.sf(chi, 1))
    fam_id = p2f.get(pid, pid)
    rows.append(dict(profile_id=pid, family_id=fam_id, n_observed=n,
                     n_hom_ref=n0, n_het=n1, n_hom_alt=n2,
                     af=p, exp_het=e1, exp_hom_alt=e2, chi2=chi, hwe_p=pv,
                     arm=arms.get(fam_id, "not_a_discovery")))
H = pd.DataFrame(rows).sort_values("hwe_p")
H.to_csv(f"{OUT}/r6_hwe.tsv", sep="\t", index=False, float_format="%.6g")

print("=== the SV that carries the distal/trans layer ===")
t = H[H.family_id.str.startswith("chr17_26788109")]
print(t[["family_id", "n_observed", "n_hom_ref", "n_het", "n_hom_alt", "af",
         "exp_het", "exp_hom_alt", "hwe_p"]].to_string(index=False))

print("\n=== HWE by discovery arm (SVs that produced >=1 discovery) ===")
d = H[H.arm != "not_a_discovery"]
for arm, g in d.groupby("arm"):
    print(f"  {arm:<18} n_SVs={len(g):3d}  median HWE p={g.hwe_p.median():.3g}  "
          f"HWE p<0.05: {int((g.hwe_p < 0.05).sum())} ({100*(g.hwe_p<0.05).mean():.0f}%)  "
          f"median hom_alt excess={(g.n_hom_alt - g.exp_hom_alt).median():+.2f}")

print("\n=== the ten strongest HWE departures among discovery SVs ===")
print(d.nsmallest(10, "hwe_p")[["family_id", "arm", "n_observed", "n_hom_ref", "n_het",
                                "n_hom_alt", "exp_het", "exp_hom_alt", "hwe_p"]]
      .to_string(index=False, float_format=lambda v: f"{v:.3g}"))
print(f"\nwrote {OUT}/r6_hwe.tsv")
