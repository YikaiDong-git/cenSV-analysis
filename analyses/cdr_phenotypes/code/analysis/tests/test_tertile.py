#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""M8 shows short and long arrays placing dips differently in RELATIVE coordinates."""
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
UNIT = ROOT + "/analyses/cdr_phenotypes"
ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

dp = pd.read_csv(UNIT + "/output/phenotypes/per_dip_positions.tsv", sep="\t")
dp = dp[dp.chm13_chrom.isin(ORDER)].copy()
dp["mb"] = dp.active_len_sum / 1e6
dp["rel"] = (dp.pos_splice + 1.0) / 2.0
dp = dp[np.isfinite(dp.rel) & np.isfinite(dp.mb)]
print("dips: %d" % len(dp))

q1, q2 = dp.mb.quantile([1 / 3, 2 / 3])
sh = dp[dp.mb < q1].rel.values
lo = dp[dp.mb >= q2].rel.values
D, p = stats.ks_2samp(sh, lo)
print("")
print("1. POOLED across chromosomes")
print("   short (<%.2f Mb) n=%d   long (>=%.2f Mb) n=%d" % (q1, len(sh), q2, len(lo)))
print("   KS D = %.4f   P = %.3g" % (D, p))
print("   median rel: short %.3f  long %.3f" % (np.median(sh), np.median(lo)))
print("   -> pooled test CONFOUNDED by chromosome: acrocentrics are both short and start-loaded.")

print("")
print("2. WITHIN CHROMOSOME (the confounder removed)")
print("   %-7s %6s %6s %9s %9s %10s" % ("chrom", "n_sh", "n_lo", "med short", "med long", "KS P"))
rows = []
for c in ORDER:
    d = dp[dp.chm13_chrom == c]
    if len(d) < 150:
        continue
    a, b = d.mb.quantile([1 / 3, 2 / 3])
    s = d[d.mb < a].rel.values
    l = d[d.mb >= b].rel.values
    if len(s) < 40 or len(l) < 40:
        continue
    D2, p2 = stats.ks_2samp(s, l)
    rows.append(dict(chrom=c, n_sh=len(s), n_lo=len(l), ms=np.median(s), ml=np.median(l),
                     D=D2, p=p2))
    print("   %-7s %6d %6d %9.3f %9.3f %10.2e" % (c, len(s), len(l), np.median(s),
                                                  np.median(l), p2))
R = pd.DataFrame(rows)
nsig = int((R.p < 0.05 / max(len(R), 1)).sum())
print("")
print("   %d of %d chromosomes differ at Bonferroni 0.05/%d" % (nsig, len(R), len(R)))
print("   median shift (long - short) = %+.3f" % np.median(R.ml - R.ms))
if nsig >= max(3, len(R) // 3):
    print("   -> the difference SURVIVES within chromosome: relative dip position is not")
    print("      distributionally independent of array length.")
else:
    print("   -> the pooled difference does NOT survive within chromosome: it was driven by")
    print("      chromosomes differing in both length and dip placement.")

print("")
print("3. SHAPE of the pooled relative-position distribution")
h, e = np.histogram(dp.rel, bins=40, range=(0, 1))
cent = 0.5 * (e[1:] + e[:-1])
print("   peak bin at rel = %.3f  (share %.3f)" % (cent[h.argmax()], h.max() / h.sum()))
sm = np.convolve(h, np.ones(3) / 3.0, mode="same")
peaks = [i for i in range(1, len(sm) - 1) if sm[i] > sm[i - 1] and sm[i] > sm[i + 1]]
peaks = sorted(peaks, key=lambda i: -sm[i])[:3]
print("   local maxima at rel = %s" % ", ".join("%.2f" % cent[i] for i in sorted(peaks)))
if len(peaks) >= 2:
    a, b = sorted(peaks[:2])
    trough = sm[a:b + 1].min()
    print("   trough between the top two = %.1f%% of the lower peak"
          % (100 * trough / min(sm[a], sm[b])))
print("   share of dips in the first 20%% of the array: %.3f" % (dp.rel < 0.2).mean())
print("   share in the last 20%%:                      %.3f" % (dp.rel > 0.8).mean())
