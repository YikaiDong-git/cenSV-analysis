#!/usr/bin/env python3
import sys, glob, csv
import numpy as np
from collections import defaultdict

diag_glob, out_land, out_cat = sys.argv[1], sys.argv[2], sys.argv[3]
files = sorted(glob.glob(diag_glob))

COLS = ["region", "trait", "arm", "mode", "sv_locus", "sv_index", "cv_r2", "sv_pip",
        "svSNP_r", "z_SV", "in_cs", "p_snps", "m_sv", "s_rss", "converged", "niter", "wall_s"]
snpr_hist = np.zeros(20)
cv_acc = defaultdict(lambda: [0, 0.0])
n_hb = n_distinct = n_tagged = n_conv = 0
wall_sum = 0.0
cand = defaultdict(dict)
CAND_CAP = 400000

for f in files:
    try:
        fh = open(f)
    except OSError:
        continue
    r = csv.reader(fh, delimiter="\t")
    try:
        hdr = next(r)
    except StopIteration:
        fh.close(); continue
    ix = {c: i for i, c in enumerate(hdr)}
    for row in r:
        if len(row) < len(COLS) or row[ix["mode"]] != "hitbundle":
            continue
        try:
            snpr = float(row[ix["svSNP_r"]]); pip = float(row[ix["sv_pip"]])
            z = float(row[ix["z_SV"]]); cv = float(row[ix["cv_r2"]])
        except ValueError:
            continue
        n_hb += 1
        if row[ix["converged"]] in ("TRUE", "True"):
            n_conv += 1
        try:
            wall_sum += float(row[ix["wall_s"]])
        except ValueError:
            pass
        snpr_hist[min(19, int(snpr * 20))] += 1
        if snpr < 0.85:
            n_distinct += 1
        if snpr >= 0.95:
            n_tagged += 1
        cb = round(cv, 1); cv_acc[cb][0] += 1; cv_acc[cb][1] += snpr
        if (pip > 0.5 or abs(z) > 5) and len(cand) < CAND_CAP:
            cand[(row[ix["sv_locus"]], row[ix["trait"]])][row[ix["arm"]]] = (
                pip, snpr, z, row[ix["in_cs"]], row[ix["s_rss"]], cv, row[ix["region"]])
    fh.close()

with open(out_cat, "w") as fh:
    fh.write("sv_locus\tregion\ttrait\tcv_r2\tpip_matched\tpip_impute\tsvSNP_r\tz_SV\t"
             "in_cs_m\tin_cs_i\ts_rss\tdivergence\tclass\n")
    ncat = 0
    for (loc, trait), arms in cand.items():
        m = arms.get("matched"); i = arms.get("impute")
        pm = m[0] if m else float("nan"); pi = i[0] if i else float("nan")
        ref = i or m
        snpr, z, cv, region = ref[1], ref[2], ref[5], ref[6]
        s = (i[4] if i else m[4]) if (i or m) else "NA"
        incs_m = m[3] if m else "NA"; incs_i = i[3] if i else "NA"
        div = abs((pm if pm == pm else 0) - (pi if pi == pi else 0)) > 0.4
        distinct = snpr < 0.85
        credible = (pm == pm and pm > 0.5) or (pi == pi and pi > 0.5)
        cls = ("DISTINCT_CREDIBLE" if (distinct and credible) else
               "TAGGED_CREDIBLE" if credible else "SUBTHRESHOLD")
        if not credible and abs(z) <= 5:
            continue
        fh.write(f"{loc}\t{region}\t{trait}\t{cv:.4f}\t{pm:.4f}\t{pi:.4f}\t{snpr:.4f}\t{z:.4f}\t"
                 f"{incs_m}\t{incs_i}\t{s}\t{int(div)}\t{cls}\n")
        ncat += 1

d = []
d.append(f"# FULL genome-wide single-SV landscape  (hit_bundle rows={n_hb:,}, files={len(files)})")
d.append(f"converged={100*n_conv/max(1,n_hb):.1f}%  mean_wall={wall_sum/max(1,n_hb):.2f}s  total_fit_core_h={wall_sum/3600:.0f}")
d.append(f"identifiability: distinct(svSNP_r<0.85)={100*n_distinct/max(1,n_hb):.1f}%  tagged(>=0.95)={100*n_tagged/max(1,n_hb):.1f}%")
d.append("\nsvSNP_r histogram (bin=0.05):")
for b in range(20):
    d.append(f"  [{b*0.05:.2f}-{(b+1)*0.05:.2f}) {'#'*int(60*snpr_hist[b]/max(1,snpr_hist.max()))} {int(snpr_hist[b]):,}")
d.append("\nmean svSNP_r by cv_r2 bin (the cv x tagging confound):")
for cb in sorted(cv_acc):
    n, s = cv_acc[cb]
    d.append(f"  cv~{cb:.1f}: n={n:,} mean_svSNP_r={s/max(1,n):.3f}")
d.append(f"\ncatalog cells (credible or |z_SV|>5): {ncat} -> {out_cat}")
open(out_land, "w").write("\n".join(d) + "\n")
print("\n".join(d))
print(f"\n[wrote] {out_land} + {out_cat} ({ncat} catalog cells)")
