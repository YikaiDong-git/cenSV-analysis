#!/usr/bin/env python3
"""chr7_sv114 -> stroke lookup (discovery-side follow-up to the AoU cerebrovascular hit)."""
import sys

WEIGHT, STR_CHR7 = sys.argv[1], sys.argv[2]
REG_S, REG_E = 57314931, 63846981
GWS = 5e-8

tags = {}
with open(WEIGHT) as fh:
    for i, line in enumerate(fh):
        if i == 0:
            continue
        f = line.rstrip("\n").split("\t")
        snp_id, w = f[9], float(f[10])
        if abs(w) <= 1e-8:
            continue
        p = snp_id.split(":")
        tags[int(p[1])] = (p[2], p[3], w)

region_max_absz, region_min_p, region_gws, n_region = 0.0, 1.0, [], 0
tag_rows = []
with open(STR_CHR7) as fh:
    for i, line in enumerate(fh):
        if i == 0:
            continue
        f = line.rstrip("\n").split("\t")
        sp = f[0].split("_")
        if len(sp) < 4:
            continue
        pos = int(sp[1])
        if pos < REG_S or pos > REG_E:
            continue
        n_region += 1
        beta, se, pval = float(f[3]), float(f[4]), float(f[5])
        z = beta / se if se > 0 else 0.0
        if abs(z) > region_max_absz:
            region_max_absz = abs(z)
        if pval < region_min_p:
            region_min_p = pval
        if pval < GWS:
            region_gws.append((pos, f[0], z, pval))
        if pos in tags:
            ref, alt, w = tags[pos]
            tag_rows.append((pos, f[0], w, z, pval))

print("# chr7_sv114 -> FinnGen I9_STR (stroke excl SAH; 34110 ca / 450023 co) lookup")
print(f"# region chr7:{REG_S}-{REG_E}  region_SNPs={n_region}  tag_SNPs_matched={len(tag_rows)}/{len(tags)}")
print(f"# REGION-WIDE   max|z|={region_max_absz:.3f}   min_p={region_min_p:.3e}   n_GWS(p<5e-8)={len(region_gws)}")
if region_gws:
    print("#   top region GWS stroke SNPs:")
    for pos, snp, z, p in sorted(region_gws, key=lambda x: x[3])[:10]:
        print(f"#     {snp}  z={z:+.2f}  p={p:.2e}")
tag_rows.sort(key=lambda x: -abs(x[3]))
wsum = sum(abs(r[2]) for r in tag_rows) or 1.0
wz = sum(abs(r[2]) * r[3] for r in tag_rows)
print(f"# TAG-SNP stroke stats ({len(tags)} predictor SNPs; sorted by |z|):")
for pos, snp, w, z, p in tag_rows:
    flag = " <== GWS" if p < GWS else (" <- nominal" if p < 5e-3 else "")
    print(f"#   pos={pos} w={w:+.4g} stroke_z={z:+.3f} p={p:.2e}{flag}")
if tag_rows:
    print(f"# TAG-SNP  max|z|={max(abs(r[3]) for r in tag_rows):.3f}  min_p={min(r[4] for r in tag_rows):.3e}"
          f"  |w|-wtd_mean_z={wz/wsum:+.4f}")
print("# INTERPRETATION: GWS in region but not at tag SNPs => signal the SV under-tagged (weak predictor);"
      " GWS AT tag SNPs => discovery under-powered a real SV-tagged stroke signal; no GWS anywhere => locus"
      " likely carries no stroke signal at FinnGen power (AoU hit stays exploratory, needs bigger stroke GWAS).")
