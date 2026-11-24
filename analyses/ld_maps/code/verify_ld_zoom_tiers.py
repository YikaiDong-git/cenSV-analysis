#!/usr/bin/env python3
"""The tier measurement must reproduce the shipped frames at its narrowest tier."""
import glob
import json
import os
import sys
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

FRAMES = (CENSV_ROOT + "/analyses/"
          "ld_maps/output/_censv_frames")
TIERS = (CENSV_ROOT + "/analyses/"
         "ld_maps/output/_ld_zoom_tiers")

files = sorted(glob.glob(os.path.join(TIERS, "*_tiers.json")))
if not files:
    sys.exit("no tier files")

bad, n_flank, ship_tot, tier_tot = [], 0, 0, 0
for f in files:
    t = json.load(open(f))
    s = json.load(open(os.path.join(FRAMES, t["chrom"], "summary.json")))
    if t["n_samples"] != s["n_samples"]:
        bad.append("%s: %d samples measured, %d shipped"
                   % (t["chrom"], t["n_samples"], s["n_samples"]))
    for side, fl in t["flanks"].items():
        n_flank += 1
        got = fl["per_tier"]["100000"]["n_snp"]
        want = s["flanks"][side]["n_snp"]
        ship_tot += want
        tier_tot += got
        if got != want:
            bad.append("%s %s: tier 100 kb holds %d, the shipped block holds %d"
                       % (t["chrom"], side, got, want))
        if fl["shipped_region"] != s["flanks"][side]["region"]:
            bad.append("%s %s: measured against region %s, shipped is %s"
                       % (t["chrom"], side, fl["shipped_region"],
                          s["flanks"][side]["region"]))

print("%d chromosomes, %d flanks" % (len(files), n_flank))
print("tier 100 kb total %d   shipped total %d" % (tier_tot, ship_tot))
tb = json.load(open(files[0]))["tiers_bp"]
tot = {str(b): 0 for b in tb}
for f in files:
    for fl in json.load(open(f))["flanks"].values():
        for b in tb:
            tot[str(b)] += fl["per_tier"][str(b)]["n_snp"]
for b in tb:
    print("  %-10s %8d SNPs   %5.1fx the shipped set"
          % (("%.0f kb" % (b / 1e3)) if b < 1e6 else ("%.0f Mb" % (b / 1e6)),
             tot[str(b)], tot[str(b)] / ship_tot))
if bad:
    print("\nFAIL (%d)" % len(bad))
    for b in bad[:20]:
        print("  " + b)
    sys.exit(1)
print("\nPASS: every flank's 100 kb tier equals the block the browser ships")
