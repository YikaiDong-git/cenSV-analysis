#!/usr/bin/env python3
import json
import os
import subprocess
import sys

import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_LAB = _censv_os.environ["CENSV_LAB"]

ROOT = CENSV_ROOT + ""
PIPE = os.path.join(ROOT, "analyses/ld_maps/code")
FRAMES = os.path.join(ROOT, "analyses/ld_maps/output/_censv_frames")
KG = CENSV_LAB + "/shared/1KG"
sys.path.insert(0, PIPE)
import plot_sv_ld_pubfig as P                                              # noqa: E402

chrom, side = sys.argv[1], sys.argv[2]
s = json.load(open(os.path.join(FRAMES, chrom, "summary.json")))
samples = s["samples"]
region = s["flanks"][side]["region"]
vcf = os.path.join(KG, "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_%s."
                        "recalibrated_variants.vcf.gz" % chrom)

gt, var, keep = P.load_snp_panel_from_vcf(vcf_path=__import__("pathlib").Path(vcf),
                                          region=region, target_samples=samples,
                                          min_maf=0.01)
loader = list(var.pos.astype(int))
print("loader: %d rows, %d distinct positions, %d samples kept"
      % (len(loader), len(set(loader)), len(keep)))
print("shipped summary says n_snp = %d" % s["flanks"][side]["n_snp"])

sf = "/tmp/_diff_%s_%s.samples" % (chrom, side)
open(sf, "w").write("\n".join(samples) + "\n")
cmd = ["bcftools", "query", "-r", region, "-S", sf,
       "-f", "%POS\t%REF\t%ALT[\t%GT]\n", vcf]
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
mine, nfields = [], set()
for line in proc.stdout:
    f = line.rstrip("\n").split("\t")
    nfields.add(len(f))
    if len(f) < 4 or len(f[1]) != 1 or len(f[2]) != 1 or "," in f[2]:
        continue
    d = np.array([P.gt_to_dosage(g) for g in f[3:]], dtype=float)
    if np.all(np.isnan(d)):
        continue
    af = float(np.nanmean(d) / 2.0)
    if min(af, 1 - af) < 0.01:
        continue
    mine.append(int(f[0]))
proc.stdout.close(); proc.wait()
os.remove(sf)
print("scan  : %d rows, %d distinct positions; field counts seen %s"
      % (len(mine), len(set(mine)), sorted(nfields)))

L, M = set(loader), set(mine)
print("only in scan  : %d  %s" % (len(M - L), sorted(M - L)[:8]))
print("only in loader: %d  %s" % (len(L - M), sorted(L - M)[:8]))
from collections import Counter
cl, cm = Counter(loader), Counter(mine)
dup_l = {p: c for p, c in cl.items() if c > 1}
dup_m = {p: c for p, c in cm.items() if c > 1}
print("duplicate positions - loader %d, scan %d" % (len(dup_l), len(dup_m)))
extra = sorted(M - L)[:3]
if extra:
    print("\nthe first records the scan keeps and the loader drops:")
    sf2 = "/tmp/_diff2.samples"
    open(sf2, "w").write("\n".join(samples) + "\n")
    for p in extra:
        out = subprocess.run(["bcftools", "query", "-r", "%s:%d-%d" % (chrom, p, p),
                              "-S", sf2, "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\n", vcf],
                             capture_output=True, text=True).stdout.strip()
        print("  " + out.replace("\n", "\n  "))
    os.remove(sf2)
