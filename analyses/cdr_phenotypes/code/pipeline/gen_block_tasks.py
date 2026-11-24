#!/usr/bin/env python
"""Task file for BLK/T2: finish the BLOCK-frame escalation."""
import glob
import os
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
UNIT = ROOT + "/analyses/cdr_phenotypes"
OUT = UNIT + "/output/tasks/escalate_block_completion.txt"

files = sorted(glob.glob(UNIT + "/output/scan/scan_1e5/block/perchrom/*_perm.tsv"))
assert files, "no block 1e5 scan files found"
df = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
floor = df.p_emp_perm.min()
sat = set(df[df.p_emp_perm <= floor * 1.0001]
          .groupby(["chrom", "phenotype"]).size().index)
print("block 1e5: %d tests, floor p = %.3e, %d floor-saturated cells"
      % (len(df), floor, len(sat)))

done = set()
for pat in ("/output/scan/scan_1e7_density/block/perchrom/*_perm.tsv",
            "/output/scan/scan_1e7/block/perchrom/*_perm.tsv"):
    for f in glob.glob(UNIT + pat):
        b = os.path.basename(f).replace("_perm.tsv", "")
        c, p = b.split("_", 1)
        done.add((c, p))
print("already escalated in block: %d" % len(done))

missing = sorted(sat - done)
print("to escalate: %d" % len(missing))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    for c, p in missing:
        fh.write("%s\t%s\tblock\n" % (c, p))
        print("   %-6s %s" % (c, p))
print("[written] %s" % OUT)
print("submit with: TASKFILE=%s NPERM=10000000 OUTBASE=%s/output/scan/scan_1e7 "
      "sbatch --array=1-%d%%6 %s/code/pipeline/run_rescan_array.sh"
      % (OUT, UNIT, len(missing), UNIT))
