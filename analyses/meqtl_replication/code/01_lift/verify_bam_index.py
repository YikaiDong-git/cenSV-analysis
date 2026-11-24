#!/usr/bin/env python3
"""Is the BAM index stale?"""
import os
import sys
import pysam
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
BAM_DIR = PROJ + "/alignments/svim_asm"
OUT = UNIT + "/output/a1_lift"

rows = []
for tag in sys.argv[1:]:
    path = BAM_DIR + "/" + tag + "_to_chm13.bam"
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    csi = path + ".csi"
    bai = path + ".bai"
    idx_path = csi if os.path.exists(csi) else bai
    if not os.path.exists(idx_path):
        raise FileNotFoundError("no index for " + path)

    bam_mtime = os.path.getmtime(path)
    idx_mtime = os.path.getmtime(idx_path)

    b = pysam.AlignmentFile(path, "rb")
    from_index = {s.contig: s.mapped for s in b.get_index_statistics()}
    idx_mapped = sum(from_index.values())

    from_file = {}
    n_unmapped = 0
    for a in b.fetch(until_eof=True):
        if a.is_unmapped:
            n_unmapped += 1
            continue
        from_file[a.reference_name] = from_file.get(a.reference_name, 0) + 1
    file_mapped = sum(from_file.values())
    b.close()

    disagree = sorted(set(from_index) | set(from_file))
    bad = [(c, from_index.get(c, 0), from_file.get(c, 0)) for c in disagree
           if from_index.get(c, 0) != from_file.get(c, 0)]

    print("=" * 72)
    print(tag)
    print("  bam mtime newer than index by: "
          + str(round((bam_mtime - idx_mtime) / 86400.0, 1)) + " days")
    print("  mapped records from INDEX: " + str(idx_mapped))
    print("  mapped records from FILE : " + str(file_mapped) + "  (unmapped " + str(n_unmapped) + ")")
    print("  VERDICT: " + ("INDEX CONSISTENT" if not bad else "INDEX STALE -- DO NOT TRUST fetch()"))
    for c, i, f in bad[:10]:
        print("    " + c + ": index " + str(i) + " vs file " + str(f))

    rows.append(dict(bam_tag=tag, idx_mapped=idx_mapped, file_mapped=file_mapped,
                     unmapped=n_unmapped, n_refs_disagreeing=len(bad),
                     consistent=(len(bad) == 0),
                     bam_newer_than_index_days=round((bam_mtime - idx_mtime) / 86400.0, 1)))

D = pd.DataFrame(rows).sort_values("bam_tag")
os.makedirs(OUT, exist_ok=True)
p = OUT + "/bam_index_verification.tsv"
tmp = p + ".tmp." + str(os.getpid())
D.to_csv(tmp, sep="\t", index=False)
os.replace(tmp, p)
print("\nwrote " + p)
if not D.consistent.all():
    raise SystemExit("at least one BAM index is stale -- fetch() results are unreliable")
print("all checked indices are consistent with their BAMs")
