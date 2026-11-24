#!/usr/bin/env python3
"""A1 step 1: write a truth VCF and an active-interval BED for each 23-chromosome simulation bundle."""
import json, os, sys, glob
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

B = (CENSV_COLLAB_ROOT + "/exploration/simulate_svim_asm_accuracy/"
     "results/multichrom_sv_v1/bundles")
REF = CENSV_ROOT + "/resources/chm13v2.0.fa.gz"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
fa = pysam.FastaFile(REF)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


header = ["##fileformat=VCFv4.2",
          '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
          '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
          '##INFO=<ID=TRUEKIND,Number=1,Type=String,Description="planted operation">',
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
header += ["##contig=<ID=%s,length=%d>" % (c, l) for c, l in zip(fa.references, fa.lengths)]

n_written = 0
for tj in sorted(glob.glob(B + "/*/alpha_mixed/truth.json")):
    d = json.load(open(tj))
    bid = os.path.basename(os.path.dirname(os.path.dirname(tj)))
    recs = []
    for t in d["truth"]:
        c, s, e, k = t["chrom"], int(t["start"]), int(t["end"]), t["kind"]
        anchor = fa.fetch(c, s - 1, s).upper()
        if k == "DEL":
            ref, alt, st, ln = fa.fetch(c, s - 1, e).upper(), anchor, "DEL", -(e - s)
        elif k in ("INS", "TANDEM", "DISPERSED"):
            ins = t["alt"].upper()
            assert len(ins) == t["inserted_length"] and len(ins) > 0, (bid, t["id"])
            ref, alt, st, ln = anchor, anchor + ins, "INS", len(ins)
        elif k == "REPLACEMENT":
            ref, alt = fa.fetch(c, s - 1, e).upper(), anchor + t["alt"].upper()
            st, ln = ("INS" if len(alt) > len(ref) else "DEL"), len(alt) - len(ref)
        elif k == "INV":
            seg = fa.fetch(c, s, e).upper()
            ref, alt, st, ln = anchor + seg, anchor + rc(seg), "INV", e - s
        else:
            raise SystemExit("unknown kind %s" % k)
        assert len(ref) >= 1 and len(alt) >= 1
        recs.append((c, s, t["id"], ref, alt, st, ln, k))
    recs.sort(key=lambda r: (r[0], r[1]))
    with open("%s/%s.truth.vcf" % (OUT, bid), "w") as fh:
        fh.write("\n".join(header) + "\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for c, pos, i, ref, alt, st, ln, k in recs:
            fh.write("%s\t%d\t%s\t%s\t%s\t.\tPASS\tSVTYPE=%s;SVLEN=%d;TRUEKIND=%s\tGT\t1/1\n"
                     % (c, pos, i, ref, alt, st, ln, k))
    with open("%s/%s.active.bed" % (OUT, bid), "w") as fh:
        for a, b in d["active_intervals"]:
            fh.write("%s\t%d\t%d\n" % (d["chrom"], a, b))
    n_written += 1
assert n_written == 69, n_written
print("wrote", n_written, "truth VCFs")
