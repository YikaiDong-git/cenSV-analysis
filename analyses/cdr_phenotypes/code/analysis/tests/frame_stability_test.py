import pandas as pd, numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
RD = (CENSV_ROOT + "/discovery"
      "/cdr_phenotypes")
v2 = pd.read_csv(RD + "/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv", sep="\t", engine="python")
mb = v2[(v2.n_blocks > 1) & v2.cdr_pos_offset_span.notna()].copy()
mb["gap_frac"] = mb.gap_bp / mb.active_len_span
mb["d_block"] = (mb.cdr_pos_offset_block - mb.cdr_pos_offset_span).abs()
mb["d_splice"] = (mb.cdr_pos_offset_splice - mb.cdr_pos_offset_span).abs()

print("Multi-block haplotype-chromosomes: %d" % len(mb))
print("\nDoes each frame react IN PROPORTION to how much gap there actually is?")
print("(a well-behaved frame should barely move when the gap is negligible)\n")
bins = [(0, .001), (.001, .01), (.01, .05), (.05, .2), (.2, 1.01)]
print("%-16s %6s %12s %14s %14s" % ("gap/hull", "n", "med gap kb", "|d| BLOCK", "|d| SPLICE"))
for lo, hi in bins:
    s = mb[(mb.gap_frac >= lo) & (mb.gap_frac < hi)]
    if not len(s):
        continue
    print("%-16s %6d %12.1f %14.4f %14.4f"
          % ("%.3f-%.3f" % (lo, hi), len(s), s.gap_bp.median()/1e3,
             s.d_block.median(), s.d_splice.median()))

triv = mb[mb.gap_frac < 0.001]
print("\nTRIVIAL splits (gap < 0.1%% of the hull; median %.2f kb of a %.2f Mb array): n=%d"
      % (triv.gap_bp.median()/1e3, triv.active_len_span.median()/1e6, len(triv)))
print("   BLOCK frame moves the phenotype by a median of %.4f  (max %.4f)"
      % (triv.d_block.median(), triv.d_block.max()))
print("   SPLICE frame moves it by                        %.4f  (max %.4f)"
      % (triv.d_splice.median(), triv.d_splice.max()))
print("\n-> a gap of ~%.0f bp in a multi-megabase array is an annotation break, not biology."
      % triv.gap_bp.median())
print("   Any frame that reacts strongly to it is measuring the annotation, not the centromere.")
