import pandas as pd, numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
RD = (CENSV_ROOT + "/discovery"
      "/cdr_phenotypes")
V1D = (CENSV_ROOT + "/discovery"
       "/cdr_qtl/output")
v2 = pd.read_csv(RD + "/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv", sep="\t", engine="python")
v1 = pd.read_csv(V1D + "/per_hap_chrom_phenotypes.tsv", sep="\t", engine="python")
F = pd.read_csv(RD + "/output/phenotypes/frame_comparison.tsv", sep="\t")

def ordk(c):
    x = str(c).replace("chr", "")
    return (0, int(x)) if x.isdigit() else (1, x)

rows = []
for c in sorted(v2.chm13_chrom.unique(), key=ordk):
    s = v2[v2.chm13_chrom == c]
    mb = s[s.n_blocks > 1]
    d2 = len(set(v1.loc[v1.chm13_chrom == c, "sample"]) - set(s["sample"]))
    fr = F[F.chrom == c]
    rows.append(dict(
        chrom=c, n=len(s),
        pct_multiblock=100.0 * (s.n_blocks > 1).mean(),
        med_gap_kb_when_split=(mb.gap_bp.median() / 1e3) if len(mb) else 0.0,
        len_inflation=s.active_len_span.median() / s.active_len_sum.median(),
        max_inflation=(s.active_len_span / s.active_len_sum).max(),
        r_span_vs_block=float(fr.r_span_vs_block.iloc[0]) if len(fr) else np.nan,
        D2_rows_dropped=d2))
R = pd.DataFrame(rows)

def verdict(r):
    tags = []
    if r.len_inflation > 1.10: tags.append("LENGTH")
    elif r.len_inflation > 1.02: tags.append("length-minor")
    if r.r_span_vs_block < 0.95: tags.append("FRAME")
    elif r.r_span_vs_block < 0.995: tags.append("frame-minor")
    if r.D2_rows_dropped: tags.append("CONTIG")
    return "+".join(tags) if tags else "clean"

R["affected_by"] = R.apply(verdict, axis=1)
pd.set_option("display.width", 200)
print(R.to_string(index=False, float_format=lambda v: "%.3f" % v))
print("\n--- summary ---")
print(R.affected_by.value_counts().to_string())
print("\nchromosomes materially affected: %s"
      % sorted(R.loc[R.affected_by.str.contains("LENGTH|FRAME|CONTIG"), "chrom"], key=ordk))
print("chromosomes fully clean       : %s"
      % sorted(R.loc[R.affected_by == "clean", "chrom"], key=ordk))
c1 = v2[(v2.chm13_chrom == "chr1") & (v2.n_blocks > 1)]
print("\nchr1 probe: %d multi-block of %d; gap kb median %.1f, p90 %.1f; "
      "median SUM %.2f Mb" % (len(c1), (v2.chm13_chrom == "chr1").sum(),
      c1.gap_bp.median()/1e3, c1.gap_bp.quantile(.9)/1e3,
      v2.loc[v2.chm13_chrom=="chr1","active_len_sum"].median()/1e6))
