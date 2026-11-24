import glob, os
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
RD = (CENSV_ROOT + "/discovery"
      "/cdr_phenotypes")
WU = (CENSV_ROOT + "/discovery"
      "/cdr_qtl_full_panel")
ALPHA = 0.05 / (11129 * 8)

def load(d, fr):
    o = []
    for f in sorted(glob.glob(d + "/perchrom/*_perm.tsv")):
        t = pd.read_csv(f, sep="\t")
        if len(t):
            t["frame"] = fr
            o.append(t)
    return pd.concat(o, ignore_index=True) if o else pd.DataFrame()

parts = [load(RD + "/output/scan/scan_1e7/span", "span"),
         load(RD + "/output/scan/scan_1e7/splice", "splice"),
         load(RD + "/output/scan/scan_1e7_density/block", "density")]
E = pd.concat([p for p in parts if len(p)], ignore_index=True)
S = load(WU + "/output/scan/scan_1e7", "shipped")

print("alpha = %.3e   |   escalated rows: corrected %s, shipped %s\n"
      % (ALPHA, format(len(E), ","), format(len(S), ",")))

def summarize(D, name):
    sig = D[D.p_emp_perm < ALPHA]
    g = (sig.groupby(["chrom", "phenotype"])
         .agg(n_sv=("sv_locus", "nunique"), best_p=("p_emp_perm", "min"),
              max_abs_t=("abs_t_obs", "max")).reset_index()
         .sort_values("max_abs_t", ascending=False))
    print("=== %s: (chrom, phenotype) cells clearing study-wide alpha ===" % name)
    print(g.to_string(index=False))
    print(">>> chromosomes: %s\n" % sorted(set(sig.chrom),
          key=lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c)))
    return g

gs = summarize(S, "SHIPPED")
ge = summarize(E, "CORRECTED")

print("=== head-to-head on cells present in both ===")
m = gs.merge(ge, on=["chrom", "phenotype"], suffixes=("_ship", "_corr"), how="outer")
m["status"] = m.apply(lambda r: "both" if pd.notna(r.max_abs_t_ship) and pd.notna(r.max_abs_t_corr)
                      else ("shipped only" if pd.notna(r.max_abs_t_ship) else "corrected only"), axis=1)
print(m[["chrom", "phenotype", "max_abs_t_ship", "max_abs_t_corr", "status"]]
      .sort_values("max_abs_t_ship", ascending=False).to_string(index=False))

os.makedirs(RD + "/output/agg", exist_ok=True)
m.to_csv(RD + "/output/agg/studywide_shipped_vs_corrected.tsv", sep="\t", index=False,
         float_format="%.6g")
E.to_csv(RD + "/output/agg/escalated_1e7_all.tsv.gz", sep="\t", index=False, compression="gzip")
print("\n[written] output/agg/studywide_shipped_vs_corrected.tsv")
