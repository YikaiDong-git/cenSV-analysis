#!/usr/bin/env python3
import hashlib
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
SRC = {
    "report": M + "/reports/chm13_observed_sv_meqtl_mac8_files/all_discovery_conditioning_consensus.tsv",
    "corrected": M + "/results/meqtl_mac8_all_discovery_conditioning_corrected/summary/conditioning_consensus.tsv",
    "uncorrected": M + "/results/meqtl_mac8_all_discovery_conditioning/summary/conditioning_consensus.tsv",
}

d = {}
for k, p in SRC.items():
    df = pd.read_csv(p, sep="\t")
    d[k] = df
    print(k + ": " + str(df.shape) + "  md5="
          + hashlib.md5(open(p, "rb").read()).hexdigest()[:10])
    print("   cols: " + str(len(df.columns)))

cols = {k: list(v.columns) for k, v in d.items()}
print("\ncolumn sets identical across all three: "
      + str(cols["report"] == cols["corrected"] == cols["uncorrected"]))
for k in ["corrected", "uncorrected"]:
    only_a = [c for c in cols["report"] if c not in cols[k]]
    only_b = [c for c in cols[k] if c not in cols["report"]]
    if only_a or only_b:
        print("   report-only vs " + k + ": " + str(only_a))
        print("   " + k + "-only vs report: " + str(only_b))

print("\npair_id sets identical: "
      + str(set(d["report"].pair_id) == set(d["corrected"].pair_id) == set(d["uncorrected"].pair_id)))

base = d["report"].set_index("pair_id").sort_index()
for k in ["corrected", "uncorrected"]:
    o = d[k].set_index("pair_id").sort_index()
    common = [c for c in base.columns if c in o.columns]
    print("\n--- report vs " + k + " ---")
    diffs = []
    for c in common:
        a, b = base[c], o[c]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            same = np.allclose(a.values.astype(float), b.values.astype(float),
                               rtol=0, atol=0, equal_nan=True)
            if not same:
                md = float(np.nanmax(np.abs(a.values.astype(float) - b.values.astype(float))))
                diffs.append((c, "numeric", md, int((a.values != b.values).sum())))
        else:
            neq = int((a.astype(str).values != b.astype(str).values).sum())
            if neq:
                diffs.append((c, "string", np.nan, neq))
    if not diffs:
        print("  IDENTICAL on every shared column -> md5 difference is formatting only")
    else:
        for c, t, md, n in diffs:
            print("  " + c + "  [" + t + "]  n_differing=" + str(n)
                  + ("  max_abs_delta=" + format(md, ".3g") if t == "numeric" else ""))

    if "expected_rint_beta" in common:
        sa = np.sign(base.expected_rint_beta.values)
        sb = np.sign(o.expected_rint_beta.values)
        print("  SIGN of expected_rint_beta differs on " + str(int((sa != sb).sum()))
              + " / " + str(len(sa)) + " pairs   <- this is what the pre-registration rests on")
        print("  frac negative: report " + str(round(float((sa < 0).mean()), 5))
              + "   " + k + " " + str(round(float((sb < 0).mean()), 5)))

r = d["report"]
print("\n=== can the 84-pair set be reconstructed from the shipped columns? ===")
for c in ["conditioning_class_scope", "conditioning_class_pooled"]:
    if c in r:
        print("  " + c + ": " + str(dict(r[c].value_counts(dropna=False))))
for c in ["max_sv_pip", "min_conditional_q_scope", "min_crossfit_q_scope",
          "median_raw_effect_retention", "n_robust_panels_scope"]:
    if c in r:
        v = pd.to_numeric(r[c], errors="coerce")
        print("  " + c + ": n_nonnull=" + str(int(v.notna().sum()))
              + "  min=" + format(float(v.min()), ".4g")
              + "  median=" + format(float(v.median()), ".4g")
              + "  max=" + format(float(v.max()), ".4g"))
if set(["max_sv_pip", "min_conditional_q_scope"]).issubset(r.columns):
    g = ((pd.to_numeric(r.max_sv_pip, errors="coerce") >= 0.9)
         & (pd.to_numeric(r.min_conditional_q_scope, errors="coerce") <= 0.05))
    print("  PIP>=0.9 & conditional q<=0.05: " + str(int(g.sum())) + " pairs  (target: 84)")
    if "median_raw_effect_retention" in r:
        g2 = g & (pd.to_numeric(r.median_raw_effect_retention, errors="coerce") >= 0.5)
        print("  + retention>=0.5             : " + str(int(g2.sum())) + " pairs")
