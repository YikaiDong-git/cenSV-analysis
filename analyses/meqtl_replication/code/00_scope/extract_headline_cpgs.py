#!/usr/bin/env python3
import glob
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
OUT = UNIT + "/output/headline_svs"
os.makedirs(OUT, exist_ok=True)

TARGETS = [("chr12_DEL333", "chr12", 37183837), ("chr8_DEL113", "chr8", 44300221)]

cands = []
for pat in ["/results/meqtl_cis/*/summary/**/*assoc*",
            "/results/meqtl_cis/*/summary/*.tsv.gz",
            "/results/meqtl_cis/*/summary/*.tsv"]:
    cands += glob.glob(M + pat, recursive=True)
cands = sorted(set(c for c in cands if os.path.getsize(c) > 5000))
print("candidate result tables: " + str(len(cands)))

rows = []
for path in cands:
    try:
        d = pd.read_csv(path, sep="\t", nrows=200000)
    except Exception:
        continue
    cc = "sv_chrom" if "sv_chrom" in d.columns else ("source_chrom" if "source_chrom" in d.columns else None)
    pc = "sv_pos" if "sv_pos" in d.columns else ("source_pos" if "source_pos" in d.columns else None)
    if cc is None or pc is None or "feature_id" not in d.columns:
        continue
    rel = path.replace(M + "/", "")
    for name, ch, pos in TARGETS:
        hit = d[(d[cc] == ch) & ((d[pc] - pos).abs() <= 500)]
        if not len(hit):
            continue
        qcols = [c for c in hit.columns if "_q" in c.lower() or c.lower().endswith("q")]
        qc = None
        for pref in ["rint_hc3_bh_q_primary_t", "bh_q", "q_value", "qval", "discovery_q"]:
            m = [c for c in qcols if pref in c]
            if m:
                qc = m[0]
                break
        nsig = int((pd.to_numeric(hit[qc], errors="coerce") < 0.01).sum()) if qc else -1
        rows.append(dict(sv=name, table=rel, n_rows=len(hit), q_col=qc or "(none)",
                         n_q_lt_01=nsig))
S = pd.DataFrame(rows)
if len(S):
    print("\n=== where each headline SV appears ===")
    print(S.sort_values(["sv", "n_q_lt_01"], ascending=[True, False]).to_string(index=False))
else:
    print("no table contains either SV with a feature_id column")

for name, ch, pos in TARGETS:
    s = S[S.sv == name]
    if not len(s):
        print("\n" + name + ": no source table found")
        continue
    best = s.sort_values(["n_q_lt_01", "n_rows"], ascending=False).iloc[0]
    d = pd.read_csv(M + "/" + best.table, sep="\t")
    cc = "sv_chrom" if "sv_chrom" in d.columns else "source_chrom"
    pc = "sv_pos" if "sv_pos" in d.columns else "source_pos"
    hit = d[(d[cc] == ch) & ((d[pc] - pos).abs() <= 500)].copy()
    qc = best.q_col if best.q_col != "(none)" else None
    if qc:
        hit["q"] = pd.to_numeric(hit[qc], errors="coerce")
        sig = hit[hit.q < 0.01]
    else:
        hit["q"] = np.nan
        sig = hit
    bcol = "rint_beta" if "rint_beta" in hit.columns else (
        "expected_rint_beta" if "expected_rint_beta" in hit.columns else None)
    keep = ["feature_id", "q"] + ([bcol] if bcol else [])
    out = (sig if len(sig) else hit)[keep].drop_duplicates("feature_id")
    out["sv_name"] = name
    out["sv_chrom"] = ch
    out["sv_pos"] = pos
    out["source_table"] = best.table
    tmp = OUT + "/" + name + "_cpgs.tsv.tmp"
    out.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, OUT + "/inputs/" + name + "_cpgs.tsv")
    print("\n" + name + ": " + str(len(out)) + " CpGs from " + best.table
          + ("  (q<0.01)" if len(sig) else "  (no q filter applied)"))
    if bcol:
        print("   frac negative beta: " + format(float((out[bcol] < 0).mean()), ".4f"))
    print("   manuscript lead present: "
          + str(any(out.feature_id.astype(str).str.startswith("chr12:131351745"))
                if name == "chr12_DEL333" else "n/a"))
