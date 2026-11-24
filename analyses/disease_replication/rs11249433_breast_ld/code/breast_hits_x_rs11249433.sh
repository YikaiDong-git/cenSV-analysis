#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
B=~/AUO_handoff${CENSV_ROOT}/analyses/disease_replication/aou_replication/results_vds
W=${CENSV_ROOT}/analyses/disease_replication/phewide_all_handoff/full_11129_ext2
OUT=${CENSV_ROOT}/analyses/disease_replication/rs11249433_breast_ld/output/breast_hits_x_rs11249433_weights.tsv
python3 - "$B" "$W" "$OUT" <<'PY'
import sys, csv, collections
B, W, OUT = sys.argv[1:4]
POS = "121538815"
hits = []
for run in ["step30_phewide", "step30_phewide_ext2"]:
    with open(f"{B}/{run}/phewide_disease_hits_enriched_flagged.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if "breast" in r["label"].lower():
                r["run"] = run
                hits.append(r)
print("breast hit rows:", len(hits), collections.Counter(h["run"] for h in hits))
print("chroms:", collections.Counter(h["chrom"] for h in hits))
need = {h["sv_locus"] for h in hits}
w = collections.defaultdict(list)
for c in sorted({h["chrom"] for h in hits}):
    with open(f"{W}/sv_weights_long_ext2_{c}.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["sv_locus"] in need:
                w[r["sv_locus"]].append(r)
missing = need - set(w)
print("hit SVs without weights in ext2 bundle:", sorted(missing))
rows = []
for h in hits:
    ws = w[h["sv_locus"]]
    absw = sorted((abs(float(x["beta_std"])) for x in ws), reverse=True)
    tot = sum(absw)
    snp = [x for x in ws if x["chrom"] == "chr1" and x["pos"] == POS]
    if snp:
        b = float(snp[0]["beta_std"])
        rank = absw.index(abs(b)) + 1
        frac = abs(b) / tot
        rows.append([h["run"], h["sv_locus"], h["phecode"], h["p_firth"], h["OR_firth"], h["cv_r2"], h["robust_class"], len(ws), "YES", f"{b:.4f}", rank, f"{frac:.3f}", snp[0]["ref"], snp[0]["alt"], snp[0]["a1_effect"]])
    else:
        mn = min(int(x["pos"]) for x in ws); mx = max(int(x["pos"]) for x in ws)
        rows.append([h["run"], h["sv_locus"], h["phecode"], h["p_firth"], h["OR_firth"], h["cv_r2"], h["robust_class"], len(ws), "no", f"span_{mn}-{mx}", "", "", "", "", ""])
hdr = ["run","sv_locus","phecode","p_firth","OR_firth","cv_r2","robust_class","n_weight_snps","has_rs11249433","beta_std_or_span","abs_rank","frac_abs_weight","ref","alt","a1"]
with open(OUT, "w") as f:
    f.write("\t".join(hdr) + "\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
yes = {r[1] for r in rows if r[8] == "YES"}; allsv = {r[1] for r in rows}
print(f"unique breast-hit SVs: {len(allsv)}; with rs11249433: {len(yes)}")
PY
column -t -s $'\t' "$OUT"
