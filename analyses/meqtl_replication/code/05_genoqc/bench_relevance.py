#!/usr/bin/env python3
import json
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
B = M + "/results/observed_sv_genotyping/hprc_truth_benchmark"
R = M + "/reports/chm13_observed_sv_meqtl_mac8_files"
UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")

ps = pd.read_csv(B + "/summary/per_sample_results.tsv", sep="\t")
pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
disc = pd.read_csv(R + "/all_discovery_conditioning_consensus.tsv", sep="\t")

bench_sv = set(ps.candidate_id.unique())
unit_sv = set(pairs.sv_id.unique())
disc_sv = set(disc.sv_id.unique())

print("Q1  benchmarked SVs: " + str(len(bench_sv))
      + "   unit source SVs: " + str(len(unit_sv))
      + "   MAC8 discovery SVs: " + str(len(disc_sv)))
print("    benchmarked INTERSECT unit-47      : " + str(len(bench_sv & unit_sv)))
print("    benchmarked INTERSECT MAC8 discovery: " + str(len(bench_sv & disc_sv)))
for s in sorted(bench_sv):
    print("      " + s + "   in_unit47=" + str(s in unit_sv)
          + "  in_disc=" + str(s in disc_sv))

bm = ps.drop_duplicates("candidate_id")[["candidate_id", "svtype", "svlen"]]
print("\n    benchmarked SV composition: " + str(dict(bm.svtype.value_counts()))
      + "   svlen median " + str(int(bm.svlen.median()))
      + " (min " + str(int(bm.svlen.min())) + ", max " + str(int(bm.svlen.max())) + ")")
up = pairs.drop_duplicates("sv_id")
uplen = (up.source_end - up.source_pos).abs()
print("    unit-47 SV composition    : " + str(dict(up.source_type.value_counts()))
      + "   span median " + str(int(uplen.median()))
      + " (min " + str(int(uplen.min())) + ", max " + str(int(uplen.max())) + ")")

print("\nQ2  recomputed confusion, per tool (truth_dosage as truth)")
rows = []
for tool in ["kanpig_dosage", "sniffles_forced_dosage", "sniffles_discovery_dosage"]:
    d = ps[["truth_dosage", tool]].copy()
    called = d[tool].notna()
    ok = d[called]
    exact = float((ok.truth_dosage == ok[tool]).mean()) if len(ok) else np.nan
    tc = ok[ok.truth_dosage > 0]
    tn = ok[ok.truth_dosage == 0]
    rows.append(dict(
        tool=tool.replace("_dosage", ""), n_rows=int(len(d)),
        call_rate=round(float(called.mean()), 4),
        exact_acc=round(exact, 4),
        n_true_carrier=int(len(tc)),
        carrier_sens=round(float((tc[tool] > 0).mean()), 4) if len(tc) else np.nan,
        n_true_noncarrier=int(len(tn)),
        noncarrier_spec=round(float((tn[tool] == 0).mean()), 4) if len(tn) else np.nan,
        mean_true_dosage=round(float(ok.truth_dosage.mean()), 4),
        mean_called_dosage=round(float(ok[tool].mean()), 4)))
cm = pd.DataFrame(rows)
print(cm.to_string(index=False))

print("\n    truth dosage distribution: " + str(dict(ps.truth_dosage.value_counts().sort_index())))
print("    distinct HPRC samples in benchmark: " + str(ps.sample_id.nunique()))
print("    rows per candidate: " + str(dict(ps.candidate_id.value_counts())))

print("\nQ3  were the reads simulated?")
sj = B + "/smoke/simulation.json"
if os.path.exists(sj):
    print("    " + json.dumps(json.load(open(sj)), indent=2)[:900])
else:
    print("    no simulation.json at " + sj)

print("\n    implied dosage attenuation (var-based), per tool:")
for tool in ["kanpig_dosage", "sniffles_forced_dosage", "sniffles_discovery_dosage"]:
    d = ps[["truth_dosage", tool]].dropna()
    if len(d) > 3 and d[tool].std() > 0:
        r = float(np.corrcoef(d.truth_dosage, d[tool])[0, 1])
        b = float(np.polyfit(d.truth_dosage, d[tool], 1)[0])
        print("      " + tool.replace("_dosage", "")
              + "  corr(truth,called)=" + str(round(r, 4))
              + "  slope called~truth=" + str(round(b, 4))
              + "  n=" + str(len(d)))
