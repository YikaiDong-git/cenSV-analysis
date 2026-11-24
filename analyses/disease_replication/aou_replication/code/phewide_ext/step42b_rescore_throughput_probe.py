#!/usr/bin/env python3
"""Step 42b — price a full 11,129 rescore with a FIRST-TOUCH throughput measurement."""
import os, sys, csv, glob, time, json
from collections import defaultdict
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step42_gate2_headline_rescore as G

WARM = {"chr11", "chr16"}
PROBE_CHROMS = ["chr1", "chr9", "chr17", "chr20"]


def main():
    per_chrom = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    fam_iids, cols, sv_fam, _ = G.build_cohort()
    Nfam = len(fam_iids)
    bpv = (Nfam + 3) // 4
    print(f"Nfam={Nfam:,}  bytes/variant-row={bpv:,}", flush=True)

    picks = defaultdict(list)
    for ch in PROBE_CHROMS:
        p = f"{G.EXT2}/sv_weights_long_ext2_{ch}.tsv"
        if not os.path.exists(p):
            continue
        rows = defaultdict(list)
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows[r["sv_locus"].strip()].append({k: v.strip() for k, v in r.items()})
        for sv in sorted(rows)[:per_chrom]:
            picks[ch].append(sv)

    results = []
    for ch in PROBE_CHROMS:
        if ch in WARM or not picks[ch]:
            continue
        bedp = f"{G.BED}/{ch}/{ch}.pericentromeric_snps.min_ac1000.bed"
        if not os.path.exists(bedp):
            print(f"  {ch}: no .bed -> skip", flush=True); continue
        want = set(picks[ch])
        W = G.load_weight_rows([f"{G.EXT2}/sv_weights_long_ext2_{ch}.tsv"], want)
        nrows = sum(len(v) for v in W.values())
        t0 = time.time()
        out = G.score_svs(W, fam_iids, cols, f"probe-{ch}")
        dt = time.time() - t0
        nsv = len(out)
        used = sum(v[2] for v in out.values())
        results.append(dict(chrom=ch, n_svs=nsv, n_weight_rows=nrows, n_snps_read=used,
                            secs=dt, secs_per_sv=dt / max(nsv, 1),
                            secs_per_weight_row=dt / max(nrows, 1),
                            MB_per_s=(used * bpv / 1e6) / dt if dt > 0 else float("nan")))
        print(f"  {ch}: {nsv} SVs / {nrows} weight rows / {used} reads in {dt:.1f}s "
              f"-> {dt/max(nrows,1)*1000:.1f} ms per weight row", flush=True)

    if not results:
        raise SystemExit("no probe chromosomes were readable")

    tot_rows = sum(r["n_weight_rows"] for r in results)
    tot_secs = sum(r["secs"] for r in results)
    per_row = tot_secs / tot_rows
    EXT2_ROWS = 328468
    est_h = per_row * EXT2_ROWS / 3600
    print(f"\n  FIRST-TOUCH: {tot_rows:,} weight rows in {tot_secs:.1f}s = {per_row*1000:.1f} ms/row",
          flush=True)
    print(f"  => full ext2 rescore ({EXT2_ROWS:,} rows, 11,129 SVs) ~= {est_h:.2f} h single-process; "
          f"~{est_h/6:.2f} h at 6-way chromosome parallelism", flush=True)
    print("  NOTE: first-touch, not cold-cache; excludes the sv.dose write, the step28b pgen import, "
          "and the downstream step30g/h refits.", flush=True)

    outp = G.guard(f"{G.THRU}/gate2_throughput_probe.tsv")
    ks = list(results[0].keys())
    with open(outp, "w") as f:
        f.write("\t".join(ks) + "\n")
        for r in results:
            f.write("\t".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]) for k in ks) + "\n")
    with open(G.guard(f"{G.THRU}/gate2_throughput_summary.json"), "w") as f:
        json.dump(dict(probe_rows=tot_rows, probe_secs=tot_secs, secs_per_weight_row=per_row,
                       ext2_weight_rows=EXT2_ROWS, est_hours_single_process=est_h,
                       bytes_per_variant_row=bpv,
                       est_total_read_GB=EXT2_ROWS * bpv / 1e9), f, indent=1)
    print(f"  -> {outp}", flush=True)


if __name__ == "__main__":
    main()
