#!/usr/bin/env python3
"""SD-drift calibration: AoU in-sample-SD vs 1kGP-EUR train-SD standardization."""
import os, csv, numpy as np, pgenlib

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
H = f"{BASE}/aou_handoff"
PGEN = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen")
OUT = f"{BASE}/aou_replication/results/step5_drift"
os.makedirs(OUT, exist_ok=True)

tm = {(r['sv_locus'], r['chrom'], r['pos'], r['ref'], r['alt']): float(r['train_maf_panel'])
      for r in csv.DictReader(open(f"{H}/train_maf.tsv"), delimiter='\t')}
weights = list(csv.DictReader(open(f"{H}/sv_weights_long.tsv"), delimiter='\t'))
by_chrom = {}
for r in weights:
    by_chrom.setdefault(r['chrom'], []).append(r)

sv_drift_num = {}; sv_drift_den = {}
sv_score_ins = {}; sv_score_trn = {}
sv_n = {}; sv_chrom = {}
N = None

for chrom in sorted(by_chrom, key=lambda c: (len(c), c)):
    wr = by_chrom[chrom]
    want = {r['pos'] for r in wr}
    offsets = [0]; pos_map = {}; vidx = 0
    with open(f"{PGEN}/acaf_threshold.{chrom}.pvar") as f:
        for line in f:
            if line[0] == '#':
                continue
            c = line.rstrip('\n').split('\t', 5)
            alts = c[4].split(',')
            offsets.append(offsets[-1] + 1 + len(alts))
            if c[1] in want:
                pos_map.setdefault(c[1], []).append((vidx, c[3], alts))
            vidx += 1
    offsets = np.array(offsets, dtype=np.uintp)
    if N is None:
        N = sum(1 for _ in open(f"{PGEN}/acaf_threshold.{chrom}.psam")) - 1
    reader = pgenlib.PgenReader(f"{PGEN}/acaf_threshold.{chrom}.pgen".encode(), raw_sample_ct=N, variant_ct=vidx, allele_idx_offsets=offsets)
    buf = np.empty(N, dtype=np.int32)
    for row in wr:
        sv, pos, ref, alt, beta = row['sv_locus'], row['pos'], row['ref'], row['alt'], float(row['beta_std'])
        ftr = tm.get((sv, chrom, pos, ref, alt))
        hit = None
        for (vi, pref, palts) in pos_map.get(pos, []):
            if pref == ref and alt in palts:
                hit = (vi, 1 + palts.index(alt)); break
            if pref == alt and ref in palts:
                hit = (vi, 0); break
        if hit is None or ftr is None:
            continue
        reader.read(hit[0], buf, allele_idx=hit[1])
        g = buf.astype(np.float64); g[g < 0] = np.nan
        faou = np.nanmean(g) / 2.0
        sd_aou = np.nanstd(g)
        if not (sd_aou > 0) or not (0 < faou < 1) or not (0 < ftr < 1):
            continue
        sd_trn = np.sqrt(2 * ftr * (1 - ftr))
        drift = np.sqrt(ftr * (1 - ftr) / (faou * (1 - faou)))
        sv_chrom[sv] = chrom
        sv_drift_num[sv] = sv_drift_num.get(sv, 0.0) + beta * beta * drift
        sv_drift_den[sv] = sv_drift_den.get(sv, 0.0) + beta * beta
        sv_n[sv] = sv_n.get(sv, 0) + 1
        gi = np.where(np.isnan(g), np.nanmean(g), g)
        if sv not in sv_score_ins:
            sv_score_ins[sv] = np.zeros(N); sv_score_trn[sv] = np.zeros(N)
        sv_score_ins[sv] += beta * gi / sd_aou
        sv_score_trn[sv] += beta * gi / sd_trn
    reader.close(); del offsets, pos_map

rows = []
for sv in sorted(sv_n):
    dfac = sv_drift_num[sv] / sv_drift_den[sv]
    a, b = sv_score_ins[sv], sv_score_trn[sv]
    corr = float(np.corrcoef(a, b)[0, 1])
    rows.append((sv, sv_chrom[sv], sv_n[sv], dfac, corr))
    del sv_score_ins[sv], sv_score_trn[sv]

rows.sort(key=lambda r: r[4])
with open(f"{OUT}/sv_sd_drift.tsv", 'w') as f:
    f.write("sv_locus\tchrom\tn_snps\tdrift_factor_betawmean\trescore_corr_insample_vs_train\n")
    for sv, ch, n, d, c in rows:
        f.write(f"{sv}\t{ch}\t{n}\t{d:.4f}\t{c:.4f}\n")

dvals = np.array([r[3] for r in rows]); cvals = np.array([r[4] for r in rows])
print(f"SVs: {len(rows)}")
print(f"drift_factor: median={np.median(dvals):.3f} range=[{dvals.min():.3f},{dvals.max():.3f}]")
print(f"rescore_corr: median={np.median(cvals):.4f} min={cvals.min():.4f} ({(cvals<0.99).sum()} SVs <0.99, {(cvals<0.95).sum()} <0.95)")
print("\nMost-affected SVs (lowest in-sample-vs-train re-score corr):")
for sv, ch, n, d, c in rows[:6]:
    print(f"  {sv}  n={n}  drift={d:.3f}  corr={c:.4f}")
print(f"\nwrote {OUT}/sv_sd_drift.tsv")
