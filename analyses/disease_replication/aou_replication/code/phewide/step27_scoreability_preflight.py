#!/usr/bin/env python3
import os, csv, glob
from collections import defaultdict

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
PH = f"{BASE}/phewide_all_handoff"
BED = f"{HOME}/workspace/working/jwein22/centromere_phewas/results/vds_pericentromeric_snps_min_ac1000"
OUT = f"{BASE}/aou_replication/results_vds/step27_phewide"
os.makedirs(OUT, exist_ok=True)


def load_weights():
    rows = []
    for f in sorted(glob.glob(f"{PH}/sv_weights_long_chr*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            rows.append((r["sv_locus"], r["chrom"], r["pos"], r["ref"], r["alt"], "signal"))
    for f in sorted(glob.glob(f"{PH}/discovery_null_ext/sv_weights_long_ext_chr*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            rows.append((r["sv_locus"], r["chrom"], r["pos"], r["ref"], r["alt"], "null_ext"))
    return rows


def index_bim(path, want_pos):
    pos_map = defaultdict(list)
    with open(path) as f:
        for line in f:
            c = line.split()
            if c[3] in want_pos:
                pos_map[c[3]].append((c[4], c[5]))
    return pos_map


def main():
    weights = load_weights()
    by_chrom = defaultdict(list)
    for w in weights:
        ch = w[1] if str(w[1]).startswith("chr") else f"chr{w[1]}"
        by_chrom[ch].append(w)
    tier = {}
    n_total = defaultdict(int)
    for sv, ch, pos, ref, alt, tr in weights:
        tier[sv] = tr
        n_total[sv] += 1
    print(f"{len(weights):,} weight-SNP rows across {len(tier):,} SVs "
          f"({sum(v=='signal' for v in tier.values())} signal / {sum(v=='null_ext' for v in tier.values())} null-ext)", flush=True)

    n_present = defaultdict(int)
    for ch in sorted(by_chrom):
        bim = f"{BED}/{ch}/{ch}.pericentromeric_snps.min_ac1000.bim"
        if not os.path.exists(bim):
            print(f"  {ch}: NO .bim -> all {len(by_chrom[ch])} weight-SNPs unscoreable", flush=True)
            continue
        want = {w[2] for w in by_chrom[ch]}
        pos_map = index_bim(bim, want)
        hit = 0
        for sv, _c, pos, ref, alt, _t in by_chrom[ch]:
            for (a1, a2) in pos_map.get(pos, []):
                if (a1 == alt and a2 == ref) or (a1 == ref and a2 == alt):
                    n_present[sv] += 1; hit += 1; break
        print(f"  {ch}: {hit:6d}/{len(by_chrom[ch]):6d} weight-SNPs matched in VDS", flush=True)
        del pos_map

    with open(f"{OUT}/scoreability_preflight.tsv", "w") as fo:
        fo.write("sv_locus\ttier\tn_weight_snps\tn_snps_in_vds\tfrac_covered\n")
        for sv in sorted(tier):
            nt = n_total[sv]; npv = n_present.get(sv, 0)
            fo.write(f"{sv}\t{tier[sv]}\t{nt}\t{npv}\t{npv/nt:.4f}\n")

    print("\n=== SCOREABILITY by tier (SV is scoreable if >=1 weight-SNP present) ===")
    for tr in ("signal", "null_ext"):
        svs = [sv for sv in tier if tier[sv] == tr]
        cov = [n_present.get(sv, 0) for sv in svs]
        frac = [n_present.get(sv, 0) / n_total[sv] for sv in svs]
        scoreable = sum(1 for c in cov if c >= 1)
        ge_half = sum(1 for fr in frac if fr >= 0.5)
        zero = sum(1 for c in cov if c == 0)
        med = sorted(frac)[len(frac) // 2] if frac else 0
        print(f"  {tr:9}: {len(svs):5d} SVs | scoreable(>=1 SNP) {scoreable:5d} | "
              f">=50% SNPs {ge_half:5d} | 0 SNPs {zero:4d} | median frac_covered {med:.3f}", flush=True)
    allsv = list(tier)
    print(f"  TOTAL    : {len(allsv):5d} SVs | scoreable {sum(1 for sv in allsv if n_present.get(sv,0)>=1):5d}")
    print(f"\nwrote {OUT}/scoreability_preflight.tsv")


if __name__ == "__main__":
    main()
