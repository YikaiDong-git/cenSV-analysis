#!/usr/bin/env python3
"""Step 41 — SCOREABILITY PRE-FLIGHT for the ext2 (11,129-predictor) handoff."""
import os, csv, glob
from collections import defaultdict

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
EXT2 = f"{BASE}/phewide_all_handoff/full_11129_ext2"
BIM = f"{HOME}/workspace/working/jwein22/centromere_phewas/results/vds_pericentromeric_snps_min_ac1000"
OUT = f"{BASE}/aou_replication/results_vds/step41_ext2_scoreability"
OLD_PREFLIGHT = f"{BASE}/aou_replication/results_vds/step27_phewide/scoreability_preflight.tsv"
os.makedirs(OUT, exist_ok=True)


def load_meta():
    m = {}
    for r in csv.DictReader(open(f"{EXT2}/sv_metadata_ext2.tsv"), delimiter="\t"):
        m[r["sv_locus"]] = r
    return m


def load_weights():
    by_chrom = defaultdict(list); n_total = defaultdict(int)
    files = sorted(glob.glob(f"{EXT2}/sv_weights_long_ext2_chr*.tsv"))
    assert files, f"no ext2 weight files under {EXT2}"
    for f in files:
        for r in csv.DictReader(open(f), delimiter="\t"):
            ch = r["chrom"] if str(r["chrom"]).startswith("chr") else f"chr{r['chrom']}"
            by_chrom[ch].append((r["sv_locus"], r["pos"], r["ref"], r["alt"]))
            n_total[r["sv_locus"]] += 1
    return by_chrom, n_total, len(files)


def index_bim(path, want_pos):
    pos_map = defaultdict(list)
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) >= 6 and c[3] in want_pos:
                pos_map[c[3]].append((c[4], c[5]))
    return pos_map


def main():
    meta = load_meta()
    by_chrom, n_total, nfiles = load_weights()
    print(f"ext2: {sum(n_total.values()):,} weight-SNP rows across {len(n_total):,} predictors "
          f"({nfiles} files); metadata rows {len(meta):,}", flush=True)
    missing_meta = set(n_total) - set(meta)
    assert not missing_meta, f"{len(missing_meta)} predictors in weights but not metadata, e.g. {list(missing_meta)[:3]}"

    n_present = defaultdict(int)
    for ch in sorted(by_chrom):
        path = f"{BIM}/{ch}/{ch}.pericentromeric_snps.min_ac1000.bim"
        if not os.path.exists(path):
            print(f"  {ch}: NO .bim -> all {len(by_chrom[ch]):,} weight-SNPs UNSCOREABLE", flush=True)
            continue
        want = {w[1] for w in by_chrom[ch]}
        pos_map = index_bim(path, want)
        hit = 0
        for sv, pos, ref, alt in by_chrom[ch]:
            for (a1, a2) in pos_map.get(pos, []):
                if (a1 == alt and a2 == ref) or (a1 == ref and a2 == alt):
                    n_present[sv] += 1; hit += 1; break
        print(f"  {ch}: {hit:7,}/{len(by_chrom[ch]):7,} weight-SNPs matched in the AoU VDS set", flush=True)
        del pos_map

    outp = f"{OUT}/ext2_scoreability.tsv"
    with open(outp, "w") as fo:
        fo.write("sv_locus\ttier\talready_shipped\tchrom\tcv_r2\tn_weight_snps\tn_snps_in_aou\tfrac_covered\n")
        for sv in sorted(n_total):
            m = meta[sv]; nt = n_total[sv]; npv = n_present.get(sv, 0)
            fo.write(f"{sv}\t{m.get('tier','')}\t{m.get('already_shipped','')}\t{m.get('chrom','')}\t"
                     f"{m.get('cv_r2','')}\t{nt}\t{npv}\t{npv/nt:.4f}\n")

    rows = list(csv.DictReader(open(outp), delimiter="\t"))
    print("\n=== SCOREABILITY by ext2 tier ===")
    print(f"{'tier':28}{'N':>7}{'>=1 SNP':>9}{'>=5 SNP':>9}{'>=50%':>8}{'ZERO':>7}{'med frac':>10}")
    for tier in sorted({r["tier"] for r in rows}):
        s = [r for r in rows if r["tier"] == tier]
        cov = [int(r["n_snps_in_aou"]) for r in s]
        frac = sorted(float(r["frac_covered"]) for r in s)
        med = frac[len(frac) // 2] if frac else 0.0
        print(f"{tier:28}{len(s):>7,}{sum(c>=1 for c in cov):>9,}{sum(c>=5 for c in cov):>9,}"
              f"{sum(f>=0.5 for f in frac):>8,}{sum(c==0 for c in cov):>7,}{med:>10.3f}")
    new = [r for r in rows if r["already_shipped"] == "no"]
    cov = [int(r["n_snps_in_aou"]) for r in new]
    print(f"\nNEW predictors (not yet scored): {len(new):,} | scoreable(>=1 SNP) {sum(c>=1 for c in cov):,} "
          f"| >=5 SNPs {sum(c>=5 for c in cov):,} | ZERO {sum(c==0 for c in cov):,}")
    hi = [r for r in new if int(r["n_snps_in_aou"]) >= 5 and r["cv_r2"] not in ("", "NA")
          and float(r["cv_r2"]) >= 0.30]
    print(f"HIGH-YIELD slice (new & cv_r2>=0.30 & >=5 AoU SNPs): {len(hi):,} predictors")

    if os.path.exists(OLD_PREFLIGHT):
        old = {r["sv_locus"]: int(r["n_snps_in_vds"]) for r in csv.DictReader(open(OLD_PREFLIGHT), delimiter="\t")}
        shipped = [r for r in rows if r["already_shipped"] == "yes"]
        both = [r for r in shipped if r["sv_locus"] in old]
        fewer = [r for r in both if int(r["n_snps_in_aou"]) < old[r["sv_locus"]]]
        same = sum(1 for r in both if int(r["n_snps_in_aou"]) == old[r["sv_locus"]])
        print(f"\n=== CONTROL vs step27 preflight (old weights) on {len(both):,} already-shipped predictors ===")
        print(f"  more AoU SNPs now: {len(both)-same-len(fewer):,} | identical: {same:,} | FEWER: {len(fewer):,}"
              f"  <- FEWER must be 0 (ext2 is a superset)")
        if fewer:
            print("  VIOLATIONS:", [(r['sv_locus'], old[r['sv_locus']], r['n_snps_in_aou']) for r in fewer[:5]])
    print(f"\nwrote {outp}")


if __name__ == "__main__":
    main()
