#!/usr/bin/env python3
"""step97d (LANE S, gap-closure) - tag-SNP SPAN + COVERAGE for ALL 123 target SVs."""
import os, csv, statistics
from collections import defaultdict, OrderedDict

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
EXT2 = f"{BASE}/phewide_all_handoff/full_11129_ext2"
BED  = f"{HOME}/workspace/working/jwein22/centromere_phewas/results/vds_pericentromeric_snps_min_ac1000"
RV   = f"{BASE}/aou_replication/results_vds"
OUT  = f"{RV}/step97_discovery_requests"
HITS = f"{RV}/step30_phewide_ext2/disease_matrices_firth/disease_hits_final.tsv"
DQ   = f"{RV}/step28_phewide_ext2/dosage_qc_ext2.tsv"

CLUSTER_GAP = 100_000
BLOCK_GAP   = 1_000_000

tiers = {}
with open(HITS, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        q = float(r["q_bh_tier"])
        if q < 0.10:
            sv = r["sv_locus"]
            t = "q<0.05" if q < 0.05 else "0.05<=q<0.10"
            if tiers.get(sv) != "q<0.05":
                tiers[sv] = t
TARGETS = sorted(tiers)
assert len(TARGETS) == 123, f"expected 123 target SVs, derived {len(TARGETS)}"

by_chrom = defaultdict(set)
for sv in TARGETS:
    by_chrom[sv.split("_")[0]].add(sv)

rows_out, bi_out, loc_out = [], [], []
for ch in sorted(by_chrom, key=lambda c: (len(c), c)):
    want_sv = by_chrom[ch]

    w = defaultdict(list)
    wpath = f"{EXT2}/sv_weights_long_ext2_{ch}.tsv"
    with open(wpath, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            sv = r["sv_locus"]
            if sv in want_sv:
                assert r["chrom"] == ch, f"chrom mismatch {sv}"
                w[sv].append((r["pos"], r["ref"], r["alt"], float(r["beta_std"])))
    assert set(w) == want_sv, f"{ch}: missing weights for {sorted(want_sv - set(w))}"

    want_pos = {p for rs in w.values() for (p, _, _, _) in rs}
    pos_map = defaultdict(list)
    nvar = 0
    blocks = []
    blk_lo = blk_hi = None
    unsorted_bim = 0
    with open(f"{BED}/{ch}/{ch}.pericentromeric_snps.min_ac1000.bim") as fh:
        for line in fh:
            c = line.split()
            if c[3] in want_pos:
                pos_map[c[3]].append((nvar, c[4], c[5], c[1]))
            p_ = int(c[3])
            if blk_lo is None:
                blk_lo = blk_hi = p_
            elif p_ < blk_hi:
                unsorted_bim += 1
            elif p_ - blk_hi > BLOCK_GAP:
                blocks.append((blk_lo, blk_hi)); blk_lo = blk_hi = p_
            else:
                blk_hi = p_
            nvar += 1
    if blk_lo is not None:
        blocks.append((blk_lo, blk_hi))
    assert unsorted_bim == 0, f"{ch}: bim not position-sorted ({unsorted_bim} inversions)"
    def in_covered(p):
        return any(a <= p <= b for (a, b) in blocks)
    cov_bp = sum(b - a for (a, b) in blocks)

    bad_id = 0
    for p, ents in pos_map.items():
        for (vi, a1, a2, sid) in ents:
            f_ = sid.split(":")
            if len(f_) != 4 or f_[1] != p or f_[2] != a2 or f_[3] != a1:
                bad_id += 1

    for sv in sorted(want_sv):
        recs = w[sv]
        tot_abs = sum(abs(b) for (_, _, _, b) in recs)
        positions_all, positions_pres = [], []
        mass_pres = 0.0
        n_pres = 0
        cls = defaultdict(int); clsm = defaultdict(float)
        absloc = defaultdict(int); abslocm = defaultdict(float)
        for (pos, ref, alt, b) in recs:
            ents = pos_map.get(pos, [])
            hit = swapped = None
            for (vi, a1, a2, sid) in ents:
                if a1 == alt and a2 == ref:
                    hit, swapped = sid, False; break
                if a1 == ref and a2 == alt:
                    hit, swapped = sid, True;  break
            ip = int(pos)
            positions_all.append(ip)
            n_at = len(ents)
            refs_at = {e[2] for e in ents}
            if n_at == 0:
                site = "NA_no_record"
            elif n_at == 1:
                site = "single_record"
            elif len(refs_at) == 1:
                site = "multi_record_shared_ref"
            else:
                site = "multi_record_ref_differs"
            if hit is not None:
                n_pres += 1
                mass_pres += abs(b)
                positions_pres.append(ip)
                k = ("multiallelic_recovered"
                     if site in ("multi_record_shared_ref", "multi_record_ref_differs")
                     else ("allele_swapped" if swapped else "exact"))
            else:
                k = "absent_pos" if n_at == 0 else "absent_no_allele_match"
                if in_covered(ip):
                    absloc["absent_inside_covered_block"] += 1
                    abslocm["absent_inside_covered_block"] += abs(b)
                else:
                    absloc["absent_in_resource_gap"] += 1
                    abslocm["absent_in_resource_gap"] += abs(b)
            cls[k] += 1; clsm[k] += abs(b)

        positions_all.sort(); positions_pres.sort()
        gaps = [positions_all[i+1] - positions_all[i] for i in range(len(positions_all)-1)]
        nclust = 1 + sum(1 for g in gaps if g > CLUSTER_GAP) if positions_all else 0

        rows_out.append(OrderedDict(
            sv_locus=sv, chrom=ch,
            n_weight_snps=len(recs),
            n_present_in_bim=n_pres,
            frac_snps_present=f"{n_pres/len(recs):.6f}",
            sum_abs_beta_total=f"{tot_abs:.10g}",
            sum_abs_beta_present=f"{mass_pres:.10g}",
            frac_weight_mass_present=(f"{mass_pres/tot_abs:.8f}" if tot_abs > 0 else "NOT_AVAILABLE_zero_total_mass"),
            span_bp=(positions_all[-1] - positions_all[0]) if positions_all else "",
            span_bp_present=((positions_pres[-1] - positions_pres[0]) if positions_pres else
                             ("NOT_AVAILABLE_no_snp_present" if n_pres == 0 else 0)),
            largest_gap_bp=(max(gaps) if gaps else 0),
            n_clusters_100kb=nclust,
            min_pos=(positions_all[0] if positions_all else ""),
            max_pos=(positions_all[-1] if positions_all else ""),
            median_inter_snp_gap_bp=(f"{statistics.median(gaps):.1f}" if gaps else "NOT_AVAILABLE_single_snp"),
        ))
        n_abs_tot = cls["absent_pos"] + cls["absent_no_allele_match"]
        loc_out.append(OrderedDict(
            sv_locus=sv, chrom=ch, n_weight_snps=len(recs), n_absent=n_abs_tot,
            n_absent_inside_covered_block=absloc["absent_inside_covered_block"],
            n_absent_in_resource_gap=absloc["absent_in_resource_gap"],
            mass_absent_inside_covered_block=f"{abslocm['absent_inside_covered_block']/tot_abs:.6f}" if tot_abs else "NA",
            mass_absent_in_resource_gap=f"{abslocm['absent_in_resource_gap']/tot_abs:.6f}" if tot_abs else "NA",
            resource_n_blocks=len(blocks),
            resource_covered_bp=cov_bp,
            resource_span=f"{blocks[0][0]}-{blocks[-1][1]}",
            resource_blocks=";".join(f"{a}-{b}" for a, b in blocks),
        ))
        bi_out.append(OrderedDict(
            sv_locus=sv, chrom=ch, n_weight_snps=len(recs), bim_id_format_violations=bad_id,
            n_exact=cls["exact"], n_allele_swapped=cls["allele_swapped"],
            n_multiallelic_recovered=cls["multiallelic_recovered"],
            n_absent_pos=cls["absent_pos"], n_absent_no_allele_match=cls["absent_no_allele_match"],
            mass_exact=f"{clsm['exact']/tot_abs:.6f}" if tot_abs else "NA",
            mass_allele_swapped=f"{clsm['allele_swapped']/tot_abs:.6f}" if tot_abs else "NA",
            mass_multiallelic_recovered=f"{clsm['multiallelic_recovered']/tot_abs:.6f}" if tot_abs else "NA",
            mass_absent=f"{(clsm['absent_pos']+clsm['absent_no_allele_match'])/tot_abs:.6f}" if tot_abs else "NA",
            tier=tiers[sv],
        ))
    del pos_map, w

rows_out.sort(key=lambda r: (r["chrom"], r["sv_locus"]))
bi_out.sort(key=lambda r: (r["chrom"], r["sv_locus"]))
loc_out.sort(key=lambda r: (r["chrom"], r["sv_locus"]))

def write(path, rows):
    with open(path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                            lineterminator="\n", extrasaction="raise")
        wr.writeheader(); wr.writerows(rows)
    return os.path.getsize(path)

s1 = write(f"{OUT}/tagsnp_coverage_123svs.tsv", rows_out)
s2 = write(f"{OUT}/tagsnp_coverage_biallelic.tsv", bi_out)
s3 = write(f"{OUT}/tagsnp_absent_location.tsv", loc_out)

dq = {}
with open(DQ, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        dq[r["sv_locus"]] = r

diffs = []
for r in rows_out:
    d = dq.get(r["sv_locus"])
    if d is None:
        diffs.append((r["sv_locus"], None, None, None)); continue
    mine = float(r["frac_weight_mass_present"])
    theirs = 1.0 - float(d["unscoreable_mass_frac"])
    diffs.append((r["sv_locus"], mine, theirs, abs(mine - theirs)))
have = [d for d in diffs if d[3] is not None]
have.sort(key=lambda x: -x[3])

cnt_diffs = []
for r in rows_out:
    d = dq.get(r["sv_locus"])
    if d:
        cnt_diffs.append((r["sv_locus"], r["n_weight_snps"] - int(d["n_total"]),
                          r["n_present_in_bim"] - int(d["n_used"])))

print(f"targets derived: {len(TARGETS)}  (q<0.05 SVs: {sum(1 for v in tiers.values() if v=='q<0.05')})")
print(f"wrote {OUT}/tagsnp_coverage_123svs.tsv  {s1} bytes")
print(f"wrote {OUT}/tagsnp_coverage_biallelic.tsv  {s2} bytes")
print(f"wrote {OUT}/tagsnp_absent_location.tsv  {s3} bytes")
print()
def q(vals, p):
    v = sorted(vals); import math
    return v[min(len(v)-1, max(0, int(round(p*(len(v)-1)))))]
fs = [float(r["frac_snps_present"]) for r in rows_out]
fm = [float(r["frac_weight_mass_present"]) for r in rows_out]
for nm, v in (("frac_snps_present", fs), ("frac_weight_mass_present", fm)):
    print(f"{nm}: min {min(v):.4f} p10 {q(v,.10):.4f} p25 {q(v,.25):.4f} median {statistics.median(v):.4f} "
          f"p75 {q(v,.75):.4f} p90 {q(v,.90):.4f} max {max(v):.4f} mean {statistics.mean(v):.4f}")
    print(f"   ==1.0: {sum(1 for x in v if x==1.0)}  >=0.95: {sum(1 for x in v if x>=0.95)}  "
          f">=0.90: {sum(1 for x in v if x>=0.90)}  <0.75: {sum(1 for x in v if x<0.75)}  <0.50: {sum(1 for x in v if x<0.50)}")
ntot = sum(r["n_weight_snps"] for r in rows_out)
npre = sum(r["n_present_in_bim"] for r in rows_out)
mtot = sum(float(r["sum_abs_beta_total"]) for r in rows_out)
mpre = sum(float(r["sum_abs_beta_present"]) for r in rows_out)
print(f"\npooled over 123 SVs: {npre}/{ntot} weight-SNP slots present = {npre/ntot:.4f}")
print(f"pooled mass-weighted (sum over SVs of sum|beta|): {mpre/mtot:.4f}  "
      f"[note: cross-SV pooling of unnormalised sum|beta| is dominated by large-|beta| SVs]")
print(f"unweighted mean of per-SV frac_weight_mass_present: {statistics.mean(fm):.4f}")

print("\nworst 12 by frac_weight_mass_present:")
for r in sorted(rows_out, key=lambda x: float(x["frac_weight_mass_present"]))[:12]:
    print(f"  {r['sv_locus']:<52} snps {r['n_present_in_bim']:>4}/{r['n_weight_snps']:<5} "
          f"fsnp {float(r['frac_snps_present']):.4f}  fmass {float(r['frac_weight_mass_present']):.4f}  "
          f"span {r['span_bp']:>10}  clus {r['n_clusters_100kb']}  tier {tiers[r['sv_locus']]}")

print("\nCROSS-CHECK frac_weight_mass_present  vs  1 - unscoreable_mass_frac (dosage_qc_ext2.tsv):")
print(f"  SVs compared: {len(have)}/123   missing from dosage_qc: {sum(1 for d in diffs if d[3] is None)}")
print(f"  MAX ABS DISCREPANCY: {have[0][3]:.3e}   on {have[0][0]}  (mine {have[0][1]:.8f} vs theirs {have[0][2]:.8f})")
print("  top 5:")
for sv, a, b, d in have[:5]:
    print(f"    {sv:<52} mine {a:.8f}  qc {b:.8f}  |diff| {d:.3e}")
print(f"  n with |diff| > 1e-6 : {sum(1 for d in have if d[3] > 1e-6)}")
print(f"  n with |diff| > 1e-9 : {sum(1 for d in have if d[3] > 1e-9)}")
bad_n = [c for c in cnt_diffs if c[1] != 0 or c[2] != 0]
print(f"  n_total / n_used count mismatches vs dosage_qc: {len(bad_n)}"
      + ("" if not bad_n else f"  e.g. {bad_n[:5]}"))

sp = [r["span_bp"] for r in rows_out]
print(f"\nspan_bp: min {min(sp):,} median {int(statistics.median(sp)):,} max {max(sp):,}")
gg = [r["largest_gap_bp"] for r in rows_out]
print(f"largest_gap_bp: min {min(gg):,} median {int(statistics.median(gg)):,} max {max(gg):,}")
cc = [r["n_clusters_100kb"] for r in rows_out]
from collections import Counter
print(f"n_clusters_100kb: {dict(sorted(Counter(cc).items()))}")
print(f"n_weight_snps: min {min(r['n_weight_snps'] for r in rows_out)} "
      f"median {int(statistics.median([r['n_weight_snps'] for r in rows_out]))} "
      f"max {max(r['n_weight_snps'] for r in rows_out)}")

tb = sum(b["n_exact"]+b["n_allele_swapped"] for b in bi_out)
tm = sum(b["n_multiallelic_recovered"] for b in bi_out)
ta = sum(b["n_absent_pos"]+b["n_absent_no_allele_match"] for b in bi_out)
print(f"\nsite-class split pooled over {ntot} weight-SNP slots: "
      f"single-record(biallelic-in-resource) {tb} ({tb/ntot:.4f}), "
      f"multi-record-recovered {tm} ({tm/ntot:.4f}), absent {ta} ({ta/ntot:.4f})")
print(f"bim_id_format_violations (should be 0): {sum(b['bim_id_format_violations'] for b in bi_out)}")
ag = sum(l["n_absent_in_resource_gap"] for l in loc_out); ai = sum(l["n_absent_inside_covered_block"] for l in loc_out)
print(f"\nWHERE the {ag+ai} absent weight-SNPs fall: {ag} ({ag/(ag+ai):.4f}) in a >1 Mb GAP with no bim record at all, "
      f"{ai} ({ai/(ag+ai):.4f}) inside a covered block (present region, variant simply not in the AC>=1000 resource)")
from collections import Counter as _C
gc = _C(l["chrom"] for l in loc_out for _ in range(l["n_absent_in_resource_gap"]))
print("  gap-absent by chrom:", dict(sorted(gc.items(), key=lambda kv: -kv[1])))
