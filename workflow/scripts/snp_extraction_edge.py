#!/usr/bin/env python3
"""Centromere-edge-outward SNP extraction."""
import argparse
import logging
import sys
from pathlib import Path

import cyvcf2
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def load_sample_list(sample_path):
    with open(sample_path) as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith('#')]


def load_centromeres(bed_path):
    df = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end'])
    return df.set_index('chrom')


def load_heterochromatin(bed_path):
    df = pd.read_csv(bed_path, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'type'])
    return df


def compute_extended_edges(cen_start, cen_end, het_df, chrom, proximity=1_000_000):
    left_edge = cen_start
    right_edge = cen_end

    chrom_het = het_df[het_df['chrom'] == chrom]

    for _, het in chrom_het.iterrows():
        h_start, h_end = het['start'], het['end']

        if het['type'] == 'short_arm':
            logger.info(f"  Skipping short_arm region {h_start:,}-{h_end:,} (would eliminate wing)")
            continue

        if h_start >= cen_start and h_end <= cen_end:
            continue

        if h_end <= cen_start:
            gap = cen_start - h_end
            if gap <= proximity:
                left_edge = min(left_edge, h_start)
                logger.info(f"  Extending left edge: het {h_start:,}-{h_end:,} "
                            f"({het['type']}, gap={gap/1e3:.0f}kb) → left_edge={left_edge:,}")

        elif h_start >= cen_end:
            gap = h_start - cen_end
            if gap <= proximity:
                right_edge = max(right_edge, h_end)
                logger.info(f"  Extending right edge: het {h_start:,}-{h_end:,} "
                            f"({het['type']}, gap={gap/1e3:.0f}kb) → right_edge={right_edge:,}")

        elif h_start < cen_start < h_end:
            left_edge = min(left_edge, h_start)
            logger.info(f"  Extending left edge: het overlaps boundary "
                        f"{h_start:,}-{h_end:,} → left_edge={left_edge:,}")

        elif h_start < cen_end < h_end:
            right_edge = max(right_edge, h_end)
            logger.info(f"  Extending right edge: het overlaps boundary "
                        f"{h_start:,}-{h_end:,} → right_edge={right_edge:,}")

    return left_edge, right_edge


def parse_vcf_snps(vcf_path, region, sample_ids, min_maf, max_missing):
    reader = cyvcf2.VCF(str(vcf_path), gts012=True, samples=sample_ids)
    samples = list(reader.samples)
    n_samples = len(samples)
    logger.info(f"  Loaded {n_samples} samples from VCF")

    gt_list = []
    variants_info = []

    for v in reader(region):
        if v.FILTER is not None:
            continue

        if len(v.ALT) != 1 or len(v.REF) != 1 or len(v.ALT[0]) != 1:
            continue

        gt = v.gt_types.astype(np.int8)
        known_mask = gt != 3
        n_known = np.sum(known_mask)

        if n_known == 0:
            continue

        ac = np.sum(gt[known_mask])
        an = 2 * n_known
        af = ac / an
        maf = min(af, 1.0 - af)

        if maf < min_maf:
            continue

        frac_missing = 1.0 - (n_known / n_samples)
        if frac_missing > max_missing:
            continue

        gt_list.append(gt)
        variants_info.append({
            'chrom': v.CHROM,
            'pos': v.POS,
            'ref': v.REF,
            'alt': v.ALT[0],
            'id': v.ID or f"{v.CHROM}_{v.POS}_{v.REF}_{v.ALT[0]}",
            'ac': int(ac),
            'an': int(an),
            'af': round(af, 6),
            'maf': round(maf, 6),
            'n_miss': int(n_samples - n_known)
        })

    reader.close()

    if not gt_list:
        return None, None, samples

    matrix = np.vstack(gt_list)
    return matrix, variants_info, samples


def select_edge_outward(variants_info, genotype_matrix, left_edge, right_edge, n_per_side):
    positions = np.array([v['pos'] for v in variants_info])

    left_mask = positions < left_edge
    left_indices = np.where(left_mask)[0]
    left_distances = left_edge - positions[left_indices]
    left_order = np.argsort(left_distances)
    left_selected = left_indices[left_order[:n_per_side]]

    right_mask = positions > right_edge
    right_indices = np.where(right_mask)[0]
    right_distances = positions[right_indices] - right_edge
    right_order = np.argsort(right_distances)
    right_selected = right_indices[right_order[:n_per_side]]

    n_left = len(left_selected)
    n_right = len(right_selected)

    inside_mask = (positions >= left_edge) & (positions <= right_edge)
    n_inside = np.sum(inside_mask)

    if n_left < n_per_side:
        logger.warning(f"  Left flank: only {n_left}/{n_per_side} SNPs available outside gap")
    if n_right < n_per_side:
        logger.warning(f"  Right flank: only {n_right}/{n_per_side} SNPs available outside gap")
    if n_inside > 0:
        logger.info(f"  Excluded {n_inside} SNPs inside gap [{left_edge:,}-{right_edge:,}]")

    selected = np.sort(np.concatenate([left_selected, right_selected]))

    filtered_variants = [variants_info[i] for i in selected]
    filtered_matrix = genotype_matrix[selected]

    logger.info(f"  Selected {n_left} left + {n_right} right = {len(selected)} total SNPs")
    if n_left > 0:
        left_span = (left_edge - positions[left_selected].min()) / 1e6
        logger.info(f"  Left span: {left_span:.2f} Mb outward from left edge")
    if n_right > 0:
        right_span = (positions[right_selected].max() - right_edge) / 1e6
        logger.info(f"  Right span: {right_span:.2f} Mb outward from right edge")

    return filtered_variants, filtered_matrix, n_left, n_right


def save_outputs(matrix, variants_info, samples, chrom, output_dir):
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    if matrix is None or len(variants_info) == 0:
        logger.warning(f"  No variants to save for {chrom}")
        return None

    df_var = pd.DataFrame(variants_info)
    np.savez_compressed(out_path / f"{chrom}_snp_genotypes.npz", genotypes=matrix)
    df_var.to_csv(out_path / f"{chrom}_snp_variants.tsv", sep='\t', index=False)
    with open(out_path / f"{chrom}_snp_samples.txt", 'w') as fh:
        fh.write('\n'.join(samples) + '\n')
    return df_var


def main():
    parser = argparse.ArgumentParser(
        description="Extract peri-centromeric SNPs using edge-outward selection"
    )
    parser.add_argument('--vcf-dir', required=True, type=Path, help="VCF directory")
    parser.add_argument('--vcf-pattern', required=True, help="VCF filename pattern with {chrom}")
    parser.add_argument('--samples', required=True, type=Path, help="Sample IDs file")
    parser.add_argument('--centromere-bed', required=True, type=Path, help="Centromere BED file")
    parser.add_argument('--het-bed', required=True, type=Path,
                        help="Heterochromatin/gap BED file for edge extension")
    parser.add_argument('--output', required=True, type=Path, help="Output directory")
    parser.add_argument('--n-snps-per-side', type=int, default=1000,
                        help="Number of SNPs to select per flank (default: 1000)")
    parser.add_argument('--max-search-distance', type=int, default=20_000_000,
                        help="Max distance from edge to search for SNPs (default: 20Mb)")
    parser.add_argument('--proximity-threshold', type=int, default=1_000_000,
                        help="Max gap between het region and centromere to consider adjacent (default: 1Mb)")
    parser.add_argument('--chromosomes', default='all',
                        help="Comma-separated chromosome list or 'all'")
    parser.add_argument('--min-maf', type=float, default=0.01, help="Min MAF (default: 0.01)")
    parser.add_argument('--max-missing', type=float, default=0.05,
                        help="Max missing fraction (default: 0.05)")
    parser.add_argument('--het-aware', action='store_true',
                        help="Extend centromere edges past adjacent heterochromatin blocks "
                             "before SNP search. Default OFF (use raw centromere boundaries).")

    args = parser.parse_args()

    if args.chromosomes.lower() == 'all':
        chroms = [str(c) for c in range(1, 23)] + ['X']
    else:
        chroms = [c.strip().replace('chr', '') for c in args.chromosomes.split(',')]

    sample_ids = load_sample_list(args.samples)
    logger.info(f"Loaded {len(sample_ids)} target sample IDs")

    centromeres = load_centromeres(args.centromere_bed)
    logger.info(f"Loaded centromere coordinates for {len(centromeres)} chromosomes")

    het_df = load_heterochromatin(args.het_bed)
    logger.info(f"Loaded {len(het_df)} heterochromatin/gap entries")

    summary = []
    for chrom_num in chroms:
        chrom = f"chr{chrom_num}"
        logger.info(f"\nProcessing {chrom}...")

        if chrom not in centromeres.index:
            logger.warning(f"  {chrom} not found in centromere BED, skipping")
            continue

        cen_start = centromeres.loc[chrom, 'start']
        cen_end = centromeres.loc[chrom, 'end']

        if args.het_aware:
            left_edge, right_edge = compute_extended_edges(
                cen_start, cen_end, het_df, chrom, args.proximity_threshold
            )
        else:
            left_edge, right_edge = cen_start, cen_end

        if left_edge != cen_start or right_edge != cen_end:
            logger.info(f"  Centromere: {cen_start:,}-{cen_end:,}")
            logger.info(f"  Extended:   {left_edge:,}-{right_edge:,}")
        else:
            mode = "het-aware no extension" if args.het_aware else "het-aware OFF"
            logger.info(f"  Centromere: {cen_start:,}-{cen_end:,} ({mode})")

        search_start = max(0, left_edge - args.max_search_distance)
        search_end = right_edge + args.max_search_distance
        region = f"{chrom}:{search_start}-{search_end}"
        logger.info(f"  Search region: {region}")

        vcf_path = args.vcf_dir / args.vcf_pattern.format(chrom=chrom_num)
        if not vcf_path.exists():
            logger.warning(f"  VCF not found: {vcf_path}")
            continue

        matrix, variants, samples = parse_vcf_snps(
            vcf_path, region, sample_ids, args.min_maf, args.max_missing
        )

        if matrix is None:
            logger.warning(f"  No SNPs found in search region for {chrom}")
            save_outputs(None, [], samples, chrom, args.output)
            continue

        logger.info(f"  Found {matrix.shape[0]} SNPs in search region after QC filters")

        sel_variants, sel_matrix, n_left, n_right = select_edge_outward(
            variants, matrix, left_edge, right_edge, args.n_snps_per_side
        )

        df_var = save_outputs(sel_matrix, sel_variants, samples, chrom, args.output)

        summary.append({
            'chromosome': chrom,
            'cen_start': cen_start,
            'cen_end': cen_end,
            'left_edge': left_edge,
            'right_edge': right_edge,
            'left_extended': left_edge != cen_start,
            'right_extended': right_edge != cen_end,
            'n_snps_search_region': matrix.shape[0],
            'n_snps_left': n_left,
            'n_snps_right': n_right,
            'n_snps_total': n_left + n_right,
            'n_samples': sel_matrix.shape[1] if sel_matrix is not None else 0,
            'mean_maf': df_var['maf'].mean() if df_var is not None else 0,
        })
        logger.info(f"  Saved {chrom}: {n_left + n_right} SNPs")

    if summary:
        df_summary = pd.DataFrame(summary)
        df_summary.to_csv(args.output / "extraction_summary.tsv", sep='\t', index=False)
        logger.info(f"\nSummary saved to {args.output / 'extraction_summary.tsv'}")

        print("\n" + "=" * 110)
        print("EDGE-OUTWARD EXTRACTION SUMMARY")
        print("=" * 110)
        print(f"{'Chr':>5}  {'Cen Start':>12}  {'Cen End':>12}  {'L Edge':>12}  {'R Edge':>12}  "
              f"{'Search':>7}  {'Left':>6}  {'Right':>6}  {'Total':>6}  {'MAF':>6}")
        for _, row in df_summary.iterrows():
            ext = ""
            if row['left_extended']:
                ext += "L"
            if row['right_extended']:
                ext += "R"
            print(f"{row['chromosome']:>5}  {row['cen_start']:>12,}  {row['cen_end']:>12,}  "
                  f"{row['left_edge']:>12,}  {row['right_edge']:>12,}  "
                  f"{row['n_snps_search_region']:>7}  "
                  f"{row['n_snps_left']:>6}  {row['n_snps_right']:>6}  "
                  f"{row['n_snps_total']:>6}  {row['mean_maf']:>6.4f}  {ext}")


if __name__ == "__main__":
    main()
