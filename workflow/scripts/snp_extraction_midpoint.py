#!/usr/bin/env python3
"""Centromere-midpoint-outward SNP extraction."""
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
    df['midpoint'] = ((df['start'] + df['end']) / 2).astype(int)
    return df.set_index('chrom')


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


def select_midpoint_outward(variants_info, genotype_matrix, midpoint, n_per_side):
    positions = np.array([v['pos'] for v in variants_info])
    left_mask = positions < midpoint
    right_mask = positions >= midpoint

    left_indices = np.where(left_mask)[0]
    left_distances = midpoint - positions[left_indices]
    left_order = np.argsort(left_distances)
    left_selected = left_indices[left_order[:n_per_side]]

    right_indices = np.where(right_mask)[0]
    right_distances = positions[right_indices] - midpoint
    right_order = np.argsort(right_distances)
    right_selected = right_indices[right_order[:n_per_side]]

    n_left = len(left_selected)
    n_right = len(right_selected)

    if n_left < n_per_side:
        logger.warning(f"  Left flank: only {n_left}/{n_per_side} SNPs available")
    if n_right < n_per_side:
        logger.warning(f"  Right flank: only {n_right}/{n_per_side} SNPs available")

    selected = np.sort(np.concatenate([left_selected, right_selected]))

    filtered_variants = [variants_info[i] for i in selected]
    filtered_matrix = genotype_matrix[selected]

    logger.info(f"  Selected {n_left} left + {n_right} right = {len(selected)} total SNPs")
    if n_left > 0:
        left_span = (midpoint - positions[left_selected].min()) / 1e6
        logger.info(f"  Left span: {left_span:.2f} Mb from midpoint")
    if n_right > 0:
        right_span = (positions[right_selected].max() - midpoint) / 1e6
        logger.info(f"  Right span: {right_span:.2f} Mb from midpoint")

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
        description="Extract peri-centromeric SNPs using midpoint-outward selection"
    )
    parser.add_argument('--vcf-dir', required=True, type=Path, help="VCF directory")
    parser.add_argument('--vcf-pattern', required=True, help="VCF filename pattern with {chrom}")
    parser.add_argument('--samples', required=True, type=Path, help="Sample IDs file")
    parser.add_argument('--centromere-bed', required=True, type=Path, help="Centromere BED file")
    parser.add_argument('--output', required=True, type=Path, help="Output directory")
    parser.add_argument('--n-snps-per-side', type=int, default=2000,
                        help="Number of SNPs to select per flank (default: 2000)")
    parser.add_argument('--max-search-distance', type=int, default=20_000_000,
                        help="Max distance from midpoint to search for SNPs (default: 20Mb)")
    parser.add_argument('--chromosomes', default='all',
                        help="Comma-separated chromosome list or 'all'")
    parser.add_argument('--min-maf', type=float, default=0.01, help="Min MAF (default: 0.01)")
    parser.add_argument('--max-missing', type=float, default=0.05,
                        help="Max missing fraction (default: 0.05)")

    args = parser.parse_args()

    if args.chromosomes.lower() == 'all':
        chroms = [str(c) for c in range(1, 23)] + ['X']
    else:
        chroms = [c.strip().replace('chr', '') for c in args.chromosomes.split(',')]

    sample_ids = load_sample_list(args.samples)
    logger.info(f"Loaded {len(sample_ids)} target sample IDs")

    centromeres = load_centromeres(args.centromere_bed)
    logger.info(f"Loaded centromere coordinates for {len(centromeres)} chromosomes")

    summary = []
    for chrom_num in chroms:
        chrom = f"chr{chrom_num}"
        logger.info(f"\nProcessing {chrom}...")

        if chrom not in centromeres.index:
            logger.warning(f"  {chrom} not found in centromere BED, skipping")
            continue

        midpoint = centromeres.loc[chrom, 'midpoint']
        cen_start = centromeres.loc[chrom, 'start']
        cen_end = centromeres.loc[chrom, 'end']
        logger.info(f"  Centromere: {cen_start:,}-{cen_end:,}  midpoint: {midpoint:,}")

        search_start = max(0, midpoint - args.max_search_distance)
        search_end = midpoint + args.max_search_distance
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

        sel_variants, sel_matrix, n_left, n_right = select_midpoint_outward(
            variants, matrix, midpoint, args.n_snps_per_side
        )

        df_var = save_outputs(sel_matrix, sel_variants, samples, chrom, args.output)

        summary.append({
            'chromosome': chrom,
            'cen_start': cen_start,
            'cen_end': cen_end,
            'midpoint': midpoint,
            'n_snps_search_region': matrix.shape[0],
            'n_snps_left': n_left,
            'n_snps_right': n_right,
            'n_snps_total': n_left + n_right,
            'n_samples': sel_matrix.shape[1],
            'mean_maf': df_var['maf'].mean() if df_var is not None else 0,
        })
        logger.info(f"  Saved {chrom}: {n_left + n_right} SNPs")

    if summary:
        df_summary = pd.DataFrame(summary)
        df_summary.to_csv(args.output / "extraction_summary.tsv", sep='\t', index=False)
        logger.info(f"\nSummary saved to {args.output / 'extraction_summary.tsv'}")

        print("\n" + "=" * 90)
        print("EXTRACTION SUMMARY")
        print("=" * 90)
        print(f"{'Chr':>5}  {'Cen Midpoint':>14}  {'Search':>7}  {'Left':>6}  {'Right':>6}  {'Total':>6}  {'MAF':>6}")
        for _, row in df_summary.iterrows():
            print(f"{row['chromosome']:>5}  {row['midpoint']:>14,}  "
                  f"{row['n_snps_search_region']:>7}  "
                  f"{row['n_snps_left']:>6}  {row['n_snps_right']:>6}  "
                  f"{row['n_snps_total']:>6}  {row['mean_maf']:>6.4f}")


if __name__ == "__main__":
    main()
