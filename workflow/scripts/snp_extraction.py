#!/usr/bin/env python3
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

def parse_vcf_snps(vcf_path, chrom, sample_ids, min_maf, max_missing):
    reader = cyvcf2.VCF(str(vcf_path), gts012=True, samples=sample_ids)
    samples = list(reader.samples)
    n_samples = len(samples)
    logger.info(f"  Loaded {n_samples} samples from VCF")

    gt_list = []
    variants_info = []

    for v in reader(chrom):
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
        logger.warning(f"No variants found in {chrom}")
        return None, None, samples

    matrix = np.vstack(gt_list)
    return matrix, variants_info, samples

def save_outputs(matrix, variants_info, samples, chrom, output_dir):
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    df_var = pd.DataFrame(variants_info)
    np.savez_compressed(out_path / f"{chrom}_snp_genotypes.npz", genotypes=matrix)
    df_var.to_csv(out_path / f"{chrom}_snp_variants.tsv", sep='\t', index=False)
    with open(out_path / f"{chrom}_snp_samples.txt", 'w') as fh:
        fh.write('\n'.join(samples) + '\n')
    return df_var

def main():
    parser = argparse.ArgumentParser(description="Extract peri-centromeric SNPs")
    parser.add_argument('--vcf-dir', required=True, type=Path, help="VCF directory")
    parser.add_argument('--vcf-pattern', required=True, help="VCF filename pattern")
    parser.add_argument('--samples', required=True, type=Path, help="Overlap sample IDs")
    parser.add_argument('--output', required=True, type=Path, help="Output directory")
    parser.add_argument('--chromosomes', default='all', help="Chromosomes list")
    parser.add_argument('--min-maf', type=float, default=0.01, help="Min MAF")
    parser.add_argument('--max-missing', type=float, default=0.05, help="Max missing fraction")

    args = parser.parse_args()
    
    if args.chromosomes.lower() == 'all':
        chroms = [str(c) for c in range(1, 23)] + ['X']
    else:
        chroms = [c.strip().replace('chr', '') for c in args.chromosomes.split(',')]

    sample_ids = load_sample_list(args.samples)
    logger.info(f"Loaded {len(sample_ids)} target sample IDs")

    summary = []
    for chrom_num in chroms:
        chrom = f"chr{chrom_num}"

        vcf_path = args.vcf_dir / args.vcf_pattern.format(chrom=chrom_num)
        matrix, variants, samples = parse_vcf_snps(vcf_path, chrom, sample_ids, args.min_maf, args.max_missing)
        df_var = save_outputs(matrix, variants, samples, chrom, args.output)
        summary.append({
            'chromosome': chrom, 'n_snps': matrix.shape[0], 
            'n_samples': matrix.shape[1], 'mean_maf': df_var['maf'].mean()
        })
        logger.info(f"Finished {chrom}: {matrix.shape[0]} SNPs")

    if summary:
        pd.DataFrame(summary).to_csv(args.output / "extraction_summary.tsv", sep='\t', index=False)
        logger.info(f"Summary saved to {args.output}")

if __name__ == "__main__":
    main()
