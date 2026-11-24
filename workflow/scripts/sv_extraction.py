#!/usr/bin/env python3
import argparse
import logging
import re
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

def get_sample_mapping(vcf_path):
    reader = cyvcf2.VCF(vcf_path)
    samples = list(reader.samples)
    reader.close()
    mapping = {}
    for idx, s in enumerate(samples):
        match = re.match(r'^([A-Za-z0-9]+)_(hap[12]|mat|pat)$', s)
        mapping.setdefault(match.group(1), {})[0 if match.group(2) in ['hap1', 'mat'] else 1] = idx
            
    return mapping

def parse_vcf_svs(vcf_path, chrom, overlap_ids, min_maf, max_missing):
    sample_mapping = get_sample_mapping(vcf_path)
    samples = [s for s in overlap_ids if s in sample_mapping]
            
    if not samples:
        return None, None, []
        
    n_samples = len(samples)
    logger.info(f"  Loaded {n_samples} diploid samples from VCF")

    reader = cyvcf2.VCF(str(vcf_path), gts012=True)
    gt_list = []
    variants_info = []

    for v in reader(chrom):        
        raw_gt = v.gt_types 
        diploid_gt = np.full(n_samples, -1, dtype=np.int8)

        for i, sid in enumerate(samples):
            idx1 = sample_mapping[sid][0]
            idx2 = sample_mapping[sid][1]
            
            diploid_gt[i] = (1 if raw_gt[idx1] == 2 else 0) + (1 if raw_gt[idx2] == 2 else 0)
                
        known_mask = diploid_gt >= 0
        n_known = np.sum(known_mask)
        if n_known == 0:
            continue
            
        frac_missing = 1.0 - (n_known / n_samples)
        if frac_missing > max_missing:
            continue
            
        known_gt = diploid_gt[known_mask]
        ac = np.sum(known_gt)
        an = 2 * n_known
        af = ac / an
        maf = min(af, 1.0 - af)

        if maf < min_maf:
            continue

        gt_list.append(diploid_gt)
        variants_info.append({
            'chrom': v.CHROM,
            'pos': v.POS,
            'end': v.INFO.get('END', v.POS),
            'id': v.ID or f"{v.CHROM}_{v.POS}",
            'svtype': v.INFO.get('SVTYPE', 'UNK'),
            'svlen': abs(v.INFO.get('SVLEN', 0) or 0),
            'ref': v.REF,
            'alt': ",".join(v.ALT) if v.ALT else '.',
            'ac': int(ac),
            'an': int(an),
            'af': round(af, 6),
            'maf': round(maf, 6),
            'n_miss': int(n_samples - n_known)
        })

    reader.close()

    if not gt_list:
        logger.warning(f"No SVs found in {chrom}")
        return None, None, samples

    matrix = np.vstack(gt_list)
    return matrix, variants_info, samples

def save_outputs(matrix, variants_info, samples, chrom, output_dir):
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    
    df_var = pd.DataFrame(variants_info)
    np.savez_compressed(out_path / f"{chrom}_sv_genotypes.npz", genotypes=matrix)
    df_var.to_csv(out_path / f"{chrom}_sv_variants.tsv", sep='\t', index=False)
    
    with open(out_path / f"{chrom}_sv_samples.txt", 'w') as fh:
        fh.write('\n'.join(samples) + '\n')
        
    return df_var

def main():
    parser = argparse.ArgumentParser(description="Extract centromeric SV genotypes")
    parser.add_argument('--vcf', required=True, type=Path, help="Cohort VCF file")
    parser.add_argument('--samples', required=True, type=Path, help="Overlap sample IDs")
    parser.add_argument('--output', required=True, type=Path, help="Output directory")
    parser.add_argument('--chromosomes', default='all', help="Chromosomes list")
    parser.add_argument('--min-maf', type=float, default=0.01, help="Min MAF")
    parser.add_argument('--max-missing', type=float, default=0.05, help="Max missing fraction")

    args = parser.parse_args()

    if args.chromosomes.lower() == 'all':
        chroms = [f"chr{c}" for c in range(1, 23)] + ["chrX"]
    else:
        chroms = [c.strip() if c.strip().startswith('chr') else f"chr{c.strip()}" for c in args.chromosomes.split(',')]

    sample_ids = load_sample_list(args.samples)
    logger.info(f"Loaded {len(sample_ids)} target sample IDs")

    summary = []
    for chrom in chroms:
        chrom_clean = chrom.strip()
        matrix, variants, samples = parse_vcf_svs(args.vcf, chrom_clean, sample_ids, args.min_maf, args.max_missing)
        
        if matrix is None:
            continue
            
        df_var = save_outputs(matrix, variants, samples, chrom_clean, args.output)
        
        sv_counts = df_var['svtype'].value_counts().to_dict()
        summary.append({
            'chromosome': chrom_clean,
            'n_svs': matrix.shape[0],
            'n_samples': matrix.shape[1],
            'mean_maf': df_var['maf'].mean(),
            'sv_types': str(sv_counts)
        })
        logger.info(f"Finished {chrom_clean}: {matrix.shape[0]} SVs")

    if summary:
        out_root = Path(args.output)
        pd.DataFrame(summary).to_csv(out_root / "extraction_summary.tsv", sep='\t', index=False)
        logger.info(f"Summary saved to {out_root}")

if __name__ == "__main__":
    main()
