#!/usr/bin/env python3
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def load_sample_list(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith('#')]


def main():
    parser = argparse.ArgumentParser(description="Select imputation target samples from 1KG")
    parser.add_argument('--metadata', required=True, type=Path)
    parser.add_argument('--training-samples', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--n-samples', type=int, default=0,
                        help="Number of samples to select (0 = all available)")
    parser.add_argument('--seed', type=int, default=42)

    args = parser.parse_args()

    training_ids = set(load_sample_list(args.training_samples))
    logger.info(f"Loaded {len(training_ids)} training samples to exclude")

    meta = pd.read_csv(args.metadata, sep=r'\s+')
    logger.info(f"Loaded {len(meta)} entries from 1KG metadata")

    training_families = set(meta[meta['SampleID'].isin(training_ids)]['FamilyID'].unique())
    candidates = meta[~meta['FamilyID'].isin(training_families)].copy()
    logger.info(f"Excluded {len(meta) - len(candidates)} entries (training + family members)")

    candidates = candidates[
        (candidates['FatherID'].astype(str).isin(['0', '0.0', 'nan', ''])) &
        (candidates['MotherID'].astype(str).isin(['0', '0.0', 'nan', '']))
    ]
    logger.info(f"After removing children: {len(candidates)} candidates")

    superpops = candidates['Superpopulation'].unique()

    if args.n_samples == 0 or args.n_samples >= len(candidates):
        result = candidates
        for sp in sorted(superpops):
            n_sp = (candidates['Superpopulation'] == sp).sum()
            logger.info(f"  {sp}: {n_sp} candidates (all)")
    else:
        n_per_pop = max(1, args.n_samples // len(superpops))
        remainder = args.n_samples - (n_per_pop * len(superpops))

        selected = []
        for sp in sorted(superpops):
            sp_cands = candidates[candidates['Superpopulation'] == sp]
            n_select = n_per_pop
            if remainder > 0 and len(sp_cands) > n_per_pop:
                n_select += 1
                remainder -= 1
            sample = sp_cands.sample(n=min(n_select, len(sp_cands)), random_state=args.seed)
            selected.append(sample)
            logger.info(f"  {sp}: selected {len(sample)}/{len(sp_cands)} candidates")
        result = pd.concat(selected)

    sample_ids = result['SampleID'].tolist()
    logger.info(f"Selected {len(sample_ids)} samples across {len(superpops)} super-populations")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as fh:
        fh.write('\n'.join(sample_ids) + '\n')

    report_path = args.output.with_suffix('.report.tsv')
    result[['SampleID', 'Population', 'Superpopulation']].to_csv(report_path, sep='\t', index=False)
    logger.info(f"Saved to {args.output} and {report_path}")


if __name__ == "__main__":
    main()
