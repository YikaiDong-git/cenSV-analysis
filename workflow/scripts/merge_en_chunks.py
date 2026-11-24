#!/usr/bin/env python3
"""Merge chunked EN training outputs into a single trained_models.pkl + sv_metrics.tsv."""
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import joblib

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chunks-dir', required=True, help="Directory containing chunk_*/")
    parser.add_argument('--output', required=True, help="Output directory for merged results")
    args = parser.parse_args()

    chunks_dir = Path(args.chunks_dir)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    chunk_dirs = sorted(chunks_dir.glob("chunk_*"))
    if not chunk_dirs:
        logger.error(f"No chunk directories found in {chunks_dir}")
        sys.exit(1)

    logger.info(f"Found {len(chunk_dirs)} chunks in {chunks_dir}")

    all_models = {}
    all_metrics = []
    scaler = None
    snp_variants = None
    sv_variants = None
    shared_samples = None

    for cd in chunk_dirs:
        pkl_f = cd / "trained_models.pkl"
        tsv_f = cd / "sv_metrics.tsv"
        if not pkl_f.exists() or not tsv_f.exists():
            logger.warning(f"  SKIP {cd.name}: missing files")
            continue

        data = joblib.load(pkl_f)
        all_models.update(data['models'])

        try:
            df = pd.read_csv(tsv_f, sep='\t')
            if len(df) == 0:
                logger.warning(f"  SKIP {cd.name}: empty metrics (header only)")
                continue
            all_metrics.append(df)
        except Exception as e:
            logger.warning(f"  SKIP {cd.name}: cannot read metrics: {e}")
            continue
        logger.info(f"  {cd.name}: {len(data['models'])} models, {len(df)} metrics rows")

        if scaler is None:
            scaler = data['scaler']
            snp_variants = data['snp_variants']
            sv_variants = data['sv_variants']
            shared_samples = data['shared_samples']

    merged_metrics = pd.concat(all_metrics, ignore_index=True)
    merged_metrics = merged_metrics.sort_values('sv_index').reset_index(drop=True)

    logger.info(f"Merged: {len(all_models)} models, {len(merged_metrics)} metrics rows")

    merged_metrics.to_csv(out_dir / "sv_metrics.tsv", sep='\t', index=False)

    save_data = {
        'models': all_models,
        'scaler': scaler,
        'snp_variants': snp_variants,
        'sv_variants': sv_variants,
        'shared_samples': shared_samples,
    }
    joblib.dump(save_data, out_dir / "trained_models.pkl")

    import numpy as np
    coefs = []
    for idx, model in all_models.items():
        row = {'sv_index': idx, 'intercept': model.intercept_}
        for i in np.where(model.coef_ != 0)[0]:
            row[f'snp_{i}'] = model.coef_[i]
        coefs.append(row)
    if coefs:
        pd.DataFrame(coefs).to_csv(out_dir / "model_coefficients.tsv", sep='\t', index=False)

    logger.info(f"Saved merged results to {out_dir}")
    logger.info(f"  trained_models.pkl: {len(all_models)} models")
    logger.info(f"  sv_metrics.tsv: {len(merged_metrics)} rows")


if __name__ == "__main__":
    main()
