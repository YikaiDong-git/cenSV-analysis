#!/usr/bin/env python3
"""elastic_net_train.py — Train elastic net models: peri-centromeric SNPs → centromeric SVs."""
import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut, KFold
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy.stats import pearsonr
import joblib

import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def load_matrix(npz_path, sample_path):
    gt = np.load(npz_path)['genotypes']
    samples = [l.strip() for l in Path(sample_path).read_text().splitlines() if l.strip()]
    return gt, samples

def align_samples(x_matrix, x_samples, y_matrix, y_samples):
    shared = sorted(set(x_samples) & set(y_samples))
    x_map = {s: i for i, s in enumerate(x_samples)}
    y_map = {s: i for i, s in enumerate(y_samples)}
    return x_matrix[:, [x_map[s] for s in shared]], \
           y_matrix[:, [y_map[s] for s in shared]], \
           shared

def train_one_sv(X, y, sv_row, snp_variants, n_folds, n_jobs):
    cv = LeaveOneOut() if n_folds == 0 else KFold(n_splits=min(n_folds, len(y)), shuffle=True, random_state=42)
    model = ElasticNetCV(l1_ratio=[0.3, 0.5, 1.0], cv=cv, alphas=50,
                         max_iter=2000, n_jobs=n_jobs, random_state=42)
    model.fit(X, y)
    y_p = np.clip(model.predict(X), 0, 2)
    mse = mean_squared_error(y, y_p)
    p_r, p_p = pearsonr(y, y_p) if np.std(y_p) > 1e-9 and np.std(y) > 1e-9 else (0.0, 1.0)
    best_cv_mse = np.min(np.mean(model.mse_path_, axis=-1))

    coef = model.coef_
    nonzero_indices = np.where(coef != 0)[0]
    nonzero_coefs = coef[nonzero_indices]
    selected_snps_pos = snp_variants.iloc[nonzero_indices]['pos'].tolist()
    selected_snps_id = snp_variants.iloc[nonzero_indices]['id'].tolist()
    sv_pos = sv_row['pos']
    sv_snp_dist = [abs(p - sv_pos) for p in selected_snps_pos]
    sv_freq = np.sum(y) / (2 * len(y))
    sv_end = sv_row['end']
    abs_coefs = np.abs(coef)
    sorted_indices = np.argsort(abs_coefs)[::-1]
    max_5_coef = coef[sorted_indices[:5]].tolist() if len(coef) > 0 else []

    return {
        'model': model,
        'r2': r2_score(y, y_p),
        'cv_r2': 1.0 - (best_cv_mse / np.var(y)) if np.var(y) > 1e-9 else 0.0,
        'max_coef': np.max(np.abs(model.coef_)) if np.any(model.coef_) else 0.0,
        'max_5_coef': max_5_coef,
        'n_nonzero_coef': len(nonzero_coefs),
        'gt_counts': "/".join([f"{int(i)}:{int(c)}" for i, c in zip(*np.unique(y, return_counts=True))]),
        'l1_ratio': model.l1_ratio_,
        'svfreq': sv_freq,
        'svpos': sv_pos,
        'svend': sv_end,
        'svtype': sv_row['svtype'],
        'svlen': sv_row['svlen'],
        'alpha': model.alpha_,
        'snp_pos': selected_snps_pos,
        'mse': mse,
        'rmse': np.sqrt(mse),
        'cv_mse': best_cv_mse,
        'cv_rmse': np.sqrt(best_cv_mse),
        'pearson_r': p_r,
        'pearson_p': p_p,
        'mae': mean_absolute_error(y, y_p),
        'selected_sv_snp_dist': sv_snp_dist,
        'selected_snp_coefs': nonzero_coefs.tolist(),
        'selected_snp_pos': selected_snps_pos,
        'selected_snp_id': selected_snps_id,
        'sv_id': sv_row['id'],
        'sv_category': sv_row.get('sv_category', 'standard'),
        'original_sv_id': sv_row.get('original_sv_id', ''),
        'mapped_haplogroup': sv_row.get('mapped_haplogroup', ''),
        'r2_score': sv_row.get('r2_score', None),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--snp-dir', required=True)
    parser.add_argument('--sv-dir', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--chromosome', required=True)
    parser.add_argument('--n-folds', type=int, default=5)
    parser.add_argument('--n-jobs', type=int, default=4)
    parser.add_argument('--min-snp-freq', type=float, default=0.1)
    parser.add_argument('--min-sv-freq', type=float, default=0.05)
    parser.add_argument('--max-snp-freq', type=float, default=0.99)
    parser.add_argument('--max-sv-freq', type=float, default=0.99)
    parser.add_argument('--sv-filter', default='all',
                        choices=['all', 'decomposed_only', 'decomposed_and_other', 'full'],
                        help="Filter SVs by category (requires sv_category column)")
    parser.add_argument('--sv-start', type=int, default=None,
                        help="First SV index (0-based, after filtering). Default: 0")
    parser.add_argument('--sv-end', type=int, default=None,
                        help="Last SV index (exclusive, after filtering). Default: all")
    args = parser.parse_args()

    chrom = args.chromosome
    snp_dir = Path(args.snp_dir)
    sv_dir = Path(args.sv_dir)

    logger.info(f"  Loading SNP genotypes...")
    X_raw, x_samples = load_matrix(snp_dir / f"{chrom}_snp_genotypes.npz", snp_dir / f"{chrom}_snp_samples.txt")
    snp_variants = pd.read_csv(snp_dir / f"{chrom}_snp_variants.tsv", sep='\t')
    logger.info(f"  X: {X_raw.shape[0]} SNPs x {X_raw.shape[1]} samples")

    logger.info(f"  Loading SV genotypes...")
    Y_raw, y_samples = load_matrix(sv_dir / f"{chrom}_sv_genotypes.npz", sv_dir / f"{chrom}_sv_samples.txt")
    sv_variants = pd.read_csv(sv_dir / f"{chrom}_sv_variants.tsv", sep='\t')
    logger.info(f"  Y: {Y_raw.shape[0]} SVs x {Y_raw.shape[1]} samples")

    if args.sv_filter != 'all' and 'sv_category' in sv_variants.columns:
        FILTER_MAP = {
            'decomposed_only': {'monophyletic', 'residual'},
            'decomposed_and_other': {'monophyletic', 'residual', 'non_decomposed'},
            'full': {'monophyletic', 'residual', 'non_decomposed', 'original_decomposed'},
        }
        mask = sv_variants['sv_category'].isin(FILTER_MAP[args.sv_filter]).values
        Y_raw = Y_raw[mask]
        sv_variants = sv_variants[mask].reset_index(drop=True)
        logger.info(f"  SV filter '{args.sv_filter}': {mask.sum()} of {len(mask)} SVs kept")

    X_aligned, Y_aligned, shared_samples = align_samples(X_raw, x_samples, Y_raw, y_samples)
    logger.info(f"  Shared samples: {len(shared_samples)}")

    def filter_freq(mat, var, min_f, max_f, label):
        freqs = np.mean((mat == 1) | (mat == 2), axis=1)
        mask = (freqs >= min_f) & (freqs <= max_f)
        logger.info(f"  Filtered {label}: {len(mask)} -> {np.sum(mask)} (freq range: [{min_f}, {max_f}])")
        return mat[mask], var[mask].reset_index(drop=True)

    X_aligned, snp_variants = filter_freq(X_aligned, snp_variants, args.min_snp_freq, args.max_snp_freq, "SNPs")
    Y_aligned, sv_variants = filter_freq(Y_aligned, sv_variants, args.min_sv_freq, args.max_sv_freq, "SVs")

    X = X_aligned.T.astype(np.float64)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    n_total = Y_aligned.shape[0]
    sv_start = args.sv_start if args.sv_start is not None else 0
    sv_end = args.sv_end if args.sv_end is not None else n_total
    sv_start = max(0, sv_start)
    sv_end = min(n_total, sv_end)
    logger.info(f"  Training SV[{sv_start}:{sv_end}) of {n_total} total (n_jobs={args.n_jobs})")

    models, metrics_list = {}, []
    for sv_idx in range(sv_start, sv_end):
        y = Y_aligned[sv_idx, :].astype(np.float64)
        sv_row = sv_variants.iloc[sv_idx]
        res = train_one_sv(X_scaled, y, sv_row, snp_variants, args.n_folds, n_jobs=args.n_jobs)
        models[sv_idx] = res.pop('model')
        res['sv_index'] = sv_idx
        res['chrom'] = chrom
        metrics_list.append(res)
        logger.info(
            f"r2={res['r2']:.3f}, cv_r2={res['cv_r2']:.3f}, "
            f"n_nonzero_coef={res['n_nonzero_coef']}, svpos={res['svpos']}, "
            f"sv_index={sv_idx}/{n_total}"
        )

    out_path = Path(args.output) / chrom
    out_path.mkdir(parents=True, exist_ok=True)
    metrics_df = pd.DataFrame(metrics_list)
    metrics_df.to_csv(out_path / "sv_metrics.tsv", sep='\t', index=False)

    save_data = {
        'models': models,
        'scaler': scaler,
        'snp_variants': snp_variants,
        'sv_variants': sv_variants,
        'shared_samples': shared_samples,
    }
    for key in ['snp_variants', 'sv_variants']:
        df = save_data[key]
        str_cols = df.select_dtypes(include=['string']).columns
        if len(str_cols) > 0:
            save_data[key] = df.astype({c: 'object' for c in str_cols})
    joblib.dump(save_data, out_path / "trained_models.pkl")

    coefs = [
        dict(sv_index=idx, intercept=model.intercept_,
             **{f"snp_{i}": model.coef_[i] for i in np.where(model.coef_ != 0)[0]})
        for idx, model in models.items()
    ]
    if coefs:
        pd.DataFrame(coefs).to_csv(out_path / "model_coefficients.tsv", sep='\t', index=False)

    n_good = len(metrics_df[metrics_df['cv_r2'] > 0]) if len(metrics_df) > 0 else 0
    median_cv_r2 = metrics_df['cv_r2'].median() if len(metrics_df) > 0 else 0
    logger.info(f"  {chrom} done: {len(metrics_list)} SVs trained, "
                f"{n_good} with cv_r2>0, median cv_r2={median_cv_r2:.4f}")
    logger.info(f"  Saved to {out_path}")


if __name__ == "__main__":
    main()
