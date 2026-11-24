#!/usr/bin/env python3

import argparse
import logging
import sys
import warnings
from pathlib import Path
import time
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut, KFold
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy.stats import pearsonr
import joblib


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("scripts/benchmark/results/confidence_check.log")
    ]
)
logger = logging.getLogger(__name__)

STATS = Counter()


def load_matrix(npz_path, sample_path):
    if not Path(npz_path).exists():
        STATS['error_file_not_found'] += 1
        return None, None
        
    data = np.load(npz_path)
    gt = data['genotypes']

    samples = []
    empty_lines = 0
    with open(sample_path) as fh:
        for line in fh:
            s = line.strip()
            if s:
                samples.append(s)
            else:
                empty_lines += 1
    
    if empty_lines > 0:
        STATS['input_empty_sample_lines'] += empty_lines
        logger.info(f"  [CHECK HIT] Found {empty_lines} empty lines in {sample_path}")

    return gt, samples


def align_samples(x_matrix, x_samples, y_matrix, y_samples):
    x_set = set(x_samples)
    y_set = set(y_samples)
    shared = sorted(x_set & y_set)

    if not shared:
        STATS['align_no_shared_samples'] += 1
        logger.info("  [CHECK HIT] No shared samples found between X and Y")
        return None, None, []

    x_idx = [x_samples.index(s) for s in shared]
    y_idx = [y_samples.index(s) for s in shared]

    X_aligned = x_matrix[:, x_idx]
    Y_aligned = y_matrix[:, y_idx]

    return X_aligned, Y_aligned, shared



def train_one_sv(X, y, sv_id, sv_row, snp_variants, n_folds, n_jobs):
    valid_mask = y >= 0
    X_valid = X[valid_mask]
    y_valid = y[valid_mask]
    call_rate = np.mean(valid_mask)
    n_valid = np.sum(valid_mask)
    y_std = np.std(y_valid)
    monomorphic = y_std == 0

    sv_info = {
        'freq': sv_row.get('maf', None),
        'pos': sv_row.get('pos', None),
        'end': sv_row.get('end', None),
        'type': sv_row.get('svtype', None),
        'len': sv_row.get('svlen', None)
    }

    if call_rate < 0.8:
        STATS['skip_low_call_rate'] += 1
        logger.info(
            f"SV {sv_id}: SKIPPED | call_rate={call_rate:.3f}, n_samples={n_valid}, std={y_std:.3f}, freq={sv_info.get('freq')}, pos={sv_info.get('pos')}, end={sv_info.get('end')}, type={sv_info.get('type')}, len={sv_info.get('len')}, reason=call rate low (<0.8)"
        )
        return None

    if n_valid < 10:
        STATS['skip_too_few_samples'] += 1
        logger.info(
            f"SV {sv_id}: SKIPPED | call_rate={call_rate:.3f}, n_samples={n_valid}, std={y_std:.3f}, freq={sv_info.get('freq')}, pos={sv_info.get('pos')}, end={sv_info.get('end')}, type={sv_info.get('type')}, len={sv_info.get('len')}, reason=too few samples (<10)"
        )
        return None

    if monomorphic:
        STATS['skip_monomorphic'] += 1
        logger.info(
            f"SV {sv_id}: SKIPPED | call_rate={call_rate:.3f}, n_samples={n_valid}, std={y_std:.3f}, freq={sv_info.get('freq')}, pos={sv_info.get('pos')}, end={sv_info.get('end')}, type={sv_info.get('type')}, len={sv_info.get('len')}, reason=monomorphic (std=0)"
        )
        return None

    STATS['pass_call_rate'] += int(call_rate >= 0.8)
    STATS['pass_sample_count'] += int(n_valid >= 10)
    STATS['pass_monomorphic'] += int(not monomorphic)

    if n_folds == 0:
        STATS['info_loo_cv_selected'] += 1
        cv = LeaveOneOut()
    else:
        actual_folds = min(n_folds, n_valid)
        if actual_folds < n_folds:
            STATS['info_folds_reduced'] += 1
        cv = KFold(n_splits=actual_folds, shuffle=True, random_state=42)

    l1_ratios = [0.3, 0.5, 1.0]
    model = ElasticNetCV(
        l1_ratio=l1_ratios,
        alphas=50,
        cv=cv,
        max_iter=2000,
        tol=1e-4,
        n_jobs=n_jobs,
        random_state=42
    )

    try:
        model.fit(X_valid, y_valid)
    except Exception as e:
        STATS['error_training_exception'] += 1
        logger.info(
            f"SV {sv_id}: SKIPPED | call_rate={call_rate:.3f}, n_samples={n_valid}, std={np.std(y_valid):.3f}, freq={sv_info.get('freq')}, pos={sv_info.get('pos')}, end={sv_info.get('end')}, type={sv_info.get('type')}, len={sv_info.get('len')}, reason=training failed ({e})"
        )
        return None

    y_pred = model.predict(X_valid)
    y_pred_clipped = np.clip(y_pred, 0, 2)

    r2 = r2_score(y_valid, y_pred_clipped)
    mse = mean_squared_error(y_valid, y_pred_clipped)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_valid, y_pred_clipped)

    if np.std(y_valid) > 0 and np.std(y_pred_clipped) > 0:
        pearson_r, pearson_p = pearsonr(y_valid, y_pred_clipped)
    else:
        STATS['warn_pearson_undefined'] += 1
        pearson_r, pearson_p = 0.0, 1.0

    n_nonzero_coef = np.sum(model.coef_ != 0)
    max_coef = np.max(np.abs(model.coef_)) if n_nonzero_coef > 0 else 0.0

    mean_cv_mse_path = np.mean(model.mse_path_, axis=-1)
    best_cv_mse = np.min(mean_cv_mse_path)
    cv_rmse = np.sqrt(best_cv_mse)
    
    y_var = np.var(y_valid)
    cv_r2 = 1.0 - (best_cv_mse / y_var) if y_var > 1e-9 else 0.0

    unique, counts = np.unique(y_valid, return_counts=True)
    gt_counts = dict(zip(unique, counts))
    gt_counts_str = f"0:{gt_counts.get(0.0,0)}/1:{gt_counts.get(1.0,0)}/2:{gt_counts.get(2.0,0)}"

    STATS['success_trained'] += 1
    selected_snps = np.where(model.coef_ != 0)[0]
    selected_coefs = model.coef_[selected_snps]
    
    if len(selected_coefs) > 0:
        abs_coefs = np.abs(selected_coefs)
        top_5_idxs = np.argsort(abs_coefs)[-5:][::-1]
        max_5_coef = selected_coefs[top_5_idxs].tolist()
        max_coef = float(np.max(abs_coefs))
        mean_coef = float(np.mean(selected_coefs))
        std_coef = float(np.std(selected_coefs))
        median_coef = float(np.median(selected_coefs))
        
        selected_snp_positions = snp_variants.iloc[selected_snps]['pos'].values.tolist()
        sv_pos = sv_info.get('pos')
        if sv_pos is not None:
            sv_snp_distances = [int(abs(p - sv_pos)) for p in selected_snp_positions]
        else:
            sv_snp_distances = []
    else:
        max_5_coef = []
        max_coef = 0.0
        mean_coef = 0.0
        std_coef = 0.0
        median_coef = 0.0
        selected_snp_positions = []
        sv_snp_distances = []

    logger.info(
        f"SV: TRAINED | "
        f"r2={r2:.3f}, cv_r2={cv_r2:.3f}, max_coef={max_coef:.3f}, max_5_coef={[round(c, 4) for c in max_5_coef]}, "
        f"n_nonzero_coef={n_nonzero_coef}, GT={gt_counts_str}, l1_ratio={model.l1_ratio_}, "
        f"mean_coef={mean_coef:.4f}, std_coef={std_coef:.4f}, median_coef={median_coef:.4f} | "
        f"n_snps={X_valid.shape[1]}, call_rate={call_rate:.3f}, std={y_std:.3f}, "
        f"freq={sv_info.get('freq')}, pos={sv_info.get('pos')}, end={sv_info.get('end')}, type={sv_info.get('type')}, len={sv_info.get('len')}, "
        f"mse={mse:.3f}, rmse={rmse:.3f}, mae={mae:.3f}, alpha={model.alpha_:.4e}, "
        f"snp_pos={selected_snp_positions}, sv_snp_dist={sv_snp_distances}, "
        f"selected_snps={selected_snps.tolist()}, selected_coefs={[round(c, 4) for c in selected_coefs.tolist()]} | "
        f"sv_id={sv_id}"
    )
    return {
        'sv_id': sv_id,
        'r2': r2,
        'cv_r2': cv_r2,
        'max_coef': max_coef,
        'max_5_coef': max_5_coef,
        'n_nonzero_coef': n_nonzero_coef,
        'GT': gt_counts_str,
        'l1_ratio': model.l1_ratio_,
        'mean_coef': mean_coef,
        'std_coef': std_coef,
        'median_coef': median_coef,
        'n_samples': n_valid,
        'n_snps': X_valid.shape[1],
        'call_rate': call_rate,
        'std': float(y_std),
        'freq': sv_info.get('freq'),
        'pos': sv_info.get('pos'),
        'end': sv_info.get('end'),
        'type': sv_info.get('type'),
        'len': sv_info.get('len'),
        'mse': mse,
        'rmse': rmse,
        'mae': mae,
        'alpha': model.alpha_,
        'snp_pos': selected_snp_positions,
        'sv_snp_dist': sv_snp_distances,
        'selected_snps': selected_snps.tolist(),
        'selected_coefs': selected_coefs.tolist()
    }


def train_chromosome(snp_dir, sv_dir, chrom, n_folds, n_jobs, min_snp_freq=0.0, min_sv_freq=0.0, max_svs=None):
    snp_dir = Path(snp_dir)
    sv_dir = Path(sv_dir)

    logger.info(f"  Checking inputs for {chrom}...")
    snp_npz = snp_dir / f"{chrom}_snp_genotypes.npz"
    if not snp_npz.exists():
        STATS['error_snp_file_missing'] += 1
        logger.warning(f"  [CHECK HIT] SNP matrix not found for {chrom}: {snp_npz}")
        return None

    sv_npz = sv_dir / f"{chrom}_sv_genotypes.npz"
    if not sv_npz.exists():
        STATS['error_sv_file_missing'] += 1
        logger.warning(f"  [CHECK HIT] SV matrix not found for {chrom}: {sv_npz}")
        return None

    X_raw, x_samples = load_matrix(snp_npz, snp_dir / f"{chrom}_snp_samples.txt")
    snp_variants = pd.read_csv(snp_dir / f"{chrom}_snp_variants.tsv", sep='\t')
    
    Y_raw, y_samples = load_matrix(sv_npz, sv_dir / f"{chrom}_sv_samples.txt")
    sv_variants = pd.read_csv(sv_dir / f"{chrom}_sv_variants.tsv", sep='\t')

    X_aligned, Y_aligned, shared_samples = align_samples(X_raw, x_samples, Y_raw, y_samples)

    if X_aligned is None or len(shared_samples) < 10:
        STATS['error_insufficient_shared_samples'] += 1
        logger.warning(f"  [CHECK HIT] Too few shared samples for {chrom}: {len(shared_samples) if shared_samples else 0}")
        return None

    if min_snp_freq > 0:
        presence = (X_aligned == 1) | (X_aligned == 2)
        freqs = np.mean(presence, axis=1)
        mask = freqs >= min_snp_freq
        n_before = X_aligned.shape[0]
        n_after = np.sum(mask)
        n_filtered = n_before - n_after
        X_aligned = X_aligned[mask]
        snp_variants = snp_variants[mask].reset_index(drop=True)
        STATS['info_snps_filtered_freq'] += n_filtered
        logger.info(f"Filtered SNPs: {n_before} -> {n_after} (filtered: {n_filtered}, min_freq={min_snp_freq})")

    if min_sv_freq > 0:
        presence = (Y_aligned == 1) | (Y_aligned == 2)
        freqs = np.mean(presence, axis=1)
        mask = freqs >= min_sv_freq
        n_before = Y_aligned.shape[0]
        n_after = np.sum(mask)
        n_filtered = n_before - n_after
        Y_aligned = Y_aligned[mask]
        sv_variants = sv_variants[mask].reset_index(drop=True)
        STATS['info_svs_filtered_freq'] += n_filtered
        logger.info(f"Filtered SVs: {n_before} -> {n_after} (filtered: {n_filtered}, min_freq={min_sv_freq})")

    X = X_aligned.T.astype(np.float64)
    Y = Y_aligned

    n_missing = np.sum(X == 3)
    if n_missing > 0:
        STATS['info_X_values_imputed'] += n_missing
    else:
        STATS['info_X_no_missing_values'] += 1

    X[X == 3] = np.nan
    col_means = np.nanmean(X, axis=0)
    col_means = np.where(np.isnan(col_means), 0, col_means)
    nan_mask = np.isnan(X)
    X[nan_mask] = np.take(col_means, np.where(nan_mask)[1])

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    n_svs = Y.shape[0]
    models = {}
    metrics_list = []
    
    if max_svs and n_svs > max_svs:
        logger.info(f"  Limiting analysis to first {max_svs} SVs for benchmarking purposes.")
        n_to_process = max_svs
    else:
        n_to_process = n_svs

    logger.info(f"  Training {n_to_process} SV models (out of {n_svs} total)...")

    for sv_idx in range(n_to_process):
        y = Y[sv_idx, :].astype(np.float64)
        sv_row = sv_variants.iloc[sv_idx]
        sv_id = sv_row.get('id', f"SV_{sv_idx}")

        result = train_one_sv(X_scaled, y, sv_id, sv_row, snp_variants, n_folds, n_jobs=n_jobs)

        if result is None:
            continue

        models[sv_idx] = None
        metrics_list.append(result)

    if not metrics_list:
        STATS['warn_no_models_trained'] += 1
        logger.warning(f"  [CHECK HIT] No SV models trained for {chrom}")
        return None

    df = pd.DataFrame(metrics_list)
    out_path = Path(sv_dir) / f"{chrom}_sv_model_metrics.tsv"
    df.to_csv(out_path, sep='\t', index=False)
    logger.info(f"Wrote per-SV metrics to {out_path}")
    return {'metrics': df, 'chrom': chrom}

    if not metrics_list:
        STATS['warn_no_models_trained'] += 1
        logger.warning(f"  [CHECK HIT] No SV models trained for {chrom}")
        return None

    return {'metrics': pd.DataFrame(metrics_list), 'chrom': chrom}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--snp-dir', required=True)
    parser.add_argument('--sv-dir', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--chromosomes', default='all')
    parser.add_argument('--n-folds', type=int, default=2)
    parser.add_argument('--n-jobs', type=int, default=4)
    parser.add_argument('--min-snp-freq', type=float, default=0.1)
    parser.add_argument('--min-sv-freq', type=float, default=0.05)
    parser.add_argument('--max-snp-freq', type=float, default=0.99)
    parser.add_argument('--max-sv-freq', type=float, default=0.99)
    parser.add_argument('--max-svs', type=int, default=None,
                        help="Limit training to first N SVs for quick testing")
    args = parser.parse_args()

    Path(args.output).mkdir(parents=True, exist_ok=True)
    
    if args.chromosomes.lower() == 'all':
        chroms = [f"chr{c}" for c in range(1, 23)]
    else:
        chroms = [c.strip() if c.strip().startswith('chr') else f"chr{c.strip()}"
                  for c in args.chromosomes.split(',')]

    for chrom in chroms:
        train_chromosome(
            args.snp_dir, args.sv_dir, chrom, args.n_folds, args.n_jobs,
            min_snp_freq=args.min_snp_freq, min_sv_freq=args.min_sv_freq,
            max_svs=args.max_svs
        )

    report_path = Path(args.output) / "confidence_report.txt"
    with open(report_path, 'w') as f:
        f.write("ELASTIC NET CONFIDENCE CHECK REPORT\n")
        f.write("===================================\n\n")
        f.write("Check Hit Counts:\n")
        for key, value in sorted(STATS.items()):
            f.write(f"{key:<30}: {value}\n")

        f.write("\nFilter Summary:\n")
        f.write(f"- SNPs filtered by min freq: {STATS.get('info_snps_filtered_freq', 0)}\n")
        f.write(f"- SVs filtered by min freq: {STATS.get('info_svs_filtered_freq', 0)}\n")
        f.write(f"- SVs skipped by call rate: {STATS.get('skip_low_call_rate', 0)}\n")
        f.write(f"- SVs skipped by sample count: {STATS.get('skip_too_few_samples', 0)}\n")
        f.write(f"- SVs skipped by monomorphic: {STATS.get('skip_monomorphic', 0)}\n")
        f.write(f"- SVs skipped by training exception: {STATS.get('error_training_exception', 0)}\n")
        f.write(f"- SVs skipped by alignment (no shared samples): {STATS.get('align_no_shared_samples', 0)}\n")
        f.write(f"- SVs skipped by file not found: {STATS.get('error_file_not_found', 0) + STATS.get('error_snp_file_missing', 0) + STATS.get('error_sv_file_missing', 0)}\n")
        f.write(f"- SNPs imputed: {STATS.get('info_X_values_imputed', 0)}\n")
        f.write(f"- SNPs with no missing values: {STATS.get('info_X_no_missing_values', 0)}\n")

        f.write("\nAnalysis:\n")
        if STATS['skip_monomorphic'] > 0:
            f.write("- 'skip_monomorphic': HIT. Monomorphic checks are necessary.\n")
        else:
            f.write("- 'skip_monomorphic': NEVER HIT. Could potentially remove this check for this dataset.\n")

        if STATS['skip_low_call_rate'] > 0:
             f.write("- 'skip_low_call_rate': HIT. Missing data handling is active.\n")

        if STATS['error_file_not_found'] > 0 or STATS['error_snp_file_missing'] > 0:
             f.write("- File existence checks: HIT. Critical for pipeline stability.\n")

        f.write("\nConclusion:\n")
        f.write("Any count > 0 implies that the corresponding 'if' statement guarded against a real occurrence.\n")

    logger.info("\nREPORT GENERATED:\n")
    with open(report_path, 'r') as f:
        report_content = f.read()
        print(report_content)
        logger.info(report_content)

if __name__ == "__main__":
    main()
