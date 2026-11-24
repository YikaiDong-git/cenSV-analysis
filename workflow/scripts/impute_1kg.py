#!/usr/bin/env python3
"""Impute centromeric SV dosages in new samples using trained elastic net models."""
import argparse
import logging
import sys
from pathlib import Path

import cyvcf2
import joblib
import numpy as np
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


def load_trained_models(model_dir):
    model_dir = Path(model_dir)
    try:
        trained = joblib.load(model_dir / "trained_models.pkl")
    except NotImplementedError:
        import pickle as _pickle
        import joblib.numpy_pickle as _jnp

        class _CompatUnpickler(_jnp.NumpyUnpickler):
            dispatch = _jnp.NumpyUnpickler.dispatch.copy()

            def load_build(self):
                stack = self.stack
                state = stack[-1]
                inst = stack[-2] if len(stack) >= 2 else None
                if (isinstance(inst, pd.arrays.StringArray)
                        and isinstance(state, tuple) and len(state) == 2):
                    _dtype_obj, arr = state
                    stack.pop()
                    stack[-1] = pd.array(np.asarray(arr, dtype=object),
                                         dtype='string')
                    return
                _jnp.NumpyUnpickler.load_build(self)

            dispatch[_pickle.BUILD[0]] = load_build

        pkl_path = model_dir / "trained_models.pkl"
        with open(pkl_path, 'rb') as f:
            trained = _CompatUnpickler(
                str(pkl_path), f, ensure_native_byte_order=True
            ).load()

    metrics = pd.read_csv(model_dir / "sv_metrics.tsv", sep='\t')
    return trained, metrics


def extract_snps_for_targets(vcf_path, chrom, snp_variants, target_samples):
    snp_lookup = {}
    for i, row in snp_variants.iterrows():
        key = (str(row['chrom']), int(row['pos']), str(row['ref']), str(row['alt']))
        snp_lookup[key] = i

    n_snps = len(snp_variants)
    logger.info(f"  Looking for {n_snps} training SNP positions in VCF")

    reader = cyvcf2.VCF(str(vcf_path), samples=target_samples, gts012=True)
    vcf_samples = list(reader.samples)
    found_samples = [s for s in target_samples if s in vcf_samples]
    logger.info(f"  Found {len(found_samples)}/{len(target_samples)} target samples in VCF")

    vcf_sample_idx = {s: i for i, s in enumerate(vcf_samples)}
    target_indices = [vcf_sample_idx[s] for s in found_samples]

    X_new = np.zeros((len(found_samples), n_snps), dtype=np.float64)
    matched_count = 0

    positions = snp_variants['pos'].values
    min_pos = int(positions.min()) - 1000
    max_pos = int(positions.max()) + 1000
    region = f"{chrom}:{min_pos}-{max_pos}"
    logger.info(f"  Scanning VCF region: {region}")

    for v in reader(region):
        if len(v.ALT) != 1 or len(v.REF) != 1 or len(v.ALT[0]) != 1:
            continue

        key = (v.CHROM, v.POS, v.REF, v.ALT[0])
        if key not in snp_lookup:
            continue

        col_idx = snp_lookup[key]
        gt = v.gt_types

        for j, vidx in enumerate(target_indices):
            g = gt[vidx]
            X_new[j, col_idx] = g if g != 3 else 0

        matched_count += 1

    reader.close()
    logger.info(f"  Matched {matched_count}/{n_snps} training SNPs in VCF "
                f"({100*matched_count/n_snps:.1f}%)")

    return X_new, matched_count, found_samples


def impute_sv_dosages(X_new, trained, metrics, r2_threshold=0.0):
    scaler = trained['scaler']
    models = trained['models']
    sv_variants = trained['sv_variants']

    X_scaled = scaler.transform(X_new)

    predictions = {}
    sv_info = []
    n_skipped_r2 = 0

    for sv_idx, model in models.items():
        cv_r2 = metrics.iloc[sv_idx].get('cv_r2', 0.0) if sv_idx < len(metrics) else 0.0

        if cv_r2 < r2_threshold:
            n_skipped_r2 += 1
            continue

        pred = np.clip(model.predict(X_scaled), 0, 2)
        predictions[sv_idx] = pred

        info = {'sv_index': sv_idx, 'cv_r2': cv_r2}
        row = sv_variants.iloc[sv_idx]
        for col in ['id', 'chrom', 'pos', 'end', 'svtype', 'svlen', 'sv_category',
                    'original_sv_id', 'mapped_haplogroup']:
            if col in row.index:
                info[col] = row[col]
        mrow = metrics.iloc[sv_idx]
        for col in ['r2', 'pearson_r', 'svfreq', 'svpos', 'svend']:
            if col in mrow.index:
                info[col] = mrow[col]

        info['imputed_mean'] = round(float(np.mean(pred)), 4)
        info['imputed_std'] = round(float(np.std(pred)), 4)
        sv_info.append(info)

    logger.info(f"  Imputed {len(predictions)} SVs "
                f"(skipped {n_skipped_r2} below cv_r2 threshold {r2_threshold})")

    return predictions, sv_info


def save_imputation_results(predictions, sv_info, samples, chrom, output_dir,
                            matched_snps, total_snps):
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    sv_indices = sorted(predictions.keys())
    dosage_data = []
    for sv_idx in sv_indices:
        row = {'sv_index': sv_idx}
        for j, sid in enumerate(samples):
            row[sid] = round(float(predictions[sv_idx][j]), 4)
        dosage_data.append(row)

    pd.DataFrame(dosage_data).to_csv(out_path / f"{chrom}_imputed_dosages.tsv", sep='\t', index=False)
    pd.DataFrame(sv_info).to_csv(out_path / f"{chrom}_imputation_summary.tsv", sep='\t', index=False)

    with open(out_path / f"{chrom}_snp_matching.tsv", 'w') as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"training_snps\t{total_snps}\n")
        fh.write(f"matched_snps\t{matched_snps}\n")
        fh.write(f"match_rate\t{matched_snps/total_snps:.4f}\n")
        fh.write(f"n_target_samples\t{len(samples)}\n")

    with open(out_path / f"{chrom}_imputed_samples.txt", 'w') as fh:
        fh.write('\n'.join(samples) + '\n')

    logger.info(f"  Saved results to {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Impute centromeric SV dosages in new samples using trained elastic net models"
    )
    parser.add_argument('--model-dir', required=True, type=Path)
    parser.add_argument('--vcf', required=True, type=Path)
    parser.add_argument('--samples', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--chromosome', required=True)
    parser.add_argument('--r2-threshold', type=float, default=0.0)

    args = parser.parse_args()

    trained, metrics = load_trained_models(args.model_dir)
    training_samples = set(trained['shared_samples'])
    snp_variants = trained['snp_variants']
    logger.info(f"Training used {len(training_samples)} samples, "
                f"{len(snp_variants)} SNPs, {len(trained['sv_variants'])} SVs")

    all_targets = load_sample_list(args.samples)
    target_samples = [s for s in all_targets if s not in training_samples]
    logger.info(f"Target samples: {len(target_samples)}")

    X_new, matched_snps, found_samples = extract_snps_for_targets(
        args.vcf, args.chromosome, snp_variants, target_samples
    )

    total_snps = len(snp_variants)
    predictions, sv_info = impute_sv_dosages(X_new, trained, metrics, args.r2_threshold)

    if not predictions:
        logger.error("No SVs passed the R² threshold")
        sys.exit(1)

    save_imputation_results(
        predictions, sv_info, found_samples, args.chromosome, args.output,
        matched_snps, total_snps
    )
    logger.info(f"Done. Imputed {len(predictions)} SVs for {len(found_samples)} samples.")


if __name__ == "__main__":
    main()
