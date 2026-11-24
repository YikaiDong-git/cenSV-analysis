#!/usr/bin/env python3
"""Unified SV-QTL association pipeline using tensorQTL."""
import argparse
import gzip
import logging
import json
import os
import time
import urllib.request
import sys
from pathlib import Path

import matplotlib
from adjustText import adjust_text
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

HG38_CHROM_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
    'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
    'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
    'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
    'chr22': 50818468, 'chrX': 156040895,
}



def load_sample_list(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith('#')]


def load_imputation_results(imputation_dir, chrom):
    imp_dir = Path(imputation_dir)
    df_dosage = pd.read_csv(imp_dir / f"{chrom}_imputed_dosages.tsv", sep='\t')
    df_summary = pd.read_csv(imp_dir / f"{chrom}_imputation_summary.tsv", sep='\t')
    samples = load_sample_list(imp_dir / f"{chrom}_imputed_samples.txt")
    logger.info(f"Loaded {len(df_dosage)} imputed SVs for {len(samples)} samples")
    return df_dosage, df_summary


def load_centromere_coords(bed_path, chrom):
    df = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end'])
    row = df[df['chrom'] == chrom].iloc[0]
    start, end = int(row['start']), int(row['end'])
    midpoint = (start + end) // 2
    logger.info(f"Centromere {chrom}: {start:,}-{end:,} (midpoint={midpoint:,})")
    return start, end, midpoint


def load_expression(expr_path, chrom):
    logger.info(f"Loading expression from {expr_path}" +
                (f" ({chrom})" if chrom else " (ALL chromosomes)"))
    df = pd.read_csv(expr_path, sep='\t', compression='gzip')

    if chrom is not None:
        chrom_str = chrom.replace('chr', '')
        mask = (df['#Chr'].astype(str) == chrom) | (df['#Chr'].astype(str) == chrom_str)
        df = df[mask].copy()
        logger.info(f"  {len(df)} genes on {chrom}")

    meta_cols = ['#Chr', 'start', 'end', 'ID']
    sample_cols = [c for c in df.columns if c not in meta_cols]

    phenotype_pos_df = pd.DataFrame({
        'chr': df['#Chr'].astype(str).apply(
            lambda x: x if x.startswith('chr') else f'chr{x}').values,
        'start': df['start'].astype(int).values,
        'end': df['end'].astype(int).values,
    }, index=df['ID'].values)

    phenotype_df = df[sample_cols].copy()
    phenotype_df.index = df['ID'].values
    phenotype_df = phenotype_df.astype(float)

    if chrom is None:
        logger.info(f"  {len(phenotype_df)} genes across all chromosomes")

    return phenotype_df, phenotype_pos_df


def load_splicing_phenotypes(splicing_bed_path, groups_path, chrom):
    logger.info(f"Loading splicing phenotypes from {splicing_bed_path}" +
                (f" ({chrom})" if chrom else " (ALL chromosomes)"))
    df = pd.read_csv(splicing_bed_path, sep='\t', compression='gzip')

    if chrom is not None:
        chrom_str = chrom.replace('chr', '')
        mask = (df['#chr'].astype(str) == chrom) | (df['#chr'].astype(str) == chrom_str)
        df = df[mask].copy()
        logger.info(f"  {len(df)} intron phenotypes on {chrom}")

    meta_cols = ['#chr', 'start', 'end', 'ID']
    sample_cols = [c for c in df.columns if c not in meta_cols]

    phenotype_pos_df = pd.DataFrame({
        'chr': df['#chr'].astype(str).apply(
            lambda x: x if x.startswith('chr') else f'chr{x}').values,
        'start': df['start'].astype(int).values,
        'end': df['end'].astype(int).values,
    }, index=df['ID'].values)

    phenotype_df = df[sample_cols].copy()
    phenotype_df.index = df['ID'].values
    phenotype_df = phenotype_df.astype(float)

    logger.info(f"Loading phenotype groups from {groups_path}")
    groups_df = pd.read_csv(groups_path, sep='\t', header=None,
                            names=['intron_id', 'gene_id'])
    groups_series = groups_df.set_index('intron_id')['gene_id']
    groups_series = groups_series.loc[groups_series.index.isin(phenotype_df.index)]

    logger.info(f"  {len(groups_series)} introns in {groups_series.nunique()} genes")

    return phenotype_df, phenotype_pos_df, groups_series


def load_covariates(cov_path, pc_sex=True, pc_geno=5, pc_peer=60):
    import re
    logger.info(f"Loading covariates from {cov_path}")
    df = pd.read_csv(cov_path, sep='\t', compression='gzip', index_col=0)
    df.loc['sex'] = df.loc['sex'].map({'XY': 1, 'XX': 0, 'male': 1, 'female': 0}).fillna(0)

    geno_rows = [idx for idx in df.index if re.match(r'^PC\d+$', str(idx))]
    peer_rows = [idx for idx in df.index if re.match(r'^PEER\d+$', str(idx))]

    drop_rows = []
    if not pc_sex:
        drop_rows.append('sex')
    drop_rows.extend(geno_rows[pc_geno:])
    drop_rows.extend(peer_rows[pc_peer:])

    if drop_rows:
        df = df.drop([r for r in drop_rows if r in df.index])

    kept_geno = len(geno_rows) - len(geno_rows[pc_geno:])
    kept_peer = len(peer_rows) - len(peer_rows[pc_peer:])
    logger.info(f"  Covariates: sex={'yes' if pc_sex else 'no'}, "
                f"geno_PCs={kept_geno}/{len(geno_rows)}, "
                f"PEER={kept_peer}/{len(peer_rows)}")

    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
    cov_df = df.T.astype(np.float32)
    logger.info(f"  {cov_df.shape[1]} covariates for {cov_df.shape[0]} samples")
    return cov_df



def reconstruct_scaler_and_predict(snp_dir, sv_dir, model_dir, chrom,
                                   imp_summary, hprc_target_samples,
                                   min_snp_freq=0.10, max_snp_freq=0.99,
                                   min_sv_freq=0.05, max_sv_freq=0.99):
    from sklearn.preprocessing import StandardScaler

    snp_dir = Path(snp_dir)
    sv_dir = Path(sv_dir)
    model_dir = Path(model_dir)

    X_raw = np.load(snp_dir / f"{chrom}_snp_genotypes.npz")['genotypes']
    x_samples = [l.strip() for l in open(snp_dir / f"{chrom}_snp_samples.txt") if l.strip()]
    snp_variants = pd.read_csv(snp_dir / f"{chrom}_snp_variants.tsv", sep='\t')
    logger.info(f"  SNPs loaded: {X_raw.shape[0]} x {X_raw.shape[1]} samples")

    Y_raw = np.load(sv_dir / f"{chrom}_sv_genotypes.npz")['genotypes']
    y_samples = [l.strip() for l in open(sv_dir / f"{chrom}_sv_samples.txt") if l.strip()]
    logger.info(f"  SVs loaded: {Y_raw.shape[0]} x {Y_raw.shape[1]} samples")

    shared = sorted(set(x_samples) & set(y_samples))
    x_map = {s: i for i, s in enumerate(x_samples)}
    y_map = {s: i for i, s in enumerate(y_samples)}
    X_aligned = X_raw[:, [x_map[s] for s in shared]]
    logger.info(f"  Shared training samples: {len(shared)}")

    def filter_freq(mat, var, min_f, max_f, label):
        freqs = np.mean((mat == 1) | (mat == 2), axis=1)
        mask = (freqs >= min_f) & (freqs <= max_f)
        logger.info(f"  Filtered {label}: {len(mask)} -> {np.sum(mask)}")
        return mat[mask], var[mask].reset_index(drop=True), mask

    X_aligned, _, snp_mask = filter_freq(
        X_aligned, snp_variants, min_snp_freq, max_snp_freq, "SNPs")

    X_train = X_aligned.T.astype(np.float64)
    scaler = StandardScaler()
    scaler.fit(X_train)
    logger.info(f"  Scaler fitted on {X_train.shape[0]} samples x {X_train.shape[1]} SNPs")

    target_in_snp = [s for s in hprc_target_samples if s in x_map]
    logger.info(f"  HPRC targets found in SNP data: {len(target_in_snp)}/{len(hprc_target_samples)}")

    if not target_in_snp:
        return {}, []

    X_hprc_raw = X_raw[:, [x_map[s] for s in target_in_snp]]
    X_hprc_filtered = X_hprc_raw[snp_mask]
    X_hprc = X_hprc_filtered.T.astype(np.float64)
    X_hprc_scaled = scaler.transform(X_hprc)
    logger.info(f"  HPRC SNP matrix: {X_hprc_scaled.shape}")

    coef_df = pd.read_csv(model_dir / "model_coefficients.tsv", sep='\t')
    logger.info(f"  Model coefficients: {len(coef_df)} SVs, up to {coef_df.shape[1]-2} SNP features")

    n_snps = X_hprc_scaled.shape[1]
    predictions = {}

    for _, imp_row in imp_summary.iterrows():
        sv_idx = int(imp_row['sv_index'])
        coef_row = coef_df[coef_df['sv_index'] == sv_idx]
        if len(coef_row) == 0:
            predictions[sv_idx] = np.full(len(target_in_snp), np.nan)
            continue

        coef_row = coef_row.iloc[0]
        intercept = coef_row['intercept']

        beta = np.zeros(n_snps)
        for col in coef_row.index:
            if col.startswith('snp_') and pd.notna(coef_row[col]):
                snp_i = int(col.replace('snp_', ''))
                if snp_i < n_snps:
                    beta[snp_i] = coef_row[col]

        dosage = np.clip(X_hprc_scaled @ beta + intercept, 0, 2)
        predictions[sv_idx] = dosage

    logger.info(f"  Predicted dosages for {len(predictions)} SVs x {len(target_in_snp)} HPRC samples")
    return predictions, target_in_snp


def augment_dosages_with_hprc(dosage_df, imp_summary, phenotype_path,
                              snp_dir, sv_dir, model_dir, chrom):
    imp_sample_cols = [c for c in dosage_df.columns if c != 'sv_index']

    with gzip.open(str(phenotype_path), 'rt') as f:
        header = f.readline().strip().split('\t')
    pheno_samples = set(header[4:])

    sv_samples = [l.strip() for l in open(Path(sv_dir) / f"{chrom}_sv_samples.txt") if l.strip()]
    hprc_with_pheno = sorted(set(sv_samples) & pheno_samples - set(imp_sample_cols))

    if not hprc_with_pheno:
        logger.info("No HPRC samples with phenotype data found; skipping augmentation")
        return dosage_df

    logger.info(f"HPRC samples with phenotype data (to add): {len(hprc_with_pheno)}")
    logger.info(f"  Samples: {hprc_with_pheno}")

    predictions, found_samples = reconstruct_scaler_and_predict(
        snp_dir=snp_dir, sv_dir=sv_dir, model_dir=model_dir,
        chrom=chrom, imp_summary=imp_summary,
        hprc_target_samples=hprc_with_pheno,
    )

    if not found_samples:
        logger.info("No HPRC samples found in SNP data; skipping augmentation")
        return dosage_df

    pred_records = []
    for _, row in dosage_df.iterrows():
        sv_idx = int(row['sv_index'])
        if sv_idx in predictions:
            pred = predictions[sv_idx]
            pred_records.append({s: round(float(pred[j]), 4)
                                 for j, s in enumerate(found_samples)})
        else:
            pred_records.append({s: np.nan for s in found_samples})

    pred_df = pd.DataFrame(pred_records, index=dosage_df.index)
    combined = pd.concat([dosage_df, pred_df], axis=1)

    all_sample_cols = imp_sample_cols + found_samples
    logger.info(f"Augmented dosage matrix: {len(combined)} SVs x {len(all_sample_cols)} samples")
    return combined



def prepare_genotype_data(dosage_df, sv_summary, centromere_midpoint, chrom,
                          min_r2, min_maf):
    valid = sv_summary.dropna(subset=['cv_r2']).copy()
    valid = valid[valid['cv_r2'] >= min_r2]
    logger.info(f"  {len(valid)} SVs pass cv_r2 >= {min_r2}")

    sample_cols = [c for c in dosage_df.columns if c != 'sv_index']
    merged = valid.merge(dosage_df, on='sv_index', how='inner')

    dosage_vals = merged[sample_cols].values.astype(float)
    af = dosage_vals.mean(axis=1) / 2
    maf = np.minimum(af, 1 - af)
    merged = merged[maf >= min_maf].copy()
    logger.info(f"  {len(merged)} SVs pass MAF >= {min_maf}")

    variant_ids = [f"sv_{int(row['sv_index'])}_{row.get('svtype', 'SV')}"
                   for _, row in merged.iterrows()]

    genotype_df = merged[sample_cols].copy()
    genotype_df.index = variant_ids
    genotype_df = genotype_df.astype(np.float32)

    variant_df = pd.DataFrame({
        'chrom': chrom,
        'pos': centromere_midpoint,
    }, index=variant_ids)
    logger.info(f"  All SVs assigned to hg38 centromere midpoint ({centromere_midpoint:,})")

    sv_r2 = {vid: row['cv_r2'] for vid, (_, row) in zip(variant_ids, merged.iterrows())}

    return genotype_df, variant_df, sv_r2



def compute_qvalues(cis_df, method='storey'):
    if method == 'storey':
        try:
            from tensorqtl import calculate_qvalues
            calculate_qvalues(cis_df, fdr=0.05)
        except Exception as e:
            logger.warning(f"  Storey q-value failed ({e}), falling back to BH FDR")
            method = 'bh'
    if method == 'bh':
        from statsmodels.stats.multitest import multipletests
        pvals = cis_df['pval_beta'].dropna().values
        if len(pvals) > 0:
            _, fdr, _, _ = multipletests(pvals, method='fdr_bh')
            cis_df['qval'] = np.nan
            cis_df.loc[cis_df['pval_beta'].notna(), 'qval'] = fdr
        else:
            cis_df['qval'] = np.nan
    logger.info(f"  {(cis_df['qval'] <= 0.05).sum()} genes with qval <= 0.05 ({method})")


def run_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
            covariates_df, window, output_dir, chrom, seed, fdr_method,
            qtl_type, group_s=None):
    from tensorqtl import cis

    label = 'sQTL' if qtl_type == 'sqtl' else 'eQTL'
    logger.info(f"Running tensorQTL cis-{label} mapping (window={window/1e6:.0f}Mb)...")
    logger.info(f"  {len(genotype_df)} variants, {len(phenotype_df)} phenotypes, "
                f"{genotype_df.shape[1]} samples")

    cis_kwargs = dict(
        covariates_df=covariates_df,
        window=window,
        seed=seed,
    )
    if group_s is not None:
        cis_kwargs['group_s'] = group_s

    try:
        cis_df = cis.map_cis(genotype_df, variant_df,
                             phenotype_df, phenotype_pos_df, **cis_kwargs)
    except ValueError as e:
        if "No phenotypes remain" in str(e):
            logger.warning("  No phenotypes remain after filters — writing empty outputs")
            out_path = Path(output_dir)
            out_path.mkdir(parents=True, exist_ok=True)
            pd.DataFrame().to_csv(out_path / f"{chrom}_cis_permutation.tsv", sep='\t')
            return None
        raise
    logger.info(f"  Permutation pass: {len(cis_df)} genes tested")

    compute_qvalues(cis_df, method=fdr_method)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    perm_file = out_path / f"{chrom}_cis_permutation.tsv"
    cis_df.to_csv(perm_file, sep='\t')
    logger.info(f"  Saved permutation results: {perm_file}")

    nominal_prefix = f"{chrom}_cis_nominal"
    nominal_kwargs = dict(
        covariates_df=covariates_df,
        prefix=nominal_prefix,
        output_dir=str(out_path),
        window=window,
    )
    if group_s is not None:
        nominal_kwargs['group_s'] = group_s

    cis.map_nominal(genotype_df, variant_df,
                    phenotype_df, phenotype_pos_df, **nominal_kwargs)
    logger.info(f"  Saved nominal results: {out_path}/{nominal_prefix}*.parquet")

    return cis_df


def run_trans(genotype_df, phenotype_df, phenotype_pos_df, covariates_df,
              sv_chrom, sv_r2, pval_threshold, output_dir):
    from tensorqtl import trans

    logger.info(f"Running trans-QTL mapping ({sv_chrom} SVs vs target genes)...")
    logger.info(f"  pval_threshold={pval_threshold}")

    trans_df = trans.map_trans(
        genotype_df, phenotype_df,
        covariates_df=covariates_df,
        pval_threshold=pval_threshold,
        maf_threshold=0.0,
        return_r2=True,
        batch_size=20000,
    )
    logger.info(f"  {len(trans_df)} significant pairs (p < {pval_threshold})")

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    if len(trans_df) == 0:
        logger.info("  No significant trans pairs found.")
        pd.DataFrame(columns=['variant_id', 'phenotype_id', 'pval', 'b', 'b_se']).to_csv(
            out_path / f"{sv_chrom}_trans_results.tsv", sep='\t', index=False)
        return trans_df

    trans_df['sv_chrom'] = sv_chrom
    trans_df['gene_chrom'] = trans_df['phenotype_id'].map(
        lambda g: phenotype_pos_df.loc[g, 'chr'] if g in phenotype_pos_df.index else 'NA')
    trans_df['gene_start'] = trans_df['phenotype_id'].map(
        lambda g: phenotype_pos_df.loc[g, 'start'] if g in phenotype_pos_df.index else 0)
    trans_df['is_cis'] = trans_df['gene_chrom'] == sv_chrom
    trans_df['sv_cv_r2'] = trans_df['variant_id'].map(sv_r2)
    trans_df = trans_df.sort_values('pval')

    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(trans_df['pval'].values, method='fdr_bh')
    trans_df['qval_bh'] = fdr

    n_cis = trans_df['is_cis'].sum()
    n_trans = (~trans_df['is_cis']).sum()
    n_fdr05 = (fdr <= 0.05).sum()
    logger.info(f"  cis pairs: {n_cis}, trans pairs: {n_trans}")
    logger.info(f"  BH-FDR q<=0.05: {n_fdr05}")

    trans_df.to_csv(out_path / f"{sv_chrom}_trans_results.tsv", sep='\t', index=False)
    logger.info(f"  Saved: {out_path / f'{sv_chrom}_trans_results.tsv'}")

    return trans_df



def build_association_table(output_dir, chrom, qtl_type):
    out_path = Path(output_dir)
    parquet_files = sorted(out_path.glob(f"{chrom}_cis_nominal*.parquet"))
    pairs = pd.concat([pd.read_parquet(f) for f in parquet_files], ignore_index=True)
    label = 'SV-intron' if qtl_type == 'sqtl' else 'SV-gene'
    logger.info(f"  {len(pairs)} {label} pairs from nominal pass")
    return pairs



def make_summary_plot(cis_df, pairs_df, sv_summary, chrom, output_dir, qtl_type):
    label = 'sQTL' if qtl_type == 'sqtl' else 'eQTL'
    pair_label = 'SV-intron' if qtl_type == 'sqtl' else 'SV-gene'

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Centromeric SV-{label} association (tensorQTL) -- {chrom}', fontsize=14)

    ax = axes[0, 0]
    pos = sv_summary['pos'].dropna().values / 1e6
    r2 = sv_summary.loc[sv_summary['pos'].notna(), 'cv_r2'].values
    scatter = ax.scatter(pos, r2, c=r2, cmap='viridis', s=15, alpha=0.7)
    ax.set_xlabel('Position (Mb, CHM13)')
    ax.set_ylabel('cv_R2')
    ax.set_title('Imputed SVs: position vs model quality')
    plt.colorbar(scatter, ax=ax, label='cv_R2')

    ax = axes[0, 1]
    pvals = cis_df['pval_beta'].dropna().values
    ax.hist(pvals, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Permutation p-value (beta-approximation)')
    ax.set_ylabel('Count')
    ax.set_title(f'Per-gene empirical p-values (n={len(pvals)})')
    ax.axvline(0.05, color='red', linestyle='--', alpha=0.5, label='p=0.05')
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    logp = -np.log10(pairs_df['pval_nominal'].clip(lower=1e-300).values)
    ax.hist(logp, bins=50, color='gray', alpha=0.7, edgecolor='black')
    ax.axvline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.5, label='p=0.05')
    ax.set_xlabel('-log10(nominal p-value)')
    ax.set_ylabel('Count')
    ax.set_title(f'{pair_label} pair p-values (n={len(pairs_df)})')
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.axis('off')
    n_genes = len(cis_df)
    n_sig_nom = (cis_df['pval_beta'].dropna() < 0.05).sum()
    n_sig_q = (cis_df['qval'].dropna() <= 0.05).sum()
    n_sig_q10 = (cis_df['qval'].dropna() <= 0.10).sum()
    summary_text = (
        f"Genes tested: {n_genes}\n"
        f"Significant (perm p<0.05): {n_sig_nom}\n"
        f"Significant (qval<=0.05): {n_sig_q}\n"
        f"Significant (qval<=0.10): {n_sig_q10}\n\n"
        f"Total {pair_label} pairs: {len(pairs_df)}\n"
        f"Nominal p<0.05: {(pairs_df['pval_nominal'] < 0.05).sum()}"
    )
    ax.text(0.1, 0.5, summary_text, transform=ax.transAxes,
            va='center', fontsize=12, family='monospace')
    ax.set_title('Summary')

    plt.tight_layout()
    plot_path = Path(output_dir) / f"{chrom}_{qtl_type}_summary_plot.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"  Saved summary plot: {plot_path}")




_PROJ_ROOT = Path(CENSV_ROOT + "")
_SHARED_CACHE = _PROJ_ROOT / "resources" / "ensembl_symbol_cache.json"


def _load_symbol_cache():
    if _SHARED_CACHE.exists():
        try:
            return json.loads(_SHARED_CACHE.read_text())
        except Exception:
            pass
    return {}


def _save_symbol_cache(cache):
    try:
        _SHARED_CACHE.parent.mkdir(parents=True, exist_ok=True)
        tmp = _SHARED_CACHE.with_suffix('.json.tmp.' + str(os.getpid()))
        tmp.write_text(json.dumps(cache, indent=0))
        os.replace(tmp, _SHARED_CACHE)
    except Exception as e:
        logger.warning(f"Could not save symbol cache: {e}")


def _query_ensembl_symbol(gene_id):
    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=0;content-type=application/json"
    try:
        req = urllib.request.Request(url)
        req.add_header("Content-Type", "application/json")
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
            return data.get("display_name") or gene_id
    except Exception:
        return None


def resolve_gene_symbols(ensg_set, cache=None):
    if cache is None:
        cache = _load_symbol_cache()
    out = {}
    new_queries = 0
    for g in ensg_set:
        if g in cache and cache[g].get("symbol"):
            out[g] = cache[g]["symbol"]
        else:
            sym = _query_ensembl_symbol(g)
            if sym:
                cache[g] = {"symbol": sym}
                out[g] = sym
                new_queries += 1
                time.sleep(0.1)
            else:
                out[g] = g
    if new_queries:
        _save_symbol_cache(cache)
        logger.info(f"  resolve_gene_symbols: {new_queries} new API hits, cache size now {len(cache)}")
    return out

def make_detail_plot(cis_df, pairs_df, chrom, output_dir, phenotype_pos_df,
                     cen_start, cen_end, qtl_type):
    label = 'sQTL' if qtl_type == 'sqtl' else 'eQTL'
    pair_label = 'SV-intron' if qtl_type == 'sqtl' else 'SV-gene'
    pos_label = 'Intron position' if qtl_type == 'sqtl' else 'Gene TSS position'

    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    fig.suptitle(f'Centromeric SV-{label} details -- {chrom}', fontsize=14)

    ax = axes[0, 0]
    plot_df = cis_df[cis_df['pval_beta'].notna()].copy()
    plot_df = plot_df.join(phenotype_pos_df[['start']], how='left')
    plot_df = plot_df.dropna(subset=['start', 'pval_beta'])

    pos_mb = plot_df['start'].values / 1e6
    logp = -np.log10(plot_df['pval_beta'].clip(lower=1e-300).values)
    qvals = plot_df['qval'].values
    is_sig = qvals <= 0.05
    is_suggestive = (qvals > 0.05) & (qvals <= 0.1)
    is_ns = ~is_sig & ~is_suggestive

    ax.scatter(pos_mb[is_ns], logp[is_ns], c='gray', s=25, alpha=0.6, label='Not significant')
    ax.scatter(pos_mb[is_suggestive], logp[is_suggestive], c='orange', s=40, alpha=0.85,
               edgecolors='black', linewidths=0.5, zorder=4, label='0.05 < q <= 0.1')
    ax.scatter(pos_mb[is_sig], logp[is_sig], c='red', s=50, alpha=0.9,
               edgecolors='black', linewidths=0.5, zorder=5, label='q <= 0.05')

    label_mask = is_sig | is_suggestive
    ensg_set = set()
    for idx in plot_df.index[label_mask]:
        if qtl_type == 'sqtl':
            parts = idx.split(':')
            ensg_set.add(parts[-1].split('.')[0] if len(parts) >= 5 else idx.split('.')[0])
        else:
            ensg_set.add(idx.split('.')[0])
    sym_map = resolve_gene_symbols(ensg_set) if ensg_set else {}
    texts = []
    for idx in plot_df.index[label_mask]:
        if qtl_type == 'sqtl':
            parts = idx.split(':')
            ensg = parts[-1].split('.')[0] if len(parts) >= 5 else idx.split('.')[0]
        else:
            ensg = idx.split('.')[0]
        gene_name = sym_map.get(ensg, ensg)
        x = plot_df.loc[idx, 'start'] / 1e6
        y = -np.log10(max(plot_df.loc[idx, 'pval_beta'], 1e-300))
        texts.append(ax.text(x, y, gene_name, fontsize=7,
                             ha='center', va='bottom',
                             bbox=dict(boxstyle='round,pad=0.15', fc='white',
                                       ec='lightgray', alpha=0.85, lw=0.3)))
    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5, lw=0.5),
                    expand_text=(1.4, 1.5), expand_points=(1.3, 1.5),
                    force_text=(0.7, 1.0), force_points=(0.5, 0.7),
                    lim=400, autoalign='xy')

    ax.axvspan(cen_start / 1e6, cen_end / 1e6, color='orange', alpha=0.2, label='Centromere')
    ax.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.4, linewidth=0.8)
    ax.set_xlabel(f'{pos_label} (Mb, hg38)')
    ax.set_ylabel('-log10(permutation p-value)')
    ax.set_title(f'Manhattan plot (per-gene {label})')
    ax.legend(fontsize=7, loc='upper right')

    ax = axes[0, 1]
    slope = pairs_df['slope'].values
    logp = -np.log10(pairs_df['pval_nominal'].clip(lower=1e-300).values)
    sig_mask = pairs_df['pval_nominal'].values < 0.05

    ax.scatter(slope[~sig_mask], logp[~sig_mask], c='gray', s=10, alpha=0.4)
    ax.scatter(slope[sig_mask], logp[sig_mask], c='red', s=15, alpha=0.7, label='p < 0.05')
    ax.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.4, linewidth=0.8)
    ax.axvline(0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)
    ax.set_xlabel('Effect size (slope)')
    ax.set_ylabel('-log10(nominal p-value)')
    ax.set_title(f'Volcano plot (all {pair_label} pairs)')
    ax.legend(fontsize=7)

    ax = axes[1, 0]
    qvals_arr = np.sort(cis_df['qval'].dropna().values)
    n_q = len(qvals_arr)
    if n_q > 0:
        exp_q = -np.log10((np.arange(1, n_q + 1)) / (n_q + 1))
        obs_q = -np.log10(qvals_arr.clip(min=1e-300))
        ax.scatter(exp_q, obs_q, c='steelblue', s=10, alpha=0.7)
        max_v = max(exp_q.max(), obs_q.max()) * 1.05
        ax.plot([0, max_v], [0, max_v], 'r--', linewidth=1, alpha=0.7, label='Uniform expectation')
        n_q05 = int((qvals_arr <= 0.05).sum())
        n_q10 = int((qvals_arr <= 0.10).sum())
        ax.text(0.05, 0.92, f'n={n_q} genes\nq\u22640.05: {n_q05}\nq\u22640.10: {n_q10}',
                transform=ax.transAxes, fontsize=9, va='top')
    ax.set_xlabel('Expected -log10(uniform)')
    ax.set_ylabel('Observed -log10(qval)')
    ax.set_title(f'QQ plot of qvals (per-gene Storey FDR)')
    ax.legend(fontsize=7)

    ax = axes[1, 1]
    dist_mb = pairs_df['start_distance'].abs().values / 1e6
    logp = -np.log10(pairs_df['pval_nominal'].clip(lower=1e-300).values)
    sig_mask = pairs_df['pval_nominal'].values < 0.05

    ax.scatter(dist_mb[~sig_mask], logp[~sig_mask], c='gray', s=10, alpha=0.4)
    ax.scatter(dist_mb[sig_mask], logp[sig_mask], c='red', s=15, alpha=0.7, label='p < 0.05')
    ax.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.4, linewidth=0.8)
    ax.set_xlabel('Distance from centromere midpoint (Mb)')
    ax.set_ylabel('-log10(nominal p-value)')
    ax.set_title(f'Distance-decay of {label} signal')
    ax.legend(fontsize=7)

    plt.tight_layout()
    plot_path = Path(output_dir) / f"{chrom}_{qtl_type}_detail_plot.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"  Saved detail plot: {plot_path}")


def make_qq_diagnostics(cis_df, pairs_df, chrom, output_dir, qtl_type):
    from scipy import stats as _sp_stats
    label = 'sQTL' if qtl_type == 'sqtl' else 'eQTL'
    pair_label = 'SV-intron' if qtl_type == 'sqtl' else 'SV-gene'

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle(f'QQ diagnostics -- {chrom} {label}', fontsize=14)

    ax = axes[0, 0]
    pbeta = np.sort(cis_df['pval_beta'].dropna().values)
    n_a = len(pbeta)
    if n_a > 0:
        exp_a = -np.log10((np.arange(1, n_a + 1)) / (n_a + 1))
        obs_a = -np.log10(pbeta.clip(min=1e-300))
        ax.scatter(exp_a, obs_a, c='steelblue', s=10, alpha=0.7)
        m = max(exp_a.max(), obs_a.max()) * 1.05
        ax.plot([0, m], [0, m], 'r--', lw=1, alpha=0.7, label='Uniform')
        chi2_a = _sp_stats.chi2.ppf(1 - pbeta.clip(min=1e-300, max=1-1e-15), df=1)
        lam_a = np.median(chi2_a) / _sp_stats.chi2.ppf(0.5, df=1)
        ax.text(0.05, 0.92, f'$\\lambda_{{GC}}$ = {lam_a:.3f}\nn={n_a} genes',
                transform=ax.transAxes, fontsize=10, va='top')
    ax.set_xlabel('Expected -log10(uniform)')
    ax.set_ylabel('Observed -log10(pval_beta)')
    ax.set_title(f'A) Gene-level perm pval QQ')
    ax.legend(fontsize=7)

    ax = axes[0, 1]
    qvals_arr = np.sort(cis_df['qval'].dropna().values)
    n_b = len(qvals_arr)
    if n_b > 0:
        exp_b = -np.log10((np.arange(1, n_b + 1)) / (n_b + 1))
        obs_b = -np.log10(qvals_arr.clip(min=1e-300))
        ax.scatter(exp_b, obs_b, c='darkgreen', s=10, alpha=0.7)
        m = max(exp_b.max(), obs_b.max()) * 1.05
        ax.plot([0, m], [0, m], 'r--', lw=1, alpha=0.7, label='Uniform')
        ax.text(0.05, 0.92, f'q\u22640.05: {(qvals_arr<=0.05).sum()}\nq\u22640.10: {(qvals_arr<=0.10).sum()}',
                transform=ax.transAxes, fontsize=10, va='top')
    ax.set_xlabel('Expected -log10(uniform)')
    ax.set_ylabel('Observed -log10(qval)')
    ax.set_title(f'B) qval QQ (Storey FDR)')
    ax.legend(fontsize=7)

    ax = axes[1, 0]
    if n_a > 0:
        ax.hist(pbeta, bins=40, color='lightsteelblue', edgecolor='gray', alpha=0.85)
        pi0 = float((pbeta > 0.5).sum()) / (0.5 * n_a) if n_a > 0 else 1.0
        pi0 = min(max(pi0, 0.0), 1.0)
        ax.axhline(n_a / 40.0, color='black', linestyle='-', alpha=0.5, lw=0.8,
                   label=f'Uniform expectation (n/40 = {n_a/40:.1f})')
        ax.axhline(n_a * pi0 / 40.0, color='red', linestyle='--', alpha=0.85, lw=1.2,
                   label=f'pi0 = {pi0:.3f} (Storey lambda=0.5)')
    ax.set_xlabel('pval_beta (per-gene)')
    ax.set_ylabel('Count')
    ax.set_title(f'C) p-value histogram + pi0 diagnostic')
    ax.legend(fontsize=7)

    ax = axes[1, 1]
    if pairs_df is not None and 'pval_nominal' in pairs_df.columns:
        pn = np.sort(pairs_df['pval_nominal'].dropna().values)
        n_d = len(pn)
        if n_d > 0:
            exp_d = -np.log10((np.arange(1, n_d + 1)) / (n_d + 1))
            obs_d = -np.log10(pn.clip(min=1e-300))
            if n_d > 50000:
                step = n_d // 50000
                idx = np.arange(0, n_d, step)
                ax.scatter(exp_d[idx], obs_d[idx], c='gray', s=5, alpha=0.4)
            else:
                ax.scatter(exp_d, obs_d, c='gray', s=5, alpha=0.4)
            m = max(exp_d.max(), obs_d.max()) * 1.05
            ax.plot([0, m], [0, m], 'r--', lw=1, alpha=0.7, label='Uniform')
            chi2_d = _sp_stats.chi2.ppf(1 - pn.clip(min=1e-300, max=1-1e-15), df=1)
            lam_d = np.median(chi2_d) / _sp_stats.chi2.ppf(0.5, df=1)
            ax.text(0.05, 0.92, f'$\\lambda_{{GC}}$ = {lam_d:.3f}\nn={n_d} {pair_label} pairs',
                    transform=ax.transAxes, fontsize=10, va='top')
    ax.set_xlabel('Expected -log10(uniform)')
    ax.set_ylabel('Observed -log10(pval_nominal)')
    ax.set_title(f'D) Per-pair pval_nominal QQ (inflation ref)')
    ax.legend(fontsize=7)

    plt.tight_layout()
    plot_path = Path(output_dir) / f"{chrom}_{qtl_type}_qq_diagnostics.png"
    fig.savefig(plot_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"  Saved QQ diagnostics: {plot_path}")


def make_trans_plots(trans_df, phenotype_pos_df, sv_chrom, output_dir,
                     sv_r2, pval_threshold, qtl_type):
    label = 'sQTL' if qtl_type == 'sqtl' else 'eQTL'

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.suptitle(f'Trans-{label}: {sv_chrom} centromeric SVs vs all genes', fontsize=14)

    ax = axes[0, 0]
    chrom_order = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    chrom_to_num = {c: i for i, c in enumerate(chrom_order)}

    plot_df = trans_df[trans_df['gene_chrom'].isin(chrom_order)].copy()
    plot_df['chrom_num'] = plot_df['gene_chrom'].map(chrom_to_num)
    plot_df['logp'] = -np.log10(plot_df['pval'].clip(lower=1e-300))

    colors = []
    for _, row in plot_df.iterrows():
        if row['gene_chrom'] == sv_chrom:
            colors.append('red')
        elif row['chrom_num'] % 2 == 0:
            colors.append('steelblue')
        else:
            colors.append('navy')

    ax.scatter(plot_df['chrom_num'], plot_df['logp'], c=colors, s=15, alpha=0.7)
    ax.set_xticks(range(len(chrom_order)))
    ax.set_xticklabels([c.replace('chr', '') for c in chrom_order], fontsize=6, rotation=45)
    ax.set_xlabel('Gene chromosome')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'Trans-Manhattan ({sv_chrom} SVs, red=cis)')

    ax = axes[0, 1]
    chrom_counts = trans_df.groupby('gene_chrom').size()
    chroms_present = [c for c in chrom_order if c in chrom_counts.index]
    counts = [chrom_counts.get(c, 0) for c in chroms_present]
    bar_colors = ['red' if c == sv_chrom else 'steelblue' for c in chroms_present]
    ax.bar(range(len(chroms_present)), counts, color=bar_colors, alpha=0.7,
           edgecolor='black', linewidth=0.5)
    ax.set_xticks(range(len(chroms_present)))
    ax.set_xticklabels([c.replace('chr', '') for c in chroms_present], fontsize=7, rotation=45)
    ax.set_xlabel('Gene chromosome')
    ax.set_ylabel(f'# pairs (p < {pval_threshold})')
    ax.set_title('Significant pairs per chromosome')

    ax = axes[0, 2]
    cis_pvals = trans_df[trans_df['is_cis']]['pval'].values
    trans_pvals = trans_df[~trans_df['is_cis']]['pval'].values
    if len(cis_pvals) > 0:
        ax.hist(-np.log10(cis_pvals), bins=30, alpha=0.7, color='red',
                label=f'cis (n={len(cis_pvals)})')
    if len(trans_pvals) > 0:
        ax.hist(-np.log10(trans_pvals), bins=30, alpha=0.7, color='steelblue',
                label=f'trans (n={len(trans_pvals)})')
    ax.set_xlabel('-log10(p-value)')
    ax.set_ylabel('Count')
    ax.set_title('p-value distribution: cis vs trans')
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    cis_mask = trans_df['is_cis']
    if cis_mask.any():
        ax.scatter(trans_df.loc[cis_mask, 'b'].abs(),
                   -np.log10(trans_df.loc[cis_mask, 'pval']),
                   c='red', s=20, alpha=0.7, label='cis')
    if (~cis_mask).any():
        ax.scatter(trans_df.loc[~cis_mask, 'b'].abs(),
                   -np.log10(trans_df.loc[~cis_mask, 'pval']),
                   c='steelblue', s=15, alpha=0.6, label='trans')
    ax.set_xlabel('|Effect size|')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Effect size: cis vs trans')
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    trans_only = trans_df[~trans_df['is_cis']]
    if len(trans_only) > 0:
        ax.scatter(trans_only['sv_cv_r2'],
                   -np.log10(trans_only['pval'].clip(lower=1e-300)),
                   c='steelblue', s=20, alpha=0.7)
        ax.set_xlabel('SV cv_R2')
        ax.set_ylabel('-log10(p-value)')
    else:
        ax.text(0.5, 0.5, 'No trans hits', ha='center', va='center', transform=ax.transAxes)
    ax.set_title('Trans hits: SV cv_R2')

    ax = axes[1, 2]
    ax.axis('off')
    n_total = len(trans_df)
    n_cis = trans_df['is_cis'].sum()
    n_trans = (~trans_df['is_cis']).sum()
    n_fdr05 = (trans_df['qval_bh'] <= 0.05).sum() if 'qval_bh' in trans_df.columns else 0
    n_fdr10 = (trans_df['qval_bh'] <= 0.10).sum() if 'qval_bh' in trans_df.columns else 0

    trans_hits = trans_df[~trans_df['is_cis']].head(10)
    top_text = ""
    for _, row in trans_hits.iterrows():
        gene = row['phenotype_id'].split('.')[0]
        top_text += f"  {row['variant_id']:15s} -> {gene} ({row['gene_chrom']}) p={row['pval']:.2e} b={row['b']:.3f}\n"

    summary = (
        f"SV chromosome: {sv_chrom}\n"
        f"Total pairs (p < {pval_threshold}): {n_total}\n"
        f"  cis (same chr): {n_cis}\n"
        f"  trans (diff chr): {n_trans}\n"
        f"BH-FDR q<=0.05: {n_fdr05}\n"
        f"BH-FDR q<=0.10: {n_fdr10}\n\n"
        f"Top trans hits:\n{top_text}"
    )
    ax.text(0.05, 0.95, summary, transform=ax.transAxes, va='top', fontsize=9, family='monospace')
    ax.set_title('Summary')

    plt.tight_layout()
    plot_path = Path(output_dir) / f"{sv_chrom}_trans_{qtl_type}_plot.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"  Saved trans plot: {plot_path}")



def main():
    parser = argparse.ArgumentParser(
        description="Unified SV-QTL association using tensorQTL"
    )
    parser.add_argument('--imputation-dir', required=True, type=Path)
    parser.add_argument('--centromere-bed', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--chromosome', required=True,
                        help='Chromosome of SVs (e.g. chr10)')
    parser.add_argument('--qtl-type', required=True, choices=['eqtl', 'sqtl'])

    parser.add_argument('--expression', type=Path,
                        help='Expression BED (required if eqtl)')
    parser.add_argument('--eqtl-covariates', type=Path,
                        help='eQTL covariates (required if eqtl)')

    parser.add_argument('--splicing-bed', type=Path,
                        help='Splicing BED (required if sqtl)')
    parser.add_argument('--phenotype-groups', type=Path,
                        help='Intron->gene groups (required if sqtl)')
    parser.add_argument('--sqtl-covariates', type=Path,
                        help='sQTL covariates (required if sqtl)')

    parser.add_argument('--regmode', default='within-chrom',
                        choices=['within-chrom', 'cross-chrom'],
                        help='Regression mode: within-chrom or cross-chrom')
    parser.add_argument('--regmode-para', default='5000000',
                        help=('For within-chrom: cis window size in bp '
                              '(0 = whole chromosome). '
                              'For cross-chrom: comma-separated target '
                              'chromosomes (e.g. 2,3,4,6)'))

    parser.add_argument('--min-r2', type=float, default=0.1)
    parser.add_argument('--min-maf', type=float, default=0.05)
    parser.add_argument('--pc-sex', action='store_true', default=True,
                        help='Include sex covariate (default: True)')
    parser.add_argument('--no-pc-sex', dest='pc_sex', action='store_false',
                        help='Exclude sex covariate')
    parser.add_argument('--pc-geno', type=int, default=5,
                        help='Number of genotype PCs to keep (default: 5, 0=none)')
    parser.add_argument('--pc-peer', type=int, default=60,
                        help='Number of PEER factors to keep (default: 60, 0=none)')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--fdr-method', choices=['storey', 'bh'], default='storey',
                        help="Cis FDR method: 'storey' or 'bh'")
    parser.add_argument('--pval-threshold', type=float, default=1e-5,
                        help='p-value threshold for cross-chrom results')

    parser.add_argument('--include-hprc', action='store_true', default=False,
                        help='Include EN-imputed HPRC samples (default: False; opt-in to prevent training leakage)')
    parser.add_argument('--no-hprc', dest='include_hprc', action='store_false',
                        help='Disable HPRC augmentation (default)')
    parser.add_argument('--model-dir', type=Path,
                        help='EN model coefficients dir (required if --include-hprc)')
    parser.add_argument('--snp-dir', type=Path,
                        help='SNP extraction dir (required if --include-hprc)')
    parser.add_argument('--sv-dir', type=Path,
                        help='SV extraction dir (required if --include-hprc)')

    args = parser.parse_args()
    chrom = args.chromosome

    if args.qtl_type == 'eqtl':
        if not args.expression or not args.eqtl_covariates:
            parser.error("--expression and --eqtl-covariates required for eqtl")
    else:
        if not args.splicing_bed or not args.phenotype_groups or not args.sqtl_covariates:
            parser.error("--splicing-bed, --phenotype-groups, and --sqtl-covariates required for sqtl")

    if args.include_hprc and (not args.model_dir or not args.snp_dir or not args.sv_dir):
        parser.error("--model-dir, --snp-dir, and --sv-dir required when --include-hprc is enabled")

    out_path = Path(args.output)
    out_path.mkdir(parents=True, exist_ok=True)

    dosage_df, sv_summary = load_imputation_results(args.imputation_dir, chrom)

    if args.include_hprc:
        phenotype_path = args.expression if args.qtl_type == 'eqtl' else args.splicing_bed
        dosage_df = augment_dosages_with_hprc(
            dosage_df, sv_summary, phenotype_path,
            args.snp_dir, args.sv_dir, args.model_dir, chrom
        )

    cen_start, cen_end, cen_mid = load_centromere_coords(args.centromere_bed, chrom)

    genotype_df, variant_df, sv_r2 = prepare_genotype_data(
        dosage_df, sv_summary, cen_mid, chrom,
        args.min_r2, args.min_maf
    )

    group_s = None
    if args.regmode == 'within-chrom':
        load_chrom = chrom
    else:
        target_chrs = [c.strip() for c in args.regmode_para.split(',')]
        target_chrs_prefixed = [f'chr{c}' if not c.startswith('chr') else c
                                for c in target_chrs]
        load_chrom = None

    if args.qtl_type == 'eqtl':
        phenotype_df, phenotype_pos_df = load_expression(args.expression, load_chrom)
    else:
        phenotype_df, phenotype_pos_df, group_s = load_splicing_phenotypes(
            args.splicing_bed, args.phenotype_groups, load_chrom
        )

    if args.regmode == 'cross-chrom':
        mask = phenotype_pos_df['chr'].isin(target_chrs_prefixed)
        phenotype_pos_df = phenotype_pos_df[mask]
        phenotype_df = phenotype_df.loc[phenotype_pos_df.index]
        if group_s is not None:
            group_s = group_s[group_s.index.isin(phenotype_df.index)]
        logger.info(f"  Filtered to {len(phenotype_df)} phenotypes on target chromosomes: "
                    f"{','.join(target_chrs_prefixed)}")

    cov_path = args.eqtl_covariates if args.qtl_type == 'eqtl' else args.sqtl_covariates
    covariates_df = load_covariates(cov_path, pc_sex=args.pc_sex,
                                     pc_geno=args.pc_geno, pc_peer=args.pc_peer)

    common_samples = sorted(
        set(genotype_df.columns) & set(phenotype_df.columns) & set(covariates_df.index)
    )
    logger.info(f"Sample overlap: genotype={genotype_df.shape[1]}, "
                f"phenotype={phenotype_df.shape[1]}, covariates={len(covariates_df)}, "
                f"common={len(common_samples)}")

    genotype_df = genotype_df[common_samples]
    phenotype_df = phenotype_df[common_samples]
    covariates_df = covariates_df.loc[common_samples]

    phenotype_pos_df['chr'] = phenotype_pos_df['chr'].astype(str).apply(
        lambda x: x if x.startswith('chr') else f'chr{x}'
    )
    variant_df['chrom'] = variant_df['chrom'].astype(str).apply(
        lambda x: x if x.startswith('chr') else f'chr{x}'
    )

    if args.regmode == 'within-chrom':
        cis_window = int(args.regmode_para)
        cen_half_width = (cen_end - cen_start) // 2

        if cis_window == 0:
            chrom_len = HG38_CHROM_SIZES[chrom]
            effective_window = max(cen_mid, chrom_len - cen_mid)
            logger.info(f"Whole-chromosome mode: chrom_len={chrom_len:,}, "
                        f"effective_window={effective_window:,}")
        else:
            effective_window = cis_window + cen_half_width
            logger.info(f"Centromere width: {(cen_end - cen_start)/1e6:.1f}Mb, "
                        f"effective cis-window: +/-{effective_window/1e6:.1f}Mb from midpoint")

        cis_df = run_cis(
            genotype_df, variant_df,
            phenotype_df, phenotype_pos_df,
            covariates_df, effective_window,
            out_path, chrom, args.seed, args.fdr_method,
            qtl_type=args.qtl_type, group_s=group_s
        )

        if cis_df is None:
            pd.DataFrame(columns=['sv_id', 'gene_id', 'pval_nominal', 'slope', 'slope_se']).to_csv(
                out_path / f"{chrom}_sv_{args.qtl_type}_association.tsv", sep='\t', index=False)
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, 'No testable phenotypes', ha='center', va='center', transform=ax.transAxes)
            fig.savefig(out_path / f"{chrom}_{args.qtl_type}_summary_plot.png", dpi=72)
            plt.close(fig)
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, 'No testable phenotypes', ha='center', va='center', transform=ax.transAxes)
            fig.savefig(out_path / f"{chrom}_{args.qtl_type}_detail_plot.png", dpi=72)
            plt.close(fig)
            logger.info("Done (no testable phenotypes).")
        else:
            pairs_df = build_association_table(out_path, chrom, args.qtl_type)
            assoc_file = out_path / f"{chrom}_sv_{args.qtl_type}_association.tsv"
            pairs_df.to_csv(assoc_file, sep='\t', index=False)
            logger.info(f"  Saved association table: {assoc_file}")

            make_summary_plot(cis_df, pairs_df, sv_summary, chrom, out_path, args.qtl_type)
            make_detail_plot(cis_df, pairs_df, chrom, out_path,
                             phenotype_pos_df, cen_start, cen_end, args.qtl_type)
            try:
                make_qq_diagnostics(cis_df, pairs_df, chrom, out_path, args.qtl_type)
            except Exception as _e:
                logger.warning(f"qq_diagnostics failed: {_e}")

            sig_df = cis_df[cis_df['qval'].notna() & (cis_df['qval'] <= 0.1)].copy()
            if len(sig_df) > 0:
                sig_df = sig_df.sort_values('qval')
                sig_csv = out_path / f"{chrom}_{args.qtl_type}_significant_q0.1.csv"
                sig_df.to_csv(sig_csv, index=True)
                logger.info(f"  Saved {len(sig_df)} significant results (q<0.1): {sig_csv}")
            else:
                logger.info(f"  No results with q < 0.1")

    else:
        trans_df = run_trans(
            genotype_df, phenotype_df, phenotype_pos_df, covariates_df,
            chrom, sv_r2, args.pval_threshold, out_path
        )

        if len(trans_df) > 0:
            make_trans_plots(trans_df, phenotype_pos_df, chrom, out_path,
                             sv_r2, args.pval_threshold, args.qtl_type)
        else:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, 'No significant trans associations',
                    ha='center', va='center', transform=ax.transAxes)
            plot_path = out_path / f"{chrom}_trans_{args.qtl_type}_plot.png"
            fig.savefig(plot_path, dpi=72)
            plt.close(fig)
            logger.info(f"  No trans results — wrote placeholder plot: {plot_path}")

    logger.info("Done.")


if __name__ == "__main__":
    main()
