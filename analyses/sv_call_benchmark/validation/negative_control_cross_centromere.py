#!/usr/bin/env python3
"""Negative control: Cross-centromere test for distal cis-QTL associations."""
import gzip
import logging
import sys
from pathlib import Path

import matplotlib
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                    handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

BD = Path(CENSV_ROOT + "")
SV_PARAM = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
OUTDIR = BD / "claude_t" / "analysis" / "distal_cis"
OUTDIR.mkdir(parents=True, exist_ok=True)

N_RANDOM_GENES = 100
MIN_DIST_FROM_CEN = 10_000_000
RNG = np.random.default_rng(42)



def load_centromeres():
    path = BD / "resources" / "grch38_centromeres.bed"
    df = pd.read_csv(path, sep='\t', header=None, names=['chrom', 'start', 'end'])
    df['midpoint'] = ((df['start'] + df['end']) / 2).astype(int)
    return df.set_index('chrom')


def load_expression_with_positions():
    path = BD / "resources" / "MAGE" / "expression_bed.bed.gz"
    df = pd.read_csv(path, sep='\t', compression='gzip')
    chr_col = '#Chr' if '#Chr' in df.columns else '#chr'
    meta_cols = [chr_col, 'start', 'end', 'ID']
    samples = [c for c in df.columns if c not in meta_cols]

    gene_info = df[[chr_col, 'start', 'end', 'ID']].copy()
    gene_info.columns = ['chrom', 'start', 'end', 'gene_id']

    expr = df[samples].copy()
    expr.index = df['ID'].values

    return expr, samples, gene_info


def load_covariates():
    path = BD / "resources" / "MAGE" / "eqtl_covariates.tab.gz"
    df = pd.read_csv(path, sep='\t', compression='gzip', index_col=0)
    df.loc['sex'] = df.loc['sex'].map({'XY': 1, 'XX': 0}).fillna(0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
    return df.T


def load_sv_dosages(sv_chrom):
    imp_dir = BD / "claude_t" / "results_v9_combined" / "imputation" / \
              "decompose" / "cosine_nj" / f"best6_{SV_PARAM}" / sv_chrom
    dosage_f = imp_dir / f"{sv_chrom}_imputed_dosages.tsv"
    if not dosage_f.exists():
        return None
    return pd.read_csv(dosage_f, sep='\t')


def load_all_genes_final():
    return pd.read_csv(BD / "claude_t" / "results_v9_combined" / "all_genes_final.tsv", sep='\t')



def classify_arm(gene_pos, cen_midpoint):
    return 'p' if gene_pos < cen_midpoint else 'q'


def get_genes_on_arm(gene_info, chrom, cen_midpoint, arm, min_dist):
    mask = gene_info['chrom'] == chrom
    genes = gene_info[mask].copy()
    if arm == 'p':
        genes = genes[genes['start'] < cen_midpoint]
        genes = genes[(cen_midpoint - genes['start']) > min_dist]
    else:
        genes = genes[genes['start'] > cen_midpoint]
        genes = genes[(genes['start'] - cen_midpoint) > min_dist]
    return genes



def test_association(sv_vals, gene_expr, cov_matrix):
    n = len(sv_vals)
    if n < 50:
        return np.nan, np.nan

    X = np.column_stack([sv_vals, cov_matrix, np.ones(n)])
    y = gene_expr

    if np.std(sv_vals) < 1e-10 or np.std(y) < 1e-10:
        return np.nan, np.nan

    try:
        beta, residuals, rank, _ = np.linalg.lstsq(X, y, rcond=None)
        if rank < X.shape[1]:
            return np.nan, np.nan
        resid = y - X @ beta
        dof = n - X.shape[1]
        if dof <= 0:
            return np.nan, np.nan
        mse = np.sum(resid ** 2) / dof
        XtX_inv = np.linalg.inv(X.T @ X)
        se = np.sqrt(mse * XtX_inv[0, 0])
        if se <= 0:
            return np.nan, np.nan
        t_stat = beta[0] / se
        p_val = 2 * scipy_stats.t.sf(abs(t_stat), dof)
        return t_stat, p_val
    except (np.linalg.LinAlgError, ValueError):
        return np.nan, np.nan


def test_sv_vs_genes(sv_row, sample_cols, gene_ids, expr_df, cov_df):
    common = sorted(set(sample_cols) & set(expr_df.columns) & set(cov_df.index))
    if len(common) < 100:
        return []

    sv_vals = sv_row[common].values.flatten().astype(float)
    C = cov_df.loc[common].values

    results = []
    for gid in gene_ids:
        gid_base = gid.split('.')[0]
        matches = [g for g in expr_df.index if g.split('.')[0] == gid_base]
        if not matches:
            continue
        y = expr_df.loc[matches[0], common].values.astype(float)
        t_stat, p_val = test_association(sv_vals, y, C)
        if not np.isnan(t_stat):
            results.append((gid, t_stat, p_val))
    return results



def identify_top_distal_svs(genes_final, n_top=10):
    cis = genes_final[genes_final['cis_trans'] == 'cis'].copy()
    cis['dist_num'] = pd.to_numeric(cis['dist_Mb'], errors='coerce')
    distal = cis[cis['dist_num'] > 10].copy()

    if len(distal) == 0:
        log.warning("No distal cis genes found!")
        return pd.DataFrame()

    distal = distal.sort_values('dist_num', ascending=False)

    seen_svs = set()
    top_svs = []
    for _, row in distal.iterrows():
        key = (row['sv_chrom'], row['sv_id'])
        if key not in seen_svs:
            seen_svs.add(key)
            top_svs.append(row)
            if len(top_svs) >= n_top:
                break

    return pd.DataFrame(top_svs)


def run_cross_centromere_test():
    log.info("=" * 70)
    log.info("CROSS-CENTROMERE NEGATIVE CONTROL")
    log.info("=" * 70)

    log.info("Loading data...")
    centromeres = load_centromeres()
    expr_df, expr_samples, gene_info = load_expression_with_positions()
    cov_df = load_covariates()
    genes_final = load_all_genes_final()

    top_svs = identify_top_distal_svs(genes_final, n_top=10)
    if len(top_svs) == 0:
        log.error("No distal SVs found. Exiting.")
        return
    log.info(f"Top {len(top_svs)} distal SVs selected:")
    for _, sv in top_svs.iterrows():
        log.info(f"  {sv['sv_chrom']} {sv['sv_id']} -> {sv['final_symbol']} "
                 f"({sv['dist_num']:.1f} Mb)")

    all_results = []

    for sv_idx_iter, (_, sv_row_info) in enumerate(top_svs.iterrows()):
        sv_chrom = sv_row_info['sv_chrom']
        sv_id = sv_row_info['sv_id']
        gene_sym = sv_row_info['final_symbol']
        gene_ensg = sv_row_info['ENSG']
        dist_mb = sv_row_info['dist_num']

        log.info(f"\n--- SV {sv_idx_iter+1}/{len(top_svs)}: {sv_chrom} {sv_id} "
                 f"(hit gene: {gene_sym}, {dist_mb:.1f} Mb) ---")

        dosage_df = load_sv_dosages(sv_chrom)
        if dosage_df is None:
            log.warning(f"  No dosage file for {sv_chrom}")
            continue

        try:
            sv_numeric_idx = int(sv_id.split('_')[1])
        except (IndexError, ValueError):
            log.warning(f"  Cannot parse sv_id: {sv_id}")
            continue

        sv_dosage_row = dosage_df[dosage_df['sv_index'] == sv_numeric_idx]
        if len(sv_dosage_row) == 0:
            log.warning(f"  SV index {sv_numeric_idx} not found in dosage matrix")
            continue
        sv_dosage_row = sv_dosage_row.iloc[0]
        sample_cols = [c for c in dosage_df.columns if c != 'sv_index']

        if sv_chrom not in centromeres.index:
            log.warning(f"  No centromere coords for {sv_chrom}")
            continue
        cen_mid = centromeres.loc[sv_chrom, 'midpoint']

        gene_matches = gene_info[gene_info['gene_id'].str.startswith(
            gene_ensg.split('.')[0])]
        if len(gene_matches) == 0:
            log.warning(f"  Gene {gene_ensg} not found in expression BED")
            continue

        hit_gene_pos = gene_matches.iloc[0]['start']
        hit_arm = classify_arm(hit_gene_pos, cen_mid)
        opp_arm = 'q' if hit_arm == 'p' else 'p'
        log.info(f"  Centromere midpoint: {cen_mid:,}  Hit gene arm: {hit_arm}")

        same_arm_genes = get_genes_on_arm(gene_info, sv_chrom, cen_mid,
                                          hit_arm, MIN_DIST_FROM_CEN)
        opp_arm_genes = get_genes_on_arm(gene_info, sv_chrom, cen_mid,
                                         opp_arm, MIN_DIST_FROM_CEN)

        hit_base = gene_ensg.split('.')[0]
        same_arm_genes = same_arm_genes[
            ~same_arm_genes['gene_id'].str.startswith(hit_base)]

        log.info(f"  Same arm ({hit_arm}): {len(same_arm_genes)} genes available")
        log.info(f"  Opposite arm ({opp_arm}): {len(opp_arm_genes)} genes available")

        n_same = min(N_RANDOM_GENES, len(same_arm_genes))
        n_opp = min(N_RANDOM_GENES, len(opp_arm_genes))

        if n_same < 10 or n_opp < 10:
            log.warning(f"  Too few genes on one arm (same={n_same}, opp={n_opp}). "
                        f"Skipping.")
            continue

        same_sample = same_arm_genes.sample(n=n_same, random_state=42 + sv_idx_iter)
        opp_sample = opp_arm_genes.sample(n=n_opp, random_state=42 + sv_idx_iter)

        log.info(f"  Testing {n_same} same-arm genes...")
        same_results = test_sv_vs_genes(sv_dosage_row, sample_cols,
                                        same_sample['gene_id'].tolist(),
                                        expr_df, cov_df)
        log.info(f"  Testing {n_opp} opposite-arm genes...")
        opp_results = test_sv_vs_genes(sv_dosage_row, sample_cols,
                                       opp_sample['gene_id'].tolist(),
                                       expr_df, cov_df)

        same_t = np.array([abs(r[1]) for r in same_results])
        opp_t = np.array([abs(r[1]) for r in opp_results])

        if len(same_t) < 5 or len(opp_t) < 5:
            log.warning(f"  Too few successful tests (same={len(same_t)}, "
                        f"opp={len(opp_t)}). Skipping.")
            continue

        mw_stat, mw_p = scipy_stats.mannwhitneyu(same_t, opp_t,
                                                   alternative='greater')

        log.info(f"  Same-arm |t|: mean={same_t.mean():.3f}, "
                 f"median={np.median(same_t):.3f}, n={len(same_t)}")
        log.info(f"  Opposite-arm |t|: mean={opp_t.mean():.3f}, "
                 f"median={np.median(opp_t):.3f}, n={len(opp_t)}")
        log.info(f"  Mann-Whitney (same > opp): U={mw_stat:.0f}, p={mw_p:.4e}")

        all_results.append({
            'sv_chrom': sv_chrom,
            'sv_id': sv_id,
            'hit_gene': gene_sym,
            'hit_ensg': gene_ensg,
            'dist_Mb': dist_mb,
            'hit_arm': hit_arm,
            'cen_midpoint': cen_mid,
            'n_same_tested': len(same_t),
            'n_opp_tested': len(opp_t),
            'mean_t_same': same_t.mean(),
            'median_t_same': np.median(same_t),
            'mean_t_opp': opp_t.mean(),
            'median_t_opp': np.median(opp_t),
            'mw_U': mw_stat,
            'mw_p': mw_p,
            'same_gt_1_frac': (same_t > 1).mean(),
            'opp_gt_1_frac': (opp_t > 1).mean(),
            'same_t_values': same_t.tolist(),
            'opp_t_values': opp_t.tolist(),
        })

    if not all_results:
        log.error("No SVs could be tested. Exiting.")
        return

    results_df = pd.DataFrame(all_results)

    log.info("\n" + "=" * 70)
    log.info("SUMMARY")
    log.info("=" * 70)

    all_same_t = np.concatenate(results_df['same_t_values'].values)
    all_opp_t = np.concatenate(results_df['opp_t_values'].values)

    mw_overall_stat, mw_overall_p = scipy_stats.mannwhitneyu(
        all_same_t, all_opp_t, alternative='greater')

    log.info(f"Pooled same-arm: mean |t| = {all_same_t.mean():.3f} "
             f"(n={len(all_same_t)})")
    log.info(f"Pooled opposite-arm: mean |t| = {all_opp_t.mean():.3f} "
             f"(n={len(all_opp_t)})")
    log.info(f"Pooled Mann-Whitney (same > opp): U={mw_overall_stat:.0f}, "
             f"p={mw_overall_p:.4e}")

    same_ttest = scipy_stats.ttest_1samp(all_same_t, 1.0, alternative='greater')
    opp_ttest = scipy_stats.ttest_1samp(all_opp_t, 1.0, alternative='greater')
    log.info(f"Same-arm mean |t| > 1: t={same_ttest.statistic:.3f}, "
             f"p={same_ttest.pvalue:.4e}")
    log.info(f"Opposite-arm mean |t| > 1: t={opp_ttest.statistic:.3f}, "
             f"p={opp_ttest.pvalue:.4e}")

    n_same_wins = (results_df['mean_t_same'] > results_df['mean_t_opp']).sum()
    log.info(f"Per-SV same > opp: {n_same_wins}/{len(results_df)}")

    out_tsv = OUTDIR / "cross_centromere.tsv"
    save_df = results_df.drop(columns=['same_t_values', 'opp_t_values'])
    save_df.to_csv(out_tsv, sep='\t', index=False)
    log.info(f"Saved: {out_tsv}")

    make_plots(results_df, all_same_t, all_opp_t, mw_overall_p)


def make_plots(results_df, all_same_t, all_opp_t, mw_overall_p):
    n_svs = len(results_df)
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    ax = axes[0, 0]
    bins = np.linspace(0, max(np.percentile(all_same_t, 99),
                              np.percentile(all_opp_t, 99)), 40)
    ax.hist(all_same_t, bins=bins, alpha=0.6, color='coral', label='Same arm',
            density=True, edgecolor='black', linewidth=0.5)
    ax.hist(all_opp_t, bins=bins, alpha=0.6, color='steelblue',
            label='Opposite arm', density=True, edgecolor='black', linewidth=0.5)
    ax.axvline(1.0, color='black', linestyle='--', alpha=0.4, label='|t|=1')
    ax.axvline(all_same_t.mean(), color='coral', linestyle='-', alpha=0.8,
               linewidth=2)
    ax.axvline(all_opp_t.mean(), color='steelblue', linestyle='-', alpha=0.8,
               linewidth=2)
    ax.set_xlabel('|t-statistic|')
    ax.set_ylabel('Density')
    ax.set_title(f'(A) Pooled |t| distribution\n'
                 f'Mann-Whitney p = {mw_overall_p:.2e}  '
                 f'(same mean={all_same_t.mean():.2f}, '
                 f'opp mean={all_opp_t.mean():.2f})')
    ax.legend(fontsize=9)

    ax = axes[0, 1]
    x_pos = np.arange(n_svs)
    width = 0.35
    ax.bar(x_pos - width/2, results_df['mean_t_same'], width, color='coral',
           label='Same arm', edgecolor='black', linewidth=0.5)
    ax.bar(x_pos + width/2, results_df['mean_t_opp'], width, color='steelblue',
           label='Opposite arm', edgecolor='black', linewidth=0.5)
    ax.axhline(1.0, color='black', linestyle='--', alpha=0.4)
    ax.set_xticks(x_pos)
    labels = [f"{r['hit_gene']}\n({r['sv_chrom']})"
              for _, r in results_df.iterrows()]
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha='right')
    ax.set_ylabel('Mean |t-statistic|')
    ax.set_title('(B) Per-SV mean |t|: same vs opposite arm')
    ax.legend(fontsize=9)

    ax = axes[1, 0]
    pvals = results_df['mw_p'].values
    colors = ['coral' if p < 0.05 else 'gray' for p in pvals]
    ax.bar(x_pos, -np.log10(pvals + 1e-300), color=colors,
           edgecolor='black', linewidth=0.5)
    ax.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5,
               label='p = 0.05')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha='right')
    ax.set_ylabel('-log10(Mann-Whitney p)')
    n_sig = sum(1 for p in pvals if p < 0.05)
    ax.set_title(f'(C) Per-SV same > opp test\n'
                 f'{n_sig}/{n_svs} significant at p<0.05')
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    ax.bar(x_pos - width/2, results_df['same_gt_1_frac'], width, color='coral',
           label='Same arm', edgecolor='black', linewidth=0.5)
    ax.bar(x_pos + width/2, results_df['opp_gt_1_frac'], width,
           color='steelblue', label='Opposite arm', edgecolor='black',
           linewidth=0.5)
    ax.axhline(0.317, color='black', linestyle='--', alpha=0.4,
               label='Expected under null (31.7%)')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha='right')
    ax.set_ylabel('Fraction of genes with |t| > 1')
    ax.set_title('(D) Fraction |t| > 1: same vs opposite arm\n'
                 '(dashed = null expectation 31.7%)')
    ax.legend(fontsize=8)

    plt.suptitle('Cross-centromere negative control:\n'
                 'do centromeric SV associations cross the centromere boundary?',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    out_png = OUTDIR / "cross_centromere.png"
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"Saved: {out_png}")


if __name__ == '__main__':
    run_cross_centromere_test()
