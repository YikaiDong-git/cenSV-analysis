#!/usr/bin/env python3
"""Proper within-ancestry replication of distal cis-QTL associations."""
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
SUPERPOPS = ['AFR', 'EUR', 'SAS', 'EAS', 'AMR']

def get_n_peer(n_samples):
    if n_samples >= 350:
        return 60
    elif n_samples >= 250:
        return 45
    elif n_samples >= 150:
        return 30
    else:
        return 15


def load_population_labels():
    f = BD / "resources" / "imputation_target_samples.report.tsv"
    df = pd.read_csv(f, sep='\t')
    return dict(zip(df['SampleID'], df['Superpopulation']))


def load_expression():
    path = BD / "resources" / "MAGE" / "expression_bed.bed.gz"
    df = pd.read_csv(path, sep='\t', compression='gzip')
    chr_col = '#Chr' if '#Chr' in df.columns else '#chr'
    meta = [chr_col, 'start', 'end', 'ID']
    samples = [c for c in df.columns if c not in meta]
    expr = df[samples].copy()
    expr.index = df['ID'].values
    return expr, samples


def load_sv_dosages(sv_chrom):
    imp_dir = BD / "claude_t" / "results_v9_combined" / "imputation" / \
              "decompose" / "cosine_nj" / f"best6_{SV_PARAM}" / sv_chrom
    f = imp_dir / f"{sv_chrom}_imputed_dosages.tsv"
    if not f.exists():
        return None
    return pd.read_csv(f, sep='\t')


def compute_peer_within_pop(expr_matrix, n_factors=15):
    X = expr_matrix.values.astype(float)
    X = X - X.mean(axis=0)
    std = X.std(axis=0)
    std[std < 1e-10] = 1
    X = X / std

    n_comp = min(n_factors, X.shape[0] - 1, X.shape[1] - 1)
    if n_comp < 1:
        return pd.DataFrame(index=expr_matrix.index)

    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    peer_factors = U[:, :n_comp] * S[:n_comp]

    cols = [f'PEER{i+1}' for i in range(n_comp)]
    return pd.DataFrame(peer_factors, index=expr_matrix.index, columns=cols)


def ols_test(sv_vals, y_vals, cov_matrix=None):
    n = len(sv_vals)
    ones = np.ones(n)
    if cov_matrix is not None and len(cov_matrix.shape) == 2 and cov_matrix.shape[1] > 0:
        X = np.column_stack([sv_vals, cov_matrix, ones])
    else:
        X = np.column_stack([sv_vals, ones])
    try:
        beta = np.linalg.lstsq(X, y_vals, rcond=None)[0]
        resid = y_vals - X @ beta
        df = n - X.shape[1]
        if df <= 0:
            return {'beta': np.nan, 'se': np.nan, 't': np.nan, 'p': np.nan, 'n': n}
        sigma2 = np.sum(resid**2) / df
        XtX_inv = np.linalg.inv(X.T @ X)
        se = np.sqrt(sigma2 * XtX_inv[0, 0])
        t = beta[0] / se if se > 1e-15 else 0
        p = 2 * scipy_stats.t.sf(abs(t), df)
        return {'beta': beta[0], 'se': se, 't': t, 'p': p, 'n': n}
    except np.linalg.LinAlgError:
        return {'beta': np.nan, 'se': np.nan, 't': np.nan, 'p': np.nan, 'n': n}


def cochrans_q(betas, ses):
    valid = ~(np.isnan(betas) | np.isnan(ses) | (ses <= 0))
    b = betas[valid]
    s = ses[valid]
    if len(b) < 2:
        return np.nan, np.nan, np.nan
    w = 1.0 / (s ** 2)
    b_fixed = np.sum(w * b) / np.sum(w)
    Q = np.sum(w * (b - b_fixed) ** 2)
    df = len(b) - 1
    p_Q = scipy_stats.chi2.sf(Q, df)
    I2 = max(0, (Q - df) / Q) if Q > 0 else 0
    return Q, p_Q, I2


def find_gene_in_expr(ensg, expr_df):
    prefix = ensg.split('.')[0]
    matches = [g for g in expr_df.index if g.startswith(prefix)]
    return matches[0] if matches else None


def load_sex_labels():
    path = BD / "resources" / "MAGE" / "eqtl_covariates.tab.gz"
    df = pd.read_csv(path, sep='\t', compression='gzip', index_col=0)
    sex_row = df.loc['sex']
    return sex_row.map({'XY': 1, 'XX': 0, 'male': 1, 'female': 0}).fillna(0).to_dict()


def load_genotype_pcs():
    path = BD / "resources" / "MAGE" / "eqtl_covariates.tab.gz"
    df = pd.read_csv(path, sep='\t', compression='gzip', index_col=0)
    pc_rows = [r for r in df.index if r.startswith('PC')]
    return df.loc[pc_rows].T


def main():
    log.info("=" * 70)
    log.info("PROPER WITHIN-ANCESTRY REPLICATION (v3)")
    log.info("Covariates: sex + 2 within-pop PCs + GTEx-scaled within-pop PEER")
    log.info("  (n>=350: 60 PEER, n>=250: 45, n>=150: 30, n<150: 15)")
    log.info("=" * 70)

    pop_labels = load_population_labels()
    sex_labels = load_sex_labels()
    genotype_pcs = load_genotype_pcs()
    expr_df, expr_samples = load_expression()

    genes_final = pd.read_csv(
        BD / "claude_t" / "results_v9_combined" / "all_genes_final.tsv", sep='\t')
    cis = genes_final[genes_final['cis_trans'] == 'cis'].copy()
    cis['dist_num'] = pd.to_numeric(cis['dist_Mb'], errors='coerce')

    eqtl_genes = cis[cis['QTL'].str.contains('eqtl', na=False)].copy()
    log.info(f"Testing {len(eqtl_genes)} cis-eQTL genes")

    log.info("Computing within-population PEER factors...")
    pop_peers = {}
    pop_samples_map = {}
    for pop in SUPERPOPS:
        pop_samps = [s for s in expr_samples if pop_labels.get(s) == pop]
        pop_samples_map[pop] = pop_samps
        if len(pop_samps) >= 30:
            n_peer = get_n_peer(len(pop_samps))
            pop_expr = expr_df[pop_samps].T
            peer_df = compute_peer_within_pop(pop_expr, n_factors=n_peer)
            pop_peers[pop] = peer_df
            log.info(f"  {pop}: n={len(pop_samps)}, {peer_df.shape[1]} PEER factors (GTEx: {n_peer})")
        else:
            log.info(f"  {pop}: n={len(pop_samps)} (too few)")

    dosage_cache = {}

    all_results = []

    for gene_i, (_, gene) in enumerate(eqtl_genes.iterrows()):
        sv_chrom = gene['sv_chrom']
        sv_id = gene['sv_id']
        gene_sym = gene['final_symbol']
        gene_ensg = gene['ENSG']
        dist_Mb = gene['dist_num']

        try:
            sv_idx = int(sv_id.split('_')[1])
        except:
            continue

        if sv_chrom not in dosage_cache:
            dosage_cache[sv_chrom] = load_sv_dosages(sv_chrom)
        dosage_df = dosage_cache[sv_chrom]
        if dosage_df is None:
            continue

        sv_row = dosage_df[dosage_df['sv_index'] == sv_idx]
        if len(sv_row) == 0:
            continue

        expr_id = find_gene_in_expr(gene_ensg, expr_df)
        if expr_id is None:
            continue

        sample_cols = [c for c in dosage_df.columns if c != 'sv_index']

        pop_results = {}
        pop_betas = []
        pop_ses = []

        for pop in SUPERPOPS:
            if pop not in pop_peers:
                continue

            pop_samps = [s for s in pop_samples_map[pop]
                         if s in sample_cols and s in expr_df.columns]
            if len(pop_samps) < 30:
                continue

            sv_vals = sv_row[pop_samps].values.flatten().astype(float)
            y_vals = expr_df.loc[expr_id, pop_samps].values.astype(float)

            if np.std(sv_vals) < 1e-10 or np.std(y_vals) < 1e-10:
                pop_results[pop] = {'n': len(pop_samps), 't': 0, 'p': 1.0,
                                     'beta': 0, 'se': np.nan}
                continue

            peer_vals = pop_peers[pop]
            common_pop = sorted(set(pop_samps) & set(peer_vals.index) &
                                set(genotype_pcs.index))
            if len(common_pop) < 30:
                continue

            sv_p = sv_row[common_pop].values.flatten().astype(float)
            y_p = expr_df.loc[expr_id, common_pop].values.astype(float)

            sex_vec = np.array([sex_labels.get(s, 0) for s in common_pop]).astype(float)
            pop_pc_vals = genotype_pcs.loc[common_pop, ['PC1', 'PC2']].values.astype(float)
            peer_p = peer_vals.loc[common_pop].values.astype(float)
            cov_p = np.column_stack([sex_vec, pop_pc_vals, peer_p]).astype(float)

            if np.std(sv_p) < 1e-10 or np.std(y_p) < 1e-10:
                continue

            res = ols_test(sv_p, y_p, cov_p)
            pop_results[pop] = res
            if not np.isnan(res['beta']) and not np.isnan(res['se']) and res['se'] > 0:
                pop_betas.append(res['beta'])
                pop_ses.append(res['se'])

        Q, p_Q, I2 = cochrans_q(np.array(pop_betas), np.array(pop_ses))

        if len(pop_betas) >= 2:
            w = 1.0 / np.array(pop_ses) ** 2
            beta_meta = np.sum(w * np.array(pop_betas)) / np.sum(w)
            se_meta = 1.0 / np.sqrt(np.sum(w))
            z_meta = beta_meta / se_meta
            p_meta = 2 * scipy_stats.norm.sf(abs(z_meta))
        else:
            beta_meta = se_meta = z_meta = p_meta = np.nan

        n_pop_sig = sum(1 for p in SUPERPOPS
                        if pop_results.get(p, {}).get('p', 1) < 0.05)
        signs = [np.sign(pop_results[p]['beta'])
                 for p in SUPERPOPS
                 if p in pop_results and not np.isnan(pop_results[p].get('beta', np.nan))
                 and pop_results[p].get('n', 0) >= 30]
        consistent_sign = len(set(signs)) <= 1 if signs else False

        row = {
            'gene': gene_sym, 'ensg': gene_ensg, 'sv_chrom': sv_chrom,
            'sv_id': sv_id, 'dist_Mb': dist_Mb,
            'risk': gene['final_risk'], 'QTL': gene['QTL'],
            'all_p': gene['qval'],
            'n_pops_tested': sum(1 for p in SUPERPOPS if p in pop_results and pop_results[p].get('n',0)>=30),
            'n_pops_sig': n_pop_sig,
            'consistent_sign': consistent_sign,
            'meta_beta': beta_meta, 'meta_se': se_meta,
            'meta_z': z_meta, 'meta_p': p_meta,
            'Q': Q, 'p_Q': p_Q, 'I2': I2,
        }
        for pop in SUPERPOPS:
            pr = pop_results.get(pop, {})
            row[f'{pop}_n'] = pr.get('n', 0)
            row[f'{pop}_beta'] = pr.get('beta', np.nan)
            row[f'{pop}_se'] = pr.get('se', np.nan)
            row[f'{pop}_t'] = pr.get('t', np.nan)
            row[f'{pop}_p'] = pr.get('p', np.nan)

        all_results.append(row)

        if (gene_i + 1) % 20 == 0:
            log.info(f"  Processed {gene_i+1}/{len(eqtl_genes)} genes...")

    results_df = pd.DataFrame(all_results)
    log.info(f"\n{'='*70}")
    log.info(f"RESULTS: {len(results_df)} genes tested")
    log.info(f"{'='*70}")

    bins = [(0, 1, '<1 Mb'), (1, 5, '1-5 Mb'), (5, 10, '5-10 Mb'), (10, 500, '>10 Mb')]
    log.info("\nWithin-ancestry survival (WITH proper PEER correction):")
    for lo, hi, label in bins:
        sub = results_df[(results_df['dist_Mb'] >= lo) & (results_df['dist_Mb'] < hi)]
        n_surv = (sub['n_pops_sig'] >= 1).sum()
        n_tot = len(sub)
        rate = 100 * n_surv / max(n_tot, 1)
        n_meta = (sub['meta_p'] < 0.05).sum()
        log.info(f"  {label:10s}: {n_surv}/{n_tot} = {rate:.1f}% (pop sig)  "
                 f"{n_meta}/{n_tot} = {100*n_meta/max(n_tot,1):.1f}% (meta sig)")

    prox = results_df[results_df['dist_Mb'] < 1]
    dist = results_df[results_df['dist_Mb'] >= 10]
    a = (dist['n_pops_sig'] >= 1).sum()
    b = len(dist) - a
    c = (prox['n_pops_sig'] >= 1).sum()
    d = len(prox) - c
    if a + b > 0 and c + d > 0:
        odds, pf = scipy_stats.fisher_exact([[a, b], [c, d]])
        log.info(f"\n  Fisher (distal vs proximal): OR={odds:.2f}, p={pf:.4f}")

    log.info("\nKey distal genes (with within-pop PEER):")
    key_genes = ['CWC22', 'PDIK1L', 'CCDC74B', 'TEX29', 'ZNF384', 'SRRM5', 'CYBB']
    for gene in key_genes:
        row = results_df[results_df['gene'] == gene]
        if len(row) == 0:
            continue
        r = row.iloc[0]
        pop_str = []
        for pop in SUPERPOPS:
            p_val = r.get(f'{pop}_p', np.nan)
            n_val = r.get(f'{pop}_n', 0)
            if n_val >= 30 and not np.isnan(p_val):
                sig = '*' if p_val < 0.05 else '~' if p_val < 0.10 else ''
                pop_str.append(f'{pop}:{p_val:.2e}{sig}')
        log.info(f"  {gene:12s} {r['dist_Mb']:5.0f}Mb  meta_p={r['meta_p']:.2e}  "
                 f"I2={r['I2']:.2f}  sign={'Y' if r['consistent_sign'] else 'N'}  "
                 f"sig={r['n_pops_sig']}/{r['n_pops_tested']}  "
                 f"{' '.join(pop_str)}")

    results_df.to_csv(OUTDIR / "ancestry_replication_proper.tsv", sep='\t', index=False)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    ax = axes[0, 0]
    labels, rates_pop, rates_meta, ns = [], [], [], []
    for lo, hi, label in bins:
        sub = results_df[(results_df['dist_Mb'] >= lo) & (results_df['dist_Mb'] < hi)]
        labels.append(f'{label}\n(n={len(sub)})')
        rates_pop.append(100*(sub['n_pops_sig']>=1).sum()/max(len(sub),1))
        rates_meta.append(100*(sub['meta_p']<0.05).sum()/max(len(sub),1))
    x = np.arange(len(labels))
    ax.bar(x - 0.2, rates_pop, 0.35, label='Any pop p<0.05', color='steelblue')
    ax.bar(x + 0.2, rates_meta, 0.35, label='Meta p<0.05', color='coral')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('% surviving')
    ax.set_title('(A) Within-ancestry survival\n(with proper PEER)')
    ax.legend()

    ax = axes[0, 1]
    cwc = results_df[results_df['gene'] == 'CWC22']
    if len(cwc) > 0:
        r = cwc.iloc[0]
        y_pos = []
        b_vals = []
        se_vals = []
        pop_labels_plot = []
        for i, pop in enumerate(SUPERPOPS):
            b = r.get(f'{pop}_beta', np.nan)
            se = r.get(f'{pop}_se', np.nan)
            if not np.isnan(b) and not np.isnan(se) and se > 0:
                y_pos.append(i)
                b_vals.append(b)
                se_vals.append(se)
                pop_labels_plot.append(f"{pop} (n={int(r.get(f'{pop}_n', 0))})")
        if b_vals:
            ax.errorbar(b_vals, y_pos, xerr=[1.96*s for s in se_vals],
                        fmt='o', color='steelblue', capsize=3)
            ax.axvline(0, color='red', linestyle='--', alpha=0.5)
            if not np.isnan(r['meta_beta']):
                ax.axvline(r['meta_beta'], color='green', linestyle='-', alpha=0.5,
                           label=f"Meta: {r['meta_beta']:.3f}")
            ax.set_yticks(y_pos)
            ax.set_yticklabels(pop_labels_plot)
            ax.set_xlabel('Effect size (beta)')
            ax.set_title(f"(B) CWC22 forest plot\nmeta p={r['meta_p']:.2e}, I2={r['I2']:.2f}")
            ax.legend(fontsize=9)

    ax = axes[1, 0]
    ccdc = results_df[results_df['gene'] == 'CCDC74B']
    if len(ccdc) > 0:
        r = ccdc.iloc[0]
        y_pos, b_vals, se_vals, pop_labels_plot = [], [], [], []
        for i, pop in enumerate(SUPERPOPS):
            b = r.get(f'{pop}_beta', np.nan)
            se = r.get(f'{pop}_se', np.nan)
            if not np.isnan(b) and not np.isnan(se) and se > 0:
                y_pos.append(i)
                b_vals.append(b)
                se_vals.append(se)
                pop_labels_plot.append(f"{pop} (n={int(r.get(f'{pop}_n', 0))})")
        if b_vals:
            ax.errorbar(b_vals, y_pos, xerr=[1.96*s for s in se_vals],
                        fmt='o', color='coral', capsize=3)
            ax.axvline(0, color='red', linestyle='--', alpha=0.5)
            if not np.isnan(r['meta_beta']):
                ax.axvline(r['meta_beta'], color='green', linestyle='-', alpha=0.5,
                           label=f"Meta: {r['meta_beta']:.3f}")
            ax.set_yticks(y_pos)
            ax.set_yticklabels(pop_labels_plot)
            ax.set_xlabel('Effect size (beta)')
            ax.set_title(f"(C) CCDC74B forest plot\nmeta p={r['meta_p']:.2e}, I2={r['I2']:.2f}")
            ax.legend(fontsize=9)

    ax = axes[1, 1]
    i2_vals = results_df['I2'].dropna()
    ax.hist(i2_vals, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
    ax.axvline(0.25, color='red', linestyle='--', alpha=0.5, label='Low heterogeneity')
    ax.axvline(0.75, color='darkred', linestyle='--', alpha=0.5, label='High heterogeneity')
    ax.set_xlabel('I-squared')
    ax.set_ylabel('Count')
    ax.set_title('(D) Effect heterogeneity\nacross ancestries')
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(OUTDIR / "ancestry_replication_proper.png", dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"\nSaved ancestry_replication_proper.png")
    log.info("Done!")


if __name__ == '__main__':
    main()
