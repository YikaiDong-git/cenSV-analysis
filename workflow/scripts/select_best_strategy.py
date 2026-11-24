#!/usr/bin/env python3
import argparse
import os
import sys
import random
import shutil
from collections import Counter

import numpy as np
import pandas as pd

ALL_STRATS = ['edge500', 'edge1000', 'edge2000', 'fksize100k', 'fksize200k', 'fksize500k',
              'midpoint500', 'midpoint1000']

BEST6_STRATS = [s for s in ALL_STRATS if 'midpoint' not in s]

SEED = 42


def read_ragged_tsv(fpath):
    with open(fpath) as f:
        header = f.readline().strip().split('\t')
        ncols = len(header)
        rows = []
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= ncols:
                rows.append(fields[:ncols])
            else:
                rows.append(fields + [''] * (ncols - len(fields)))
    return pd.DataFrame(rows, columns=header)


def read_dosage_dict(fpath):
    dosage_rows = {}
    with open(fpath) as f:
        header_line = f.readline()
        for line in f:
            parts = line.split('\t', 1)
            if parts:
                dosage_rows[parts[0].strip()] = line
    return header_line, dosage_rows


def read_summary_dict(fpath):
    summary_rows = {}
    with open(fpath) as f:
        header_line = f.readline()
        for line in f:
            parts = line.split('\t', 1)
            if parts:
                summary_rows[parts[0].strip()] = line
    return header_line, summary_rows


def process_mode(mode_name, strats, en_base_dir, imp_base_dir, output_dir, suffix, chrom):
    sv_data = {}
    for strat in strats:
        metrics_path = os.path.join(en_base_dir, f'{strat}{suffix}', chrom, 'sv_metrics.tsv')
        if not os.path.exists(metrics_path):
            continue
        try:
            df = read_ragged_tsv(metrics_path)
            df['cv_r2'] = pd.to_numeric(df['cv_r2'], errors='coerce')
            df['n_nonzero_coef'] = pd.to_numeric(df['n_nonzero_coef'], errors='coerce')
            df['sv_index'] = pd.to_numeric(df['sv_index'], errors='coerce')
            for _, row in df.iterrows():
                si = row['sv_index']
                if pd.isna(si):
                    continue
                si = int(si)
                if si not in sv_data:
                    sv_data[si] = {}
                cvr2 = row['cv_r2'] if not pd.isna(row['cv_r2']) else 0.0
                nnz = int(row['n_nonzero_coef']) if not pd.isna(row['n_nonzero_coef']) else 999
                sv_data[si][strat] = (cvr2, nnz)
        except Exception as e:
            print(f"  WARN: Failed to read {metrics_path}: {e}", file=sys.stderr)

    if not sv_data:
        return [], {}

    dosage_data = {}
    summary_data = {}
    dosage_header = None
    summary_header = None
    samples_file = None

    for strat in strats:
        imp_dir = os.path.join(imp_base_dir, f'{strat}{suffix}', chrom)
        dosage_path = os.path.join(imp_dir, f'{chrom}_imputed_dosages.tsv')
        summary_path = os.path.join(imp_dir, f'{chrom}_imputation_summary.tsv')
        samples_path = os.path.join(imp_dir, f'{chrom}_imputed_samples.txt')

        if os.path.exists(dosage_path):
            hdr, rows = read_dosage_dict(dosage_path)
            dosage_data[strat] = rows
            if dosage_header is None:
                dosage_header = hdr
        if os.path.exists(summary_path):
            hdr, rows = read_summary_dict(summary_path)
            summary_data[strat] = rows
            if summary_header is None:
                summary_header = hdr
        if samples_file is None and os.path.exists(samples_path):
            samples_file = samples_path

    random.seed(SEED)
    selection = []
    dosage_selection = []

    for si in sorted(sv_data.keys()):
        strat_dict = sv_data[si]
        candidates = [(cvr2, -nnz, random.random(), s) for s, (cvr2, nnz) in strat_dict.items()]
        candidates.sort(reverse=True)

        best_cvr2, neg_nnz, _, best_s = candidates[0]

        dosage_strat = None
        dosage_cvr2 = None
        si_str = str(si)
        for cvr2, neg_n, _, s in candidates:
            if s in dosage_data and si_str in dosage_data[s]:
                dosage_strat = s
                dosage_cvr2 = cvr2
                break

        selection.append((si, best_s, best_cvr2, int(-neg_nnz), len(strat_dict), dosage_strat is not None))
        if dosage_strat is not None:
            dosage_selection.append((si, dosage_strat, dosage_cvr2))

    out_dir = os.path.join(output_dir, f'{mode_name}{suffix}', chrom)
    os.makedirs(out_dir, exist_ok=True)

    if dosage_header and dosage_selection:
        with open(os.path.join(out_dir, f'{chrom}_imputed_dosages.tsv'), 'w') as fout:
            fout.write(dosage_header)
            for si, ds, _ in dosage_selection:
                fout.write(dosage_data[ds][str(si)])
        print(f"    Dosages: {len(dosage_selection)} SVs written")

    if summary_header:
        summary_written = 0
        with open(os.path.join(out_dir, f'{chrom}_imputation_summary.tsv'), 'w') as fout:
            fout.write(summary_header)
            for si, ds, _ in dosage_selection:
                si_str = str(si)
                if ds in summary_data and si_str in summary_data[ds]:
                    fout.write(summary_data[ds][si_str])
                    summary_written += 1
        print(f"    Summary: {summary_written} SVs written")

    if samples_file:
        shutil.copy2(samples_file, os.path.join(out_dir, f'{chrom}_imputed_samples.txt'))

    log_rows = []
    dosage_map = {si: (ds, dc) for si, ds, dc in dosage_selection}
    for si, best_s, best_cvr2, nnz, n_strats, has_dosage in selection:
        row = {
            'sv_index': si,
            'best_strategy_by_cvr2': best_s,
            'best_cv_r2': best_cvr2,
            'n_nonzero_coef': nnz,
            'n_strategies_available': n_strats,
            'has_dosage': has_dosage,
        }
        if si in dosage_map:
            row['dosage_strategy'] = dosage_map[si][0]
            row['dosage_cv_r2'] = dosage_map[si][1]
        else:
            row['dosage_strategy'] = ''
            row['dosage_cv_r2'] = ''
        for s in strats:
            row[f'cv_r2_{s}'] = sv_data[si][s][0] if si in sv_data and s in sv_data[si] else ''
        log_rows.append(row)

    pd.DataFrame(log_rows).to_csv(
        os.path.join(out_dir, f'{chrom}_best_strategy_selection.tsv'), sep='\t', index=False)

    strat_wins = Counter(ds for _, ds, _ in dosage_selection)
    dosage_cvr2_vals = [dc for _, _, dc in dosage_selection]
    n_fallback = sum(1 for si, best_s, _, _, _, _ in selection
                     if si in dosage_map and dosage_map[si][0] != best_s)

    return selection, {
        'total_svs_in_metrics': len(selection),
        'total_svs_with_dosage': len(dosage_selection),
        'n_fallback': n_fallback,
        'mean_dosage_cvr2': np.mean(dosage_cvr2_vals) if dosage_cvr2_vals else 0,
        'strat_wins': strat_wins,
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Select best strategy per SV based on cv_r2.')
    parser.add_argument('--en-base-dir', required=True,
                        help='Base directory for elastic net metrics '
                             '(e.g. results_v9/elastic_net/decompose/cosine_nj)')
    parser.add_argument('--imp-base-dir', required=True,
                        help='Base directory for imputation dosages '
                             '(e.g. results_v9/imputation/decompose/cosine_nj)')
    parser.add_argument('--suffix', required=True,
                        help='Strategy suffix '
                             '(e.g. _mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200)')
    parser.add_argument('--output', required=True,
                        help='Base output directory for merged results')
    parser.add_argument('--chromosomes', default='all',
                        help='Comma-separated chromosome list, or "all" for chr1-22,X (default: all)')
    parser.add_argument('--mode', default='best8', choices=['best8', 'best6', 'both'],
                        help='Which mode to run (default: best8)')
    return parser.parse_args(argv)


def main():
    args = parse_args()

    if args.chromosomes == 'all':
        chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    else:
        chroms = [c.strip() for c in args.chromosomes.split(',')]

    if args.mode == 'both':
        modes = [
            ('best8', ALL_STRATS),
            ('best6', BEST6_STRATS),
        ]
    elif args.mode == 'best6':
        modes = [('best6', BEST6_STRATS)]
    else:
        modes = [('best8', ALL_STRATS)]

    print(f"EN base dir:  {args.en_base_dir}")
    print(f"IMP base dir: {args.imp_base_dir}")
    print(f"Output dir:   {args.output}")
    print(f"Suffix:       {args.suffix}")
    print(f"Processing {len(chroms)} chromosomes\n")

    for mode_name, strats in modes:
        print(f"\n{'='*80}")
        print(f"MODE: {mode_name} ({len(strats)} strategies: {', '.join(strats)})")
        print(f"{'='*80}")

        total_dosage = 0
        total_fallback = 0
        all_wins = Counter()

        for chrom in chroms:
            print(f"\n  {chrom}:")
            selection, stats = process_mode(
                mode_name, strats,
                args.en_base_dir, args.imp_base_dir, args.output, args.suffix,
                chrom)

            if stats:
                total_dosage += stats['total_svs_with_dosage']
                total_fallback += stats['n_fallback']
                all_wins.update(stats['strat_wins'])
                print(f"    {stats['total_svs_with_dosage']} SVs with dosage, "
                      f"{stats['n_fallback']} fallbacks, "
                      f"mean cv_r2={stats['mean_dosage_cvr2']:.4f}")
                top3 = stats['strat_wins'].most_common(3)
                print(f"    Top winners: {', '.join(f'{s}:{c}' for s, c in top3)}")

        total_all_wins = sum(all_wins.values())
        print(f"\n{'='*80}")
        print(f"TOTAL {mode_name}: {total_dosage} SVs with dosage, {total_fallback} fallbacks")
        print(f"Strategy dosage wins:")
        for s in strats:
            cnt = all_wins.get(s, 0)
            pct = 100 * cnt / total_all_wins if total_all_wins > 0 else 0
            print(f"  {s:15s}: {cnt:5d} ({pct:5.1f}%)")
        print(f"{'='*80}")


if __name__ == '__main__':
    main()
