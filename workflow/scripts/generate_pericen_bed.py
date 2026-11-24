#!/usr/bin/env python3
"""Generate pericentromeric BED regions with heterochromatin awareness."""
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--centromere-bed', required=True)
    parser.add_argument('--het-bed', required=True)
    parser.add_argument('--fksize', required=True, type=int)
    parser.add_argument('--output', required=True)
    parser.add_argument('--proximity', type=int, default=1_000_000,
                        help="Max distance to consider het block as adjacent (default: 1 Mb)")
    parser.add_argument('--het-chroms', default='all',
                        help="Comma-separated chromosomes that get het-aware extension, "
                             "or 'all' (default). Chromosomes not listed use raw centromere edges.")
    args = parser.parse_args()

    if args.het_chroms.strip().lower() == 'all':
        het_chroms = None
    else:
        het_chroms = {c.strip() for c in args.het_chroms.split(',') if c.strip()}

    cen = pd.read_csv(args.centromere_bed, sep='\t', header=None,
                      names=['chrom', 'start', 'end'])
    het = pd.read_csv(args.het_bed, sep='\t', header=None,
                      names=['chrom', 'start', 'end', 'type'])

    regions = []
    for _, c in cen.iterrows():
        ch = c['chrom']
        cs, ce = c['start'], c['end']

        left_edge = cs
        right_edge = ce
        ch_het = het[het['chrom'] == ch] if (het_chroms is None or ch in het_chroms) else het.iloc[0:0]

        for _, h in ch_het.iterrows():
            hs, he, ht = h['start'], h['end'], h['type']
            if ht == 'short_arm':
                continue
            if hs >= cs and he <= ce:
                continue
            if he <= cs:
                gap = cs - he
                if gap <= args.proximity:
                    left_edge = min(left_edge, hs)
            elif hs >= ce:
                gap = hs - ce
                if gap <= args.proximity:
                    right_edge = max(right_edge, he)
            else:
                left_edge = min(left_edge, hs)
                right_edge = max(right_edge, he)

        left_start = max(0, left_edge - args.fksize)
        left_end = left_edge
        right_start = right_edge
        right_end = right_edge + args.fksize

        regions.append((ch, left_start, left_end))
        regions.append((ch, right_start, right_end))

    with open(args.output, 'w') as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")


if __name__ == '__main__':
    main()
