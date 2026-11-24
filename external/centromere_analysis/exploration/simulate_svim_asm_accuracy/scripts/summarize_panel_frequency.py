#!/usr/bin/env python3
"""Describe MAF filtering of the frozen observed catalog, without labeling errors."""
import argparse
from collections import Counter
import csv
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path


THRESHOLDS = ('0', '0.005', '0.01', '0.05', '0.1')


def maf_from_export(row, haplotypes=392):
    ac = int(row['ac'])
    if not 0 <= ac <= haplotypes:
        raise ValueError('AC outside the specified haplotype denominator')
    if abs(Decimal(row['af']) - Decimal(ac) / haplotypes) > Decimal('0.00000050001'):
        raise ValueError('Exported AF disagrees with AC / haplotype denominator')
    return Fraction(min(ac, haplotypes - ac), haplotypes)


def summarize(rows, haplotypes=392):
    counts = {t: Counter() for t in THRESHOLDS}
    bins = Counter(); types = Counter(); minor_counts = Counter(); n = 0; seen = set()
    for row in rows:
        key = row['catalog_row_id']
        if key in seen:
            raise ValueError('Duplicate catalog row identifier')
        seen.add(key)
        maf = maf_from_export(row, haplotypes)
        n += 1; minor_counts[str(maf * haplotypes)] += 1
        label = ('<1%' if maf < Fraction('0.01') else '1--<5%' if maf < Fraction('0.05')
                 else '5--<10%' if maf < Fraction('0.1') else '>=10%')
        bins[label] += 1
        for threshold in THRESHOLDS:
            if maf < Fraction(threshold):
                continue
            c = counts[threshold]; c['retained'] += 1
            c['all_alpha_rich'] += row.get('all_alleles_alpha_rich') == 'True'
            c['any_TE_rich'] += int(row.get('n_TE_rich_alleles') or 0) > 0
            size = abs(int(float(row['svlen'])))
            if row['svtype'] != 'BND' and size > 0:
                c['sized_records'] += 1
                c['size_le_7kb'] += size <= 7000
            types[threshold, row['svtype']] += 1
    retention = [dict(min_maf=t, retained=counts[t]['retained'], removed=n-counts[t]['retained'],
        fraction_retained=counts[t]['retained']/n if n else None,
        **{k:v for k,v in counts[t].items() if k!='retained'}) for t in THRESHOLDS]
    return dict(n_catalog_rows=n, haplotype_denominator=haplotypes,
        frequency_definition=f'MAF=min(AF,1-AF); exact min(AC,{haplotypes}-AC)/{haplotypes} after agreement with rounded exported AF',
        unit='One retained catalog row with aggregate alternate-presence frequency; not a sequence-resolved alternate allele or an independent mutation',
        baseline_filters='Previously retained AC>=2, AF<=0.99 and within-chromosome coordinate/type/length deduplication; additional MAF filters applied here',
        error_rates_by_maf=None,
        interpretation='Observed frequency distribution and retention only. No truth label exists for these natural catalog rows; the simulations do not assign population MAF.',
        retention=retention, bins=[dict(maf_bin=k, records=v, fraction=v/n if n else None) for k,v in bins.items()],
        minor_allele_counts=dict(minor_counts),
        by_threshold_type=[dict(min_maf=t, svtype=k, records=v) for (t,k),v in sorted(types.items())])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=Path('results/panel_repeat_composition.tsv'))
    p.add_argument('--output-prefix',type=Path,default=Path('results/panel_frequency'))
    p.add_argument('--haplotypes',type=int,default=392)
    args=p.parse_args()
    with args.input.open() as f:
        result=summarize(csv.DictReader(f,delimiter='\t'),args.haplotypes)
    expected=json.loads(Path('results/panel_repeat_composition.json').read_text())['n_catalog_rows']
    assert result['n_catalog_rows']==expected
    result['source']=dict(path=str(args.input.resolve()),sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),bytes=args.input.stat().st_size)
    args.output_prefix.parent.mkdir(parents=True,exist_ok=True)
    args.output_prefix.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    for name in ['retention','bins','by_threshold_type']:
        rows=result[name];fields=list(dict.fromkeys(k for r in rows for k in r))
        with args.output_prefix.with_name(args.output_prefix.name+'_'+name).with_suffix('.tsv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    print(json.dumps(dict(rows=result['n_catalog_rows'],retention=result['retention']),indent=2))


if __name__=='__main__':
    main()
