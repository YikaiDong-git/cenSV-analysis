#!/usr/bin/env python3
"""Describe frequency-error patterns within reported SV type and size strata."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

from cohort_frequency_metrics import _write_tsv, frequency_bin


def error_bounds(rows):
    n=len(rows)
    result=dict(rows=n,all_false_lower_count=sum(r['n_tp']==0 and r['n_unresolved']==0 for r in rows),
        all_false_upper_count=sum(r['n_tp']==0 for r in rows),
        any_false_lower_count=sum(r['n_fp']>0 for r in rows),
        any_false_upper_count=sum(r['n_fp']>0 or r['n_unresolved']>0 for r in rows),
        true_carrier_assertions=sum(r['n_tp'] for r in rows),false_carrier_assertions=sum(r['n_fp'] for r in rows),
        unresolved_carrier_assertions=sum(r['n_unresolved'] for r in rows))
    for k in ('all_false_lower','all_false_upper','any_false_lower','any_false_upper'):
        result[k]=result[k+'_count']/n if n else None
    return result


def size_bin(row):
    size=abs(row['svlen'])
    return 'unsized' if size==0 else '<=1kb' if size<=1000 else '1-7kb' if size<=7000 else '>7kb'


def stratify(manifest_path):
    manifest_path=Path(manifest_path);root=manifest_path.parent;m=json.loads(manifest_path.read_text())
    joint=(root/'joint_analysis/analysis.json').exists()
    analysis_root=root/('joint_analysis' if joint else 'analysis')
    if joint:
        rows=[r for r in json.loads((analysis_root/'analysis.json').read_text())['rows'] if r['dedup_retained']]
    else:
        rows=[r for chrom in m['config']['chromosomes'] for r in
              json.loads((analysis_root/chrom/'analysis.json').read_text())['rows'] if r['dedup_retained']]
    output=[];contrasts=[]
    for population in ('all','af_le_0.5'):
        subset=[r for r in rows if population=='all' or r['af']<=.5]
        for dimensions in [('svtype',),('size_bin',),('svtype','size_bin')]:
            grouped=defaultdict(list)
            for row in subset:
                key=tuple(size_bin(row) if k=='size_bin' else row[k] for k in dimensions)
                grouped[key].append(row)
            for key,group in sorted(grouped.items()):
                label=dict(zip(dimensions,key))
                for band in ('<1%','1-<5%','5-<10%','>=10%'):
                    selected=[r for r in group if frequency_bin(r['maf'])==band]
                    output.append(dict(population=population,stratification='+'.join(dimensions),**label,
                                       frequency_bin=band,**error_bounds(selected)))
                common=[r for r in group if r['maf']>=.05]
                for cutoff in (.01,.05):
                    rare=[r for r in group if r['maf']<cutoff];a=error_bounds(rare);b=error_bounds(common)
                    result=dict(population=population,stratification='+'.join(dimensions),**label,
                        rare_cutoff=cutoff,common_cutoff=.05,rare_rows=len(rare),common_rows=len(common))
                    for endpoint in ('all_false','any_false'):
                        al,au=a[endpoint+'_lower'],a[endpoint+'_upper'];bl,bu=b[endpoint+'_lower'],b[endpoint+'_upper']
                        result.update({endpoint+'_ratio_lower':al/bu if al is not None and bu else None,
                                       endpoint+'_ratio_upper':au/bl if au is not None and bl else None,
                                       endpoint+'_difference_lower':al-bu if al is not None and bu is not None else None,
                                       endpoint+'_difference_upper':au-bl if au is not None and bl is not None else None})
                    contrasts.append(result)
    out=analysis_root/'merged'
    _write_tsv(out/'frequency_by_type_size.tsv',output)
    _write_tsv(out/'type_size_enrichment.tsv',contrasts)
    print(f'Wrote {len(output)} stratum estimates and {len(contrasts)} descriptive comparisons')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,default=Path('results/cohort_frequency_v1/manifest.json'))
    stratify(parser.parse_args().manifest)
