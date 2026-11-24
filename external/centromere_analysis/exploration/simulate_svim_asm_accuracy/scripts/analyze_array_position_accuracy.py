#!/usr/bin/env python3
"""Postprocess frozen simulations by active-HOR array position and CDR overlap."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
ABS_EDGES = (0, 50_000, 100_000, 250_000, 500_000, 1_000_000, float('inf'))
ABS_LABELS = ('0-50 kb', '50-100 kb', '100-250 kb', '250-500 kb', '500-1000 kb', '>=1000 kb')
REL_EDGES = (0, .2, .4, .6, .8, 1.00000001)
REL_LABELS = ('0-20%', '20-40%', '40-60%', '60-80%', '80-100%')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows):
    if not rows:
        Path(path).write_text('')
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w') as f:
        w = csv.DictWriter(f, fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)


def read_bed(path):
    result = defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith(('#', 'track', 'browser')):
            continue
        col = line.split()
        if len(col) < 3 or not col[0].startswith('chr'):
            raise ValueError(f'Not a BED interval: {line[:100]}')
        a, b = int(col[1]), int(col[2])
        if not 0 <= a < b:
            raise ValueError('Invalid BED bounds')
        result[col[0]].append((a,b))
    merged = {}
    for chrom, intervals in result.items():
        out = []
        for a,b in sorted(intervals):
            if out and a <= out[-1][1]:
                out[-1] = (out[-1][0], max(b,out[-1][1]))
            else:
                out.append((a,b))
        merged[chrom] = out
    return merged


def bin_value(value, edges, labels):
    return next(label for a,b,label in zip(edges,edges[1:],labels) if a <= value < b)


def position(chrom, start, end, arrays, cdrs, convention='midpoint'):
    start, end = int(start), int(end)
    if end < start:
        raise ValueError('Reversed reference edit')
    point = (start + end) / 2 if convention == 'midpoint' else float(start)
    candidates = [(a,b) for a,b in arrays.get(chrom,[]) if a <= point < b]
    if not candidates and start == end:
        candidates = [(a,b) for a,b in arrays.get(chrom,[]) if point == b]
    if not candidates:
        return dict(in_array=False, absolute_bin='Outside active arrays', relative_bin='Outside active arrays',
                    cdr_point='Unannotated' if chrom not in cdrs else 'Outside CDR',
                    cdr_overlap='Unannotated' if chrom not in cdrs else 'Outside CDR')
    if len(candidates) != 1:
        raise ValueError('Overlapping merged arrays')
    a,b = candidates[0]
    distance = min(point-a, b-point)
    if convention == 'nearest_breakpoint':
        distance = min(max(0,min(start-a,b-start)), max(0,min(end-a,b-end)))
    depth = min(1.,2*distance/(b-a))
    overlaps = any(x <= start < y if start == end else max(start,x) < min(end,y)
                   for x,y in cdrs.get(chrom,[]))
    inside_cdr = any(x <= point < y for x,y in cdrs.get(chrom,[]))
    return dict(in_array=True,array_start=a,array_end=b,array_length=b-a,point=point,
                distance_to_end_bp=distance,relative_depth=depth,
                array_boundary_crossing=start<a or end>b,
                absolute_bin=bin_value(distance,ABS_EDGES,ABS_LABELS),
                relative_bin=bin_value(depth,REL_EDGES,REL_LABELS),
                cdr_point=('Inside CDR' if inside_cdr else 'Outside CDR') if chrom in cdrs else 'Unannotated',
                cdr_overlap=('Overlaps CDR' if overlaps else 'No CDR overlap') if chrom in cdrs else 'Unannotated')


def size_bin(size):
    n=abs(int(size))
    return 'Unsized' if n==0 else '<=1 kb' if n<=1000 else '1-7 kb' if n<=7000 else '>7 kb'


def annotate(row, arrays, cdrs):
    result = dict(row, **position(row['chrom'],row['start'],row['end'],arrays,cdrs))
    for mode in ('start','nearest_breakpoint'):
        p=position(row['chrom'],row['start'],row['end'],arrays,cdrs,mode)
        result[mode+'_absolute_bin']=p['absolute_bin']
        result[mode+'_relative_bin']=p['relative_bin']
    return result


def inputs(root,arrays,cdrs,pilot_output=None):
    calls,events=[],[]
    base=root/'results/multichrom_diagnostics_v2/summary'
    for row in read_tsv(base/'calls.tsv'):
        if row['scenario'] not in ('alpha_mixed','control'):
            continue
        start=int(row['start']);end=int(row['end'])
        calls.append(annotate(dict(dataset='23-chromosome',case_id=row['case_id'],chrom=row['call_chrom'],
            source_chrom=row['chrom'],rate=row['array_rate'],svtype=row['called_type'],
            size_bin=size_bin(end-start if row['called_type']=='DEL' else len(json.loads(row['edit'])[2]) if row['edit'] and json.loads(row['edit']) else 0),
            start=start,end=end,status=row['status'],record_id=row['call_id']),arrays,cdrs))
    for row in read_tsv(base/'events.tsv'):
        if row['scenario']!='alpha_mixed':
            continue
        events.append(annotate(dict(dataset='23-chromosome',case_id=row['case_id'],chrom=row['chrom'],rate=row['array_rate'],
            svtype=row['true_type'],size_bin=size_bin(row['size_bp']),start=int(row['start']),end=int(row['end']),
            status=row['status'],recovered=row['recovered']=='True',site_id=row['chrom']+'__'+row['event_id']),arrays,cdrs))
    assert len(calls)==2055 and Counter(r['status'] for r in calls)==dict(tp=600,fp=1354,unresolved=101)
    assert len(events)==1017 and sum(r['recovered'] for r in events)==598
    old=root/'results/complete_truth_v1'
    verification={r['case']:r for r in read_tsv(root/'results/chr1_call_error_v1/strict_metric_verification.tsv')}
    old_counts=Counter()
    for row in read_tsv(old/'calls.tsv'):
        if row['scope']!='production_active':
            continue
        if 'fragmented' in row['scenario']:
            continue
        old_counts[row['case']] += 1
        calls.append(annotate(dict(dataset='Chr1 intact',case_id=row['case'],chrom=row['chrom'],source_chrom='chr1',
            rate=row['array_rate'],svtype=row['kind'],size_bin=size_bin(row['called_edit_size_bp']),
            start=int(row['start']),end=int(row['end']),status=row['status'],record_id=row['id']),arrays,cdrs))
    assert all(count==int(verification[case]['n_calls']) for case,count in old_counts.items())
    for row in read_tsv(old/'events.tsv'):
        if row['scope']!='production_active' or 'fragmented' in row['scenario']:
            continue
        events.append(annotate(dict(dataset='Chr1 intact',case_id=row['case'],chrom=row['chrom'],rate=row['array_rate'],
            svtype=row['kind'],size_bin=size_bin(row['size_bp']),start=int(row['start']),end=int(row['end']),
            status=row['status'],recovered=row['recovered']=='True',site_id=row['chrom']+'__'+row['id']+'__'+row['seed']),arrays,cdrs))
    import pysam
    from parse_complete_truth_calls import parse_complete_calls
    from multichrom_truth_scoring import score_complete_truth
    pilot=root/'results/chm13_chr1_pilot';manifest=load(pilot/'manifest.json')
    with pysam.FastaFile(manifest['reference']) as reference:
        assert hashlib.sha256(reference.fetch('chr1').upper().encode()).hexdigest()==manifest['reference_chr_sha256']
        for e in manifest['events']:
            if e['context'] not in ('array_edge','array_interior'):
                continue
            start=e['position'];end=start+e['size'] if e['kind']=='DEL' else start
            alt='' if e['kind']=='DEL' else e['allele']
            if e['kind']=='DEL':
                assert reference.fetch('chr1',start,end).upper()==e['allele']
            assert hashlib.sha256(e['allele'].encode()).hexdigest()==e['allele_sha256']
            event=dict(id=e['id'],chrom='chr1',start=start,end=end,alt=alt,kind='DEL' if e['kind']=='DEL' else 'INS')
            source=pilot/'cases'/e['id']/'active.vcf.gz'
            parsed=parse_complete_calls(source,reference)
            scored=score_complete_truth(reference,[event],parsed,[])
            if pilot_output is not None:
                pilot_output.mkdir(parents=True,exist_ok=True)
                scored['metadata']=dict(source_vcf=str(source),source_vcf_sha256=sha(source),pilot_manifest_sha256=sha(pilot/'manifest.json'),reference_chr1_sha256=manifest['reference_chr_sha256'],strict_scorer_sha256=sha(ROOT/'scripts/multichrom_truth_scoring.py'))
                (pilot_output/(e['id']+'.json')).write_text(json.dumps(scored,indent=2)+'\n')
            assert scored['metrics']['unresolved_calls']==0 and scored['metrics']['unresolved_events']==0
            for c,z in zip(parsed,scored['calls']):
                assert c['id']==z['id']
                a,b=c['edit'][:2] if c['edit'] else (c['pos'],c['pos'])
                calls.append(annotate(dict(dataset='CHM13 single-site pilot',case_id=e['id'],chrom=c['chrom'],source_chrom='chr1',rate='0.0',
                    svtype=c['kind'],size_bin=size_bin(c['length']),start=a,end=b,status=z['status'],record_id=c['id']),arrays,cdrs))
            ev=scored['events'][0]
            events.append(annotate(dict(dataset='CHM13 single-site pilot',case_id=e['id'],chrom='chr1',rate='0.0',
                svtype=event['kind'],size_bin=size_bin(e['size']),start=start,end=end,status=ev['status'],recovered=ev['recovered'],
                site_id=e['id'],pilot_context=e['context']),arrays,cdrs))
    a=load(root/'results/cohort_frequency_v1/joint_analysis/analysis.json')
    members=defaultdict(list)
    for r in a['memberships']:members[r['haplotype_index']].append(r)
    truths={}
    for h in range(256):
        d=root/f'results/cohort_frequency_v1/joint_cohort/haps/H{h:03d}'
        parsed={r['id']:r for r in load(d/'scoring/parsed_calls_production_active.json')}
        truths[h]={r['id']:r for r in load(d/'truth.json')['truth']}
        for member in members[h]:
            c=parsed[member['source_id'].split('__',1)[1]]
            edit=c.get('edit')
            start,end=(edit[0],edit[1]) if edit else (c['pos'],c['pos'])
            for cut in (0,5):
                status=member[f'status_maf_{cut}']
                if status is None:
                    continue
                calls.append(annotate(dict(dataset='Cohort production' if cut==0 else 'Cohort MAF>=5%',case_id=f'H{h:03d}',
                    chrom=c['chrom'],source_chrom=c['source_query_chrom'],rate='founder+private',svtype=c['kind'],
                    size_bin=size_bin(c['length']),start=start,end=end,status=status,record_id=c['id']),arrays,cdrs))
    for row in a['event_carriers']:
        h=row['haplotype_index'];e=truths[h][row['event_id']]
        for cut in (0,5):
            events.append(annotate(dict(dataset='Cohort production' if cut==0 else 'Cohort MAF>=5%',case_id=f'H{h:03d}',
                chrom=e['chrom'],rate='founder+private',svtype=e['kind'],size_bin=size_bin(max(e['end']-e['start'],len(e['alt']))),
                start=e['start'],end=e['end'],status=row[f'status_maf_{cut}'],recovered=row[f'recovered_maf_{cut}'],
                site_id=e['id'],true_af=row['true_af']),arrays,cdrs))
    for label,n,tp in [('Cohort production',7659,2277),('Cohort MAF>=5%',3905,2081)]:
        selected=[r for r in calls if r['dataset']==label]
        assert len(selected)==n and sum(r['status']=='tp' for r in selected)==tp
    return calls,events


def aggregate(rows, endpoint):
    result=[]
    for field in ('absolute_bin','relative_bin','cdr_point','cdr_overlap',
                  'start_absolute_bin','start_relative_bin','nearest_breakpoint_absolute_bin','nearest_breakpoint_relative_bin'):
        for extra in ((),('svtype',),('size_bin',),('chrom',),('rate',),('svtype','size_bin')):
            groups=defaultdict(list)
            for r in rows:groups[(r['dataset'],r[field],*(r[k] for k in extra))].append(r)
            for keys,group in sorted(groups.items()):
                n=len(group)
                correct=sum(r['status']=='tp' if endpoint=='precision' else bool(r['recovered']) for r in group)
                unresolved=sum(r['status']=='unresolved' for r in group)
                result.append(dict(endpoint=endpoint,dataset=keys[0],position_measure=field,position_bin=keys[1],
                    stratification='+'.join(extra) or 'All',**dict(zip(extra,keys[2:])),n=n,correct=correct,
                    false_or_missed=n-correct-unresolved,unresolved=unresolved,lower=correct/n,upper=(correct+unresolved)/n,
                    unique_cases=len({r['case_id'] for r in group}),unique_sites=len({r['site_id'] for r in group}) if endpoint=='sensitivity' else None))
    return result


def exposure_bp(arrays, chromosomes, measure, lo, hi):
    intervals=[(a,b) for c,v in arrays.items() if c in chromosomes for a,b in v]
    if measure=='relative_bin':
        return sum((b-a)*max(0,min(1,hi)-max(0,lo)) for a,b in intervals)
    return sum(2*max(0,min((b-a)/2,hi)-lo) for a,b in intervals)


def observed(root,arrays,cdrs):
    source=root/'results/panel_repeat_composition.tsv'
    rows=[]
    for r in read_tsv(source):
        start=int(r['pos'])-1
        rows.append(dict(chrom=r['chrom'],svtype=r['svtype'],size_bin=size_bin(r['svlen']),
            **position(r['chrom'],start,start,arrays,cdrs,'start')))
    assert len(rows)==188591
    result=[]
    observed_chromosomes={r['chrom'] for r in rows}
    if not observed_chromosomes <= set(arrays):
        raise ValueError('Observed chromosome lacks array annotation')
    for measure in ('absolute_bin','relative_bin'):
        labels=ABS_LABELS if measure=='absolute_bin' else REL_LABELS
        edges=ABS_EDGES if measure=='absolute_bin' else REL_EDGES
        for chrom in ('ALL',*sorted(observed_chromosomes)):
            selected=[r for r in rows if chrom=='ALL' or r['chrom']==chrom]
            for i,label in enumerate(labels):
                lo,hi=edges[i:i+2]
                exposed=exposure_bp(arrays,observed_chromosomes if chrom=='ALL' else {chrom},measure,lo,hi)
                count=sum(r[measure]==label for r in selected)
                result.append(dict(chrom=chrom,position_measure=measure,position_bin=label,n_rows=count,
                    reference_exposure_bp=exposed,rows_per_Mb=count/exposed*1e6 if exposed else None))
    return rows,result


def figures(output,precision,sensitivity,density):
    labels=['23-chromosome','Chr1 intact','Cohort production','Cohort MAF>=5%']
    colors=['#386f91','#805b9c','#b75c40','#6b8540']
    for measure,bins,name in [('relative_bin',REL_LABELS,'relative_position'),('absolute_bin',ABS_LABELS,'distance_to_end')]:
        fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        for table,ax,title in [(precision,axes[0],'Precision at the reported edit'),(sensitivity,axes[1],'Sensitivity at the true edit')]:
            for label,color in zip(labels,colors):
                indexed={r['position_bin']:r for r in table if r['dataset']==label and r['position_measure']==measure and r['stratification']=='All'}
                y=[100*indexed[b]['lower'] if b in indexed else np.nan for b in bins]
                u=[100*indexed[b]['upper'] if b in indexed else np.nan for b in bins]
                ax.plot(range(len(bins)),y,'o-',label=label,color=color)
                ax.fill_between(range(len(bins)),y,u,color=color,alpha=.12)
            ax.set_title(title,loc='left',fontweight='bold');ax.set_ylim(0,105);ax.set_xticks(range(len(bins)),bins,rotation=25,ha='right')
            ax.set_ylabel('Correct / total (%)');ax.grid(axis='y',alpha=.15);ax.spines[['top','right']].set_visible(False)
        axes[0].legend(frameon=False,fontsize=8)
        fig.supxlabel('Depth from either end: 0% = boundary, 100% = midpoint' if measure=='relative_bin' else 'Distance from nearest active-HOR array end')
        for suffix in ('pdf','png'):fig.savefig(output/f'{name}.{suffix}',dpi=200,bbox_inches='tight')
        plt.close(fig)
    selected={r['position_bin']:r for r in density if r['chrom']=='ALL' and r['position_measure']=='relative_bin'}
    fig,ax=plt.subplots(figsize=(7,3.5),layout='constrained')
    ax.bar(REL_LABELS,[selected[b]['rows_per_Mb'] for b in REL_LABELS],color='#386f91')
    ax.set_ylabel('Observed catalog rows / reference Mb');ax.set_xlabel('Depth from nearest end (100% = array midpoint)')
    ax.set_title('Observed HPRC call density; correctness is not known',loc='left',fontsize=12,fontweight='bold')
    for suffix in ('pdf','png'):fig.savefig(output/f'observed_density.{suffix}',dpi=200,bbox_inches='tight')
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--array-bed',type=Path,default=Path(CENSV_ROOT + '/resources/chm13v2.0_censat_v2.1.active.bed'))
    parser.add_argument('--cdr-bed',type=Path)
    parser.add_argument('--output',type=Path,default=Path('results/array_position_v1'))
    args=parser.parse_args();csv.field_size_limit(10_000_000)
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    arrays=read_bed(args.array_bed);cdrs=read_bed(args.cdr_bed) if args.cdr_bed else {}
    calls,events=inputs(args.root,arrays,cdrs,out/'pilot_strict_recheck')
    p=aggregate(calls,'precision');s=aggregate(events,'sensitivity')
    obs,d=observed(args.root,arrays,cdrs)
    for name,rows in [('calls_annotated.tsv',calls),('events_annotated.tsv',events),('precision.tsv',p),('sensitivity.tsv',s),('observed_density.tsv',d)]:write_tsv(out/name,rows)
    figures(out,p,s,d)
    info=dict(array_bed=str(args.array_bed),array_bed_sha256=sha(args.array_bed),cdr_bed=str(args.cdr_bed) if args.cdr_bed else None,
        cdr_bed_sha256=sha(args.cdr_bed) if args.cdr_bed else None,script_sha256=sha(__file__),
        coordinate_convention='zero-based half-open; midpoint of normalized reference edit; insertion point for INS',
        relative_depth='2 * min(distance to left end, distance to right end) / array length',
        arrays=arrays,observed_chromosomes=sorted({r['chrom'] for r in obs}),n_calls=len(calls),n_truth_occurrences=len(events),n_observed_rows=len(obs),
        observed_outside_arrays=sum(not r['in_array'] for r in obs),call_counts=dict(Counter(r['dataset'] for r in calls)),
        truth_counts=dict(Counter(r['dataset'] for r in events)),source_sha256={str(p):sha(p) for p in [args.root/'results/multichrom_diagnostics_v2/summary/calls.tsv',args.root/'results/multichrom_diagnostics_v2/summary/events.tsv',args.root/'results/cohort_frequency_v1/joint_analysis/analysis.json',args.root/'results/panel_repeat_composition.tsv',args.root/'results/complete_truth_v1/calls.tsv',args.root/'results/complete_truth_v1/events.tsv',args.root/'results/chm13_chr1_pilot/manifest.json',args.root/'results/chr1_call_error_v1/strict_metric_verification.tsv']},
        scope='Posthoc analysis of complete-truth outcomes; no new alignment or calling; original pilot rescored with current strict core. No R2/monophyly filtering. Spatial trends may reflect chromosome, SV-type, size and event-placement composition. Repeated carriers are not independent trials.')
    (out/'provenance.json').write_text(json.dumps(info,indent=2)+'\n')
    for table in (p,s):
        for r in table:
            if r['stratification']=='All' and r['position_measure'] in ('absolute_bin','relative_bin','cdr_point'):print(json.dumps(r))


if __name__=='__main__':main()
