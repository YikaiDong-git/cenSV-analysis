#!/usr/bin/env python3
"""Export complete-truth event/call metrics, negative-window rates, and figures."""
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys
import argparse
import os

import numpy as np
import pysam

from complete_truth_reference import CompositeReference
from complete_truth_scoring import score_complete_truth
from complete_truth_relaxed import score_relaxed
from parse_complete_truth_calls import parse_complete_calls
from realistic_sv_simulation import ROOT
from simulate_sv_recovery import write_json


def table(path,rows):
    rows=list(rows)
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w') as handle:
        writer=csv.DictWriter(handle,keys,delimiter='\t');writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v) if isinstance(v,(list,dict,tuple)) else v for k,v in row.items()})


def diagnostics(events,calls,rows):
    for event,row in zip(events,rows):
        expected=(event['start'],event['end'],event['alt'])
        candidate=[c for c in calls if c['chrom']==event['chrom'] and c.get('edit') and
                   abs(c['edit'][0]-event['start'])<=5000]
        if not candidate:continue
        closest=min(candidate,key=lambda c:(abs(c['edit'][0]-expected[0])+abs(c['edit'][1]-expected[1]),
                                            abs(c['length']-(len(expected[2])-expected[1]+expected[0])),str(c['id'])))
        a,b,alt=closest['edit']
        row.update(nearest_call_id=closest['id'],nearest_call_type=closest['kind'],
                   nearest_start_shift=a-event['start'],nearest_end_shift=b-event['end'],
                   nearest_net_length_error=closest['length']-(len(event['alt'])-event['end']+event['start']),
                   nearest_allele_interpretation=closest['allele_interpretation'])


def figures(metrics,events,out,suffix):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    scenarios=['snv_only','isolated','standing','compound','snv_only_fragmented','standing_fragmented']
    rates=sorted({r['array_rate'] for r in metrics})
    for ax,key,title in zip(axes,['sensitivity_lower','precision_lower','specificity_lower'],
                            ['Strict event sensitivity','Strict record precision','Active-array window specificity']):
        for scenario in scenarios:
            means=[];low=[];high=[]
            for rate in rates:
                rows=[r for r in metrics if r['scope']=='production_active' and r['scenario']==scenario and r['array_rate']==rate]
                if key=='specificity_lower':
                    values=[r.get('active_array_specificity_lower') for r in rows]
                else:values=[r.get(key) for r in rows]
                values=[x for x in values if x is not None]
                means.append(np.mean(values) if values else np.nan)
                low.append(min(values) if values else np.nan);high.append(max(values) if values else np.nan)
            ax.plot(np.asarray(rates)*100,means,'o-',label=scenario)
            ax.fill_between(np.asarray(rates)*100,low,high,alpha=.10)
        ax.set(xlabel='Nominal array/flank substitution rate (%)',ylabel='Fraction',title=title,ylim=(-.03,1.03))
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=3,frameon=False,fontsize=8)
    fig.suptitle('Complete sequence truth: means and seed ranges; unresolved matches kept separate',fontsize=11)
    fig.tight_layout(rect=(0,.10,1,.92))
    fig.savefig(out/f'metrics{suffix}.png',dpi=180);fig.savefig(out/f'metrics{suffix}.pdf');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,3.8))
    for ax,scope in zip(axes,['raw_full','production_active']):
        for scenario in ['snv_only','snv_only_fragmented']:
            for seed in sorted({r['seed'] for r in metrics}):
                group=sorted([r for r in metrics if r['scope']==scope and r['scenario']==scenario and r['seed']==seed],key=lambda r:r['array_rate'])
                if group:ax.plot([r['array_rate']*100 for r in group],[r['fp_calls'] for r in group],
                                  'o-' if scenario=='snv_only' else 's--',alpha=.65,label=f'{scenario}, seed {seed}')
        ax.set(xlabel='Nominal array/flank substitution rate (%)',ylabel='False SV records / query',title=scope.replace('_',' '))
    axes[1].legend(fontsize=6,frameon=False)
    fig.suptitle('Known zero-SV controls: all reported SV records are false under the simulation model',fontsize=10)
    fig.tight_layout(rect=(0,0,1,.93));fig.savefig(out/f'negative_controls{suffix}.png',dpi=180)
    fig.savefig(out/f'negative_controls{suffix}.pdf');plt.close(fig)
    labels=['panel_DEL_340','panel_INS_340','panel_DEL_680','panel_INS_680','panel_DEL_mid','panel_INS_long',
            'panel_DEL_tail','copy_340','copy_680','copy_200kb','inversion_10kb','dispersed_680']
    titles=['Panel DEL 339–341 bp','Panel INS 339–341 bp','Panel DEL 678–682 bp','Panel INS 678–682 bp',
            'Panel DEL 1–3 kb','Panel INS 5–15 kb','Panel DEL 50–200 kb','Tandem 340 bp','Tandem 680 bp',
            'Tandem 200 kb','Inversion 10 kb','Dispersed copy 680 bp']
    fig,axes=plt.subplots(1,3,figsize=(12,6.8),sharey=True)
    for ax,scenario in zip(axes,['isolated','standing','standing_fragmented']):
        matrix=np.full((len(labels),len(rates)),np.nan)
        for y,label in enumerate(labels):
            for x,rate in enumerate(rates):
                group=[e for e in events if e['scope']=='production_active' and e['scenario']==scenario and
                       e['array_rate']==rate and e['label']==label]
                if not group:continue
                tp=sum(e['recovered'] for e in group);matrix[y,x]=tp/len(group)
                ax.text(x,y,f'{tp}/{len(group)}',ha='center',va='center',fontsize=8,
                        color='white' if tp/len(group)>.65 else 'black')
        image=ax.imshow(matrix,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto')
        ax.set(title=scenario.replace('_',' '),xlabel='Nominal array/flank substitution rate (%)')
        ax.set_xticks(range(len(rates)),[f'{r*100:g}' for r in rates]);ax.set_yticks(range(len(labels)),titles)
    fig.suptitle('Focal event sensitivity by class: class-balanced scenarios, not a panel-wide accuracy estimate',fontsize=10)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(out/f'focal_classes{suffix}.png',dpi=180)
    fig.savefig(out/f'focal_classes{suffix}.pdf');plt.close(fig)
    accepted=[e for e in events if e.get('breakpoint_shift_bp') is not None and e.get('scope')=='production_active']
    if accepted:
        fig,axes=plt.subplots(1,2,figsize=(8,3))
        axes[0].hist([abs(e['breakpoint_shift_bp']) for e in accepted],bins=30)
        axes[0].set(xlabel='Absolute breakpoint displacement (bp)',ylabel='Certified events')
        axes[1].hist([e['size_error_bp'] for e in accepted if e.get('size_error_bp') is not None],bins=20)
        axes[1].set(xlabel='Size error (bp)',ylabel='Certified events')
        fig.tight_layout();fig.savefig(out/f'localization_size{suffix}.png',dpi=180)
        fig.savefig(out/f'localization_size{suffix}.pdf');plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',default='results/complete_truth_v1/manifest.json')
    p.add_argument('--allow-partial',action='store_true');a=p.parse_args()
    m=json.loads(Path(a.manifest).read_text());out=Path(a.manifest).resolve().parent
    metrics=[];event_rows=[];call_rows=[];window_rows=[];components=[];missing=[];relaxed_curves=[];near_events=[]
    for case in m['cases']:
        folder=out/'cases'/case['id']
        if not (folder/'complete.json').exists():missing.append(case['id']);continue
        bg=m['backgrounds'][case['background']];template=CompositeReference(bg['template'],m['reference'])
        calls=parse_complete_calls(folder/'svim/variants.vcf',template,folder/'svim/candidate_trace.json')
        with pysam.VariantFile(str(folder/'active.vcf.gz')) as vcf:retained={r.id for r in vcf}
        for scope,selected in [('raw_full',calls),('production_active',[c for c in calls if c['id'] in retained])]:
            result=score_complete_truth(template,case['truth'],selected,case['negative_windows'])
            relaxed=score_relaxed(template,case['truth'],selected,result)
            meta=dict(case=case['id'],background=case['background'],array_rate=bg['nominal_array_rate'],
                      seed=bg['seed'],scenario=case['scenario'],scope=scope)
            metric=dict(meta,**result['metrics'])
            chosen_curve=next(c for c in relaxed['curves'] if c['breakpoint_tolerance_bp']==500 and c['relative_dna_edit_threshold']==.01)
            metric.update(exploratory_near_recovery_fraction=chosen_curve['approximate_recovery_fraction'],
                          exploratory_additional_near_events=chosen_curve['additional_approximate_events'])
            metric['zero_SV_control_with_reported_calls']=(not case['truth'] and metric['n_calls']>0)
            metric['zero_SV_control_sequence_neutral_callset']=(not case['truth'] and metric['n_calls']>0 and
                                                              metric['full_sequence_reconstruction'] is True)
            windows=result['negative_windows']
            window_metadata={w['id']:w for w in case['negative_windows']}
            for window in windows:
                original=window_metadata[window['id']]
                window['stratum']=original['stratum']
                window['contains_input_contig_cut']=original.get('contains_input_contig_cut',False)
            for stratum in ['active_array','left_flank','right_flank','left_arm','right_arm']:
                group=[w for w in windows if w.get('stratum')==stratum and w.get('eligible_negative')]
                n=len(group);fp=sum(w['false_positive_window'] for w in group)
                uncertain=sum(w['unresolved_window'] and not w['false_positive_window'] for w in group)
                metric[stratum+'_negative_windows']=n
                metric[stratum+'_specificity_lower']=(n-fp-uncertain)/n if n else None
                metric[stratum+'_specificity_upper']=(n-fp)/n if n else None
            truth_lookup={e['id']:e for e in case['truth']};selected_lookup={c['id']:c for c in selected}
            for row in result['events']:
                truth=truth_lookup[row['id']]
                row.update(origin=truth['origin'],label=truth['label'],
                           representation_class='not_single_call_certified')
                if row.get('breakpoint_shift_bp') is not None:
                    shift=row['breakpoint_shift_bp']
                    row['representation_class']=('exact' if row.get('exact') else
                                                 'intrinsic_tandem_copy_placement' if truth['kind']=='TANDEM' and shift==-len(truth['alt'])
                                                 else 'other_sequence_equivalent')
                matched=row.get('matched_call_ids',[])
                if len(matched)==1 and truth.get('source_start') is not None:
                    trace=selected_lookup[matched[0]]['candidate_trace']
                    if trace.get('source_start') is not None:
                        row['source_start_shift']=trace['source_start']-truth['source_start']
                        row['source_end_shift']=trace['source_end']-truth['source_end']
            for row,parsed in zip(result['calls'],selected):
                edit=parsed.get('edit')
                row.update(allele_interpretation=parsed['allele_interpretation'],
                           reported_chrom=parsed['reported_chrom'],reported_pos=parsed['reported_pos'],
                           candidate_group=parsed['candidate_trace'].get('candidate_group'),
                           called_edit_size_bp=max(edit[1]-edit[0],len(edit[2])) if edit else None,
                           called_net_length_change_bp=parsed['length'],reported_svlen=parsed['reported_svlen'])
            single=[e for e in result['events'] if e.get('breakpoint_shift_bp') is not None]
            metric.update(single_call_certified_events=len(single),
                          n_panel_derived_truth=sum(e['origin']=='observed_panel_allele' for e in result['events']),
                          tp_panel_derived_events=sum(e['origin']=='observed_panel_allele' and e['recovered'] for e in result['events']),
                          single_call_type_compatible_fraction=(sum(e['type_compatible'] is True for e in single)/len(single) if single else None),
                          single_call_literal_type_fraction=(sum(e['literal_type_match'] is True for e in single)/len(single) if single else None),
                          single_call_within_100bp_fraction=(sum(e['localization_within_100bp'] is True for e in single)/len(single) if single else None),
                          intrinsic_tandem_placement_alias_events=sum(e['representation_class']=='intrinsic_tandem_copy_placement' for e in single),
                          fp_off_chr1_records=sum(c['status']=='fp' and c['chrom']!='chr1' for c in result['calls']))
            false_records=[c for c in result['calls'] if c['status']=='fp']
            metric['fp_records_size_10_to_49']=sum(c['called_edit_size_bp'] is not None and c['called_edit_size_bp']<50 for c in false_records)
            metric['fp_records_size_ge50']=sum(c['called_edit_size_bp'] is not None and c['called_edit_size_bp']>=50 for c in false_records)
            metric['fp_records_size_undefined']=sum(c['called_edit_size_bp'] is None for c in false_records)
            metrics.append(metric);event_rows.extend(dict(meta,**r) for r in result['events'])
            call_rows.extend(dict(meta,**r) for r in result['calls'])
            window_rows.extend(dict(meta,**r) for r in windows)
            components.extend(dict(meta,**r) for r in result['components'])
            relaxed_curves.extend(dict(meta,**r) for r in relaxed['curves'])
            near_events.extend(dict(meta,**r) for r in relaxed['events'])
            write_json(folder/f'evaluation_{scope}.json',result)
            write_json(folder/f'exploratory_closeness_{scope}.json',relaxed)
        template.close()
    if missing and not a.allow_partial:raise RuntimeError(f'{len(missing)} queries incomplete')
    suffix='_partial' if missing else ''
    for name,rows in [('metrics',metrics),('events',event_rows),('calls',call_rows),('negative_windows',window_rows),('components',components)]:
        table(out/f'{name}{suffix}.tsv',rows)
    table(out/f'exploratory_tolerance_curves{suffix}.tsv',relaxed_curves)
    table(out/f'exploratory_near_events{suffix}.tsv',near_events)
    grouped=defaultdict(list)
    for row in metrics:grouped[(row['array_rate'],row['scenario'],row['scope'])].append(row)
    aggregated=[]
    for (rate,scenario,scope),rows in grouped.items():
        item=dict(array_rate=rate,scenario=scenario,scope=scope,n_queries=len(rows))
        for key in ['n_truth','n_calls','tp_events','fn_events','unresolved_events','tp_calls','fp_calls','unresolved_calls']:
            item[key]=sum(r[key] for r in rows)
        for key in ['sensitivity_lower','sensitivity_upper','precision_lower','precision_upper','active_array_specificity_lower','active_array_specificity_upper','exploratory_near_recovery_fraction']:
            values=[r[key] for r in rows if r.get(key) is not None]
            item[key+'_mean']=float(np.mean(values)) if values else None
            item[key+'_min']=min(values) if values else None;item[key+'_max']=max(values) if values else None
        item['pooled_sensitivity_lower']=item['tp_events']/item['n_truth'] if item['n_truth'] else None
        item['pooled_sensitivity_upper']=(item['tp_events']+item['unresolved_events'])/item['n_truth'] if item['n_truth'] else None
        item['pooled_precision_lower']=item['tp_calls']/item['n_calls'] if item['n_calls'] else None
        item['pooled_precision_upper']=(item['tp_calls']+item['unresolved_calls'])/item['n_calls'] if item['n_calls'] else None
        aggregated.append(item)
    table(out/f'summary{suffix}.tsv',aggregated)
    report=dict(planned_queries=len(m['cases']),completed_queries=len(metrics)//2,missing=missing,
                planned_truth_events=sum(len(c['truth']) for c in m['cases']),
                scored_truth_events=sum(r['n_truth'] for r in metrics if r['scope']=='production_active'),
                evaluation_job_id=os.getenv('SLURM_JOB_ID'),
                manifest_sha256=hashlib.sha256(Path(a.manifest).read_bytes()).hexdigest(),
                truth_definition=m['truth_definition'],scope=m['scope'],
                uncertainty='Lower/upper are representation/evaluability bounds, not statistical confidence intervals',
                exploratory_tolerances='Post-calibration sensitivity analysis; primary strict labels unchanged. Main near-recovery display: 500 bp boundaries, 1% event-scaled DNA edit distance, 10% net-size error',
                independence=m['independence'],query_summaries=metrics)
    write_json(out/f'run_summary{suffix}.json',report)
    figures(metrics,event_rows,out,suffix)
    html_report(out,suffix,report,aggregated)
    source={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in (ROOT/'scripts').glob('*.py')}
    write_json(out/f'scoring_provenance{suffix}.json',dict(source_sha256=source,
                                                        evaluation_job_id=os.getenv('SLURM_JOB_ID'),completed=True,
                                                        manifest_sha256=hashlib.sha256(Path(a.manifest).read_bytes()).hexdigest()))
    print(json.dumps(dict(completed_queries=report['completed_queries'],missing=len(missing)),indent=2))
    if missing and os.getenv('SLURM_JOB_ID'):
        from recover_complete_truth import maybe_schedule_retry
        maybe_schedule_retry(m,out)


def html_report(out,suffix,report,summary):
    import html
    headers=['Array rate','Scenario','Queries','Mean strict sensitivity bounds','Exploratory near-recovery','Mean strict precision bounds','Mean active-array specificity bounds']
    body=[]
    fmt=lambda value:'undefined' if value is None else f'{100*value:.1f}%'
    for row in summary:
        if row['scope']!='production_active':continue
        pairs=[fmt(row[k+'_lower_mean'])+' – '+fmt(row[k+'_upper_mean'])
               for k in ['sensitivity','precision','active_array_specificity']]
        cells=[f"{row['array_rate']*100:.1f}%",row['scenario'],str(row['n_queries']),pairs[0],
               fmt(row['exploratory_near_recovery_fraction_mean']),pairs[1],pairs[2]]
        body.append('<tr>'+''.join('<td>'+html.escape(x)+'</td>' for x in cells)+'</tr>')
    status='COMPLETE' if not report['missing'] else 'PARTIAL'
    doc='''<!doctype html><meta charset="utf-8"><title>Complete-truth SV benchmark</title>
<style>body{font:16px system-ui;max-width:1150px;margin:35px auto;padding:0 20px;color:#18222e}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}img{max-width:100%}p{line-height:1.5}</style>
<h1>Complete-truth minimap2 / svim-asm benchmark</h1>'''
    doc+=f'<p><strong>{status}: {report["completed_queries"]}/{report["planned_queries"]} queries.</strong> '
    doc+='All simulated substitutions and structural transformations are logged and exactly replayable. '
    doc+='The table summarizes the production active-array filter; raw full-chromosome results are in the TSV files.</p>'
    doc+='<p>Lower and upper bounds reflect unresolved representations or bounded matching searches; '
    doc+='they are not confidence intervals. Rates are conditional on these scenarios, not estimates of reference-panel accuracy. '
    doc+='Precision depends on the planted-event density and class mix. '
    doc+='Exact known-SNP-restored structural transformations are scored; this is not literal genotyping of every background SNP. '
    doc+='SV-type agreement does not prove a uniquely identifiable biological mechanism.</p>'
    doc+='<p>Zero-SV controls contain substitutions but no simulated structural events. Reported SVs are false under that event definition; '
    doc+='a balanced set of such calls can nevertheless be sequence-neutral on the known background. '
    doc+='Full-sequence reconstruction is reported separately so event errors are not conflated with incorrect DNA reconstruction. '
    doc+='These results cannot by themselves label existing reference-panel calls biologically true or false.</p>'
    doc+='<table><thead><tr>'+''.join('<th>'+x+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join(body)+'</tbody></table>'
    doc+=f'<img src="metrics{suffix}.png"><img src="negative_controls{suffix}.png">'
    doc+=f'<img src="focal_classes{suffix}.png">'
    doc+='<p>Specificity counts fixed negative reference windows. Precision counts VCF records; one event can have multiple records. '
    doc+='High negative-window specificity can coexist with lower event precision when errors cluster around real SVs. '
    doc+='Standing SVs are included in truth, not subtracted using baseline calls. Contig splits are assembly boundaries, not biological SVs. '
    doc+='Full-versus-fragmented comparisons jointly change contiguity, query names, and query-length-dependent mapping ties. '
    doc+='The zero-divergence no-SV replicates are duplicate controls, not independent biological samples.</p>'
    doc+='<p>Strictly unmatched calls near true events may be partial or sequence-inexact representations, rather than invented events. '
    doc+='Constructed inversions and large copies are intentionally overrepresented relative to the observed catalog; '
    doc+='use the class-stratified results rather than interpreting the overall class-balanced sensitivity as panel sensitivity.</p>'
    doc+='<p>The exploratory near-recovery column permits 500 bp boundary displacement, 1% DNA edit distance scaled by true event size, '
    doc+='and 10% net-size error. It can combine split calls, but does not certify exact sequence recovery. '
    doc+='This tolerance analysis was added after calibration and is explicitly exploratory. '
    doc+='Full 100/500/1000 bp and 0/0.1/1/5% DNA-error curves and censored-distance/search bounds are in '
    doc+='<a href="exploratory_tolerance_curves'+suffix+'.tsv">tolerance curves</a> and '
    doc+='<a href="exploratory_near_events'+suffix+'.tsv">near-event diagnostics</a>.</p>'
    doc+='<p>Files: <a href="metrics'+suffix+'.tsv">query metrics</a>, <a href="events'+suffix+'.tsv">events</a>, '
    doc+='<a href="calls'+suffix+'.tsv">calls</a>, <a href="negative_windows'+suffix+'.tsv">negative windows</a>, '
    doc+='<a href="run_summary'+suffix+'.json">run summary</a>.</p>'
    doc+='<p>Compiled report: <a href="latex_report'+suffix+'/report.pdf">PDF</a> and '
    doc+='<a href="latex_report'+suffix+'/report.tex">LaTeX source</a>. '
    doc+='The scoring job builds these after exporting the tables and figures.</p>'
    (out/f'report{suffix}.html').write_text(doc)


if __name__=='__main__':main()
