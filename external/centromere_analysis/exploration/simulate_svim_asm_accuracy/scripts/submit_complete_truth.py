#!/usr/bin/env python3
"""Submit calibration and then the remaining fully logged SV challenge matrix."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from realistic_sv_simulation import ROOT,CALLER
from simulate_sv_recovery import write_json

OUT=ROOT/'results/complete_truth_v1'


def submit(extra,action,*args,small=False,memory_gb=None,time_limit_minutes=None):
    memory_gb=memory_gb if memory_gb is not None else 2 if small else 24
    duration=time_limit_minutes if time_limit_minutes is not None else 10 if small else 30 if action=='run' else 60
    nodes=json.loads((OUT/'manifest.json').read_text())['config']['node_allowlist']
    command=['sbatch','--parsable','--partition=nodes','--nodes=1','--ntasks=1',
             '--nodelist='+','.join(nodes),
             '--cpus-per-task='+('1' if small else '4'),f'--mem={memory_gb}G',
             f'--time={duration//60:02d}:{duration%60:02d}:00',f'--chdir={ROOT}',
             f'--output={OUT}/slurm_logs/%x_%A_%a.out',f'--error={OUT}/slurm_logs/%x_%A_%a.err']
    command+=extra+['scripts/run_complete_truth_task.sh',action]+list(args)
    result=subprocess.check_output(command,text=True).strip().split(';')[0]
    print(action,result,flush=True);return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['calibration','remaining']);a=parser.parse_args()
    manifest=OUT/'manifest.json';m=json.loads(manifest.read_text())
    config=m['config'];validation=json.loads((OUT/'oracle_validation.json').read_text())
    if not validation['passed'] or validation['manifest_sha256']!=hashlib.sha256(manifest.read_bytes()).hexdigest():
        raise RuntimeError('Final manifest has not passed oracle validation')
    if validation['scorer_sha256']!=hashlib.sha256((ROOT/'scripts/complete_truth_scoring.py').read_bytes()).hexdigest():
        raise RuntimeError('Scorer changed after oracle validation')
    if config['cpus_per_job']>4 or len(config['node_allowlist'])>4 or len(m['cases'])+6>150:
        raise ValueError('Workflow exceeds authorized limits')
    if set(config['node_allowlist'])-{'node03','node04','node08','node12'}:
        raise ValueError('Use the same four nodes as existing workflow')
    receipt=OUT/'submission.json'
    r=json.loads(receipt.read_text()) if receipt.exists() else dict(preparation_jobs=['4789923','4789943'],
        nodes=config['node_allowlist'],cpus_per_query=4,max_concurrent_queries=8,
        planned_total_jobs=len(m['cases'])+5)
    if a.stage in r and r[a.stage].get('score_job'):raise RuntimeError('Stage already submitted')
    seed=config['seeds'][0]
    calibration=[i for i,c in enumerate(m['cases']) if m['backgrounds'][c['background']]['seed']==seed and
                 m['backgrounds'][c['background']]['nominal_array_rate'] in [0,.01] and
                 c['scenario'] in ['snv_only','isolated']]
    if len(calibration)!=4:raise ValueError('Expected four implementation calibration queries')
    selected=calibration if a.stage=='calibration' else [i for i in range(len(m['cases'])) if i not in calibration]
    if a.stage=='remaining':
        summary=json.loads((OUT/'run_summary_partial.json').read_text())
        scoring=json.loads((OUT/'scoring_provenance_partial.json').read_text())
        if (not scoring.get('completed') or scoring.get('evaluation_job_id')!=r['calibration']['score_job'] or
            summary.get('evaluation_job_id')!=r['calibration']['score_job'] or
            scoring['manifest_sha256']!=validation['manifest_sha256']):
            raise RuntimeError('Calibration scoring receipt is stale or incomplete')
        if not all((OUT/'cases'/m['cases'][i]['id']/'complete.json').exists() for i in calibration):
            raise RuntimeError('A calibration job failed or is incomplete')
        rows=[x for x in summary['query_summaries'] if x['scope']=='raw_full' and
              x['case'] in {m['cases'][i]['id'] for i in calibration}]
        if len(rows)!=4:raise RuntimeError('Calibration evaluation missing')
        for row in rows:
            if row['tp_events']+row['fn_events']+row['unresolved_events']!=row['n_truth']:
                raise RuntimeError('Inconsistent event accounting')
            if row['tp_calls']+row['fp_calls']+row['unresolved_calls']!=row['n_calls']:
                raise RuntimeError('Inconsistent call accounting')
            for key in ['sensitivity_lower','sensitivity_upper','precision_lower','precision_upper','specificity_lower','specificity_upper']:
                if row[key] is not None and not 0<=row[key]<=1:raise ValueError('Metric out of range')
        write_json(OUT/'calibration_gate.json',dict(passed=True,criterion='Completion, truth oracle, and accounting consistency; no performance threshold',
                                                  queries=4))
    if a.stage in r:
        job=r[a.stage]['query_job'];print('Reusing query array',job,flush=True)
    else:
        job=submit(['--job-name=svtruth_'+a.stage,'--array='+','.join(map(str,selected))+'%8'],'run')
        r[a.stage]=dict(query_job=job,indices=selected);write_json(receipt,r)
    score=submit(['--job-name=svtruth_score','--dependency=afterany:'+job],'score','--allow-partial')
    r[a.stage]['score_job']=score;write_json(receipt,r)
    if a.stage=='calibration':
        continuation=submit(['--job-name=svtruth_continue','--dependency=afterany:'+score],'continue',small=True)
        r['continuation_job']=continuation;write_json(receipt,r)
    sources={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in (ROOT/'scripts').glob('*.py')}
    caller_sources={str(f.relative_to(CALLER)):hashlib.sha256(f.read_bytes()).hexdigest() for f in CALLER.rglob('*.py')}
    write_json(OUT/f'provenance_{a.stage}.json',dict(source_sha256=sources,caller_source_sha256=caller_sources,
        manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),reference_resolved=str(Path(m['reference']).resolve())))


if __name__=='__main__':main()
