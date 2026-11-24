"""One targeted resource/transient retry; never retry a scientific/software failure."""
import json
import os
from pathlib import Path
import subprocess

from simulate_sv_recovery import write_json


def retryable_indices(manifest,out,states):
    allowed={'TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL'}
    return [i for i,c in enumerate(manifest['cases'])
            if not (out/'cases'/c['id']/'complete.json').exists() and states.get(i) in allowed]


def maybe_schedule_retry(manifest,out):
    receipt=out/'submission.json';record=json.loads(receipt.read_text())
    if os.getenv('SLURM_JOB_ID')!=record.get('remaining',{}).get('score_job') or 'resource_retry' in record:
        return
    array=record['remaining']['query_job']
    missing=[i for i,c in enumerate(manifest['cases']) if not (out/'cases'/c['id']/'complete.json').exists()]
    if not missing:return
    result=subprocess.run(['sacct','-n','-X','-j',array,'--format=JobID%40,State%30','-P'],capture_output=True,text=True)
    states={}
    for line in result.stdout.splitlines():
        values=line.split('|')
        if len(values)>1 and values[0].startswith(array+'_'):
            try:states[int(values[0].split('_')[-1])]=values[1].split()[0].rstrip('+')
            except ValueError:pass
    errors={str(i):[dict(path=str(p),tail=p.read_text(errors='replace')[-2000:])
                    for p in (out/'slurm_logs').glob(f'*_{array}_{i}.err')] for i in missing}
    for i in missing:
        if states.get(i) not in [None,'FAILED']:continue
        message='\n'.join(e['tail'] for e in errors[str(i)]).lower()
        if any(token in message for token in ['oom-kill','oom_kill','out of memory','memoryerror']):
            states[i]='OUT_OF_MEMORY'
        elif 'time limit' in message:states[i]='TIMEOUT'
    selected=retryable_indices(manifest,out,states)
    write_json(out/'failure_audit.json',dict(missing_indices=missing,states=states,
                                            retryable_indices=selected,slurm_errors=errors,
                                            accounting_error=result.stderr if result.returncode else None))
    if not selected:return
    if record['planned_total_jobs']+len(selected)+1>150:raise RuntimeError('Retry would exceed job limit')
    from submit_complete_truth import submit
    job=submit(['--job-name=svtruth_resource_retry','--array='+','.join(map(str,selected))+'%8'],
               'run',memory_gb=32,time_limit_minutes=60)
    record['resource_retry']=dict(query_job=job,indices=selected,memory_gb=32,time_limit_minutes=60)
    record['planned_total_jobs']+=len(selected)+1;write_json(receipt,record)
    score=submit(['--job-name=svtruth_score_retry','--dependency=afterany:'+job],'score','--allow-partial')
    record['resource_retry']['score_job']=score;write_json(receipt,record)
