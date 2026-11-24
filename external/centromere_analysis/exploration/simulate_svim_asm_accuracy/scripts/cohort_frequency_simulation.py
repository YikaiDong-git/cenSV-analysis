#!/usr/bin/env python3
"""Known-frequency chromosome cohort; per-query production calling and exact truth."""
import argparse
from collections import Counter
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
import pysam

from classify_sv_call_errors_v2 import classify_call_errors
from complete_truth_benchmark import apply_truth
from complete_truth_reference import CompositeReference
from evaluate_multichrom_sv import _write_json, validate_truth
from ensure_scipy_assignment import ensure_assignment_available
from multichrom_sv_simulation import mutate_background, read_active_bed, windows_for
from multichrom_truth_scoring import score_complete_truth
from native_truth import reverse_complement
from parse_complete_truth_calls import parse_complete_calls
from realistic_sv_simulation import BED, BIN, CALLER, INDEX, REF, ROOT


RECIPES = [
    ('alpha_299', 'INS', 299, 'Alu_A_alpha'),
    ('alpha_6087', 'INS', 6087, 'long_L1_B_alpha'),
    ('TE_299', 'INS', 299, 'Alu_A_TE'),
    ('TE_6087', 'INS', 6087, 'long_L1_B_TE'),
    ('deletion_340', 'DEL', 340, None),
    ('deletion_680', 'DEL', 680, None),
    ('deletion_6087', 'DEL', 6087, None),
    ('deletion_100000', 'DEL', 100000, None),
    ('tandem_340', 'TANDEM', 340, None),
    ('tandem_680', 'TANDEM', 680, None),
    ('tandem_6087', 'TANDEM', 6087, None),
    ('tandem_100000', 'TANDEM', 100000, None),
    ('donor_50000', 'INS', 50000, None),
    ('dispersed_680', 'DISPERSED', 680, None),
    ('inversion_10000', 'INV', 10000, None),
    ('replacement_680_340', 'REPLACEMENT', 680, None),
]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def sequence_sha(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def sample_name(index):
    return f'SIM{index//2:03d}_hap{index%2+1}'


def write_sequence(handle, name, sequence):
    handle.write('>' + name + '\n')
    for start in range(0, len(sequence), 1000000):
        handle.write(sequence[start:start+1000000] + '\n')


def carriers(count, n, seed, founders=8, controls_per_founder=2):
    assert n % founders == 0 and 0 < count <= n//2
    rng = np.random.default_rng(seed); allocations = np.full(founders, count//founders)
    allocations[rng.permutation(founders)[:count % founders]] += 1
    result = []
    for founder, k in enumerate(allocations):
        result.extend((rng.choice(n//founders-controls_per_founder, int(k), replace=False) + founder*(n//founders)).tolist())
    assert len(set(result)) == count
    return sorted(result)


def geometry(chrom, intervals, config, chrom_index):
    rng = np.random.default_rng(config['base_seed'] + chrom_index*100003)
    free = [[a+15000, b-15000] for a,b in intervals if b-a>30000]
    padding = config['event_separation_bp']
    def reserve(size):
        available = [(i,a,b) for i,(a,b) in enumerate(free) if b-a >= size]
        if not available:
            raise ValueError(f'{chrom}: no remaining space for {size}-bp support')
        weights = np.array([b-a-size+1 for _,a,b in available], dtype=float)
        i,a,b = available[int(rng.choice(len(available),p=weights/weights.sum()))]
        start = a if rng.integers(2)==0 else b-size
        left,right=free.pop(i)
        free.extend([[left,start-padding],[start+size+padding,right]])
        free[:]=[x for x in free if x[0]<x[1]]
        return start
    events = []
    for recipe_index,(recipe,kind,size,payload) in enumerate(RECIPES):
        for replicate in range(config['placements_per_recipe_per_chromosome']):
            group = chrom_index*2 + replicate
            count = config['carrier_counts'][(recipe_index + group) % len(config['carrier_counts'])]
            events.append(dict(id=f'{recipe}_site{replicate+1}',recipe=recipe,kind=kind,size=size,
                chrom=chrom,placement=replicate,payload_id=payload,carrier_count=count,
                carriers=carriers(count,config['n_haplotypes'],config['base_seed']+chrom_index*100003+recipe_index*103+replicate*17,
                                  controls_per_founder=config['guaranteed_controls_per_founder'])))
    for e in sorted(events,key=lambda x:(-x['size'],x['id'])):
        if e['kind']=='TANDEM':
            start=reserve(e['size']+1)
            e.update(start=start+e['size'],end=start+e['size'],source_start=start,source_end=start+e['size'])
        elif e['kind'] in ('DEL','INV','REPLACEMENT'):
            start=reserve(e['size']);e.update(start=start,end=start+e['size'])
        elif e['payload_id']:
            start=reserve(1);e.update(start=start,end=start)
        else:
            source=reserve(e['size']);start=reserve(1)
            e.update(start=start,end=start,source_start=source,source_end=source+e['size'])
    return sorted(events,key=lambda e:e['start'])


def prepare(config_path):
    config=json.loads(Path(config_path).read_text());out=ROOT/config['output_dir']
    assert config['n_haplotypes']==256 and len(config['chromosomes'])*2==len(config['carrier_counts'])
    assert len(config['founder_array_substitution_rates'])==8
    out.mkdir(parents=True,exist_ok=True)
    if (out/'manifest.json').exists():raise FileExistsError(out/'manifest.json')
    original=json.loads((ROOT/'results/multichrom_sv_v1/manifest.json').read_text())
    active=read_active_bed(BED)
    plans={chrom:dict(events=geometry(chrom,active[chrom],config,k),active_intervals=active[chrom])
           for k,chrom in enumerate(config['chromosomes'])}
    balance={recipe:Counter(e['carrier_count'] for p in plans.values() for e in p['events'] if e['recipe']==recipe)
             for recipe,_,_,_ in RECIPES}
    assert all(v==Counter(config['carrier_counts']) for v in balance.values())
    m=dict(config=config,plans=plans,reference=str(REF),active_bed=str(BED),index=str(INDEX),
        caller_source=str(CALLER),bin_dir=str(BIN),payloads=original['payloads'],
        caller_parameters=original['caller_parameters'],reference_lengths=original['reference_lengths'],
        alignment_arguments=['-a','-x','asm5','--cs','-r2k','-t',str(config['mapping_threads'])],
        identity=dict(config_sha256=sha(config_path),active_bed_sha256=sha(BED),
                      source_manifest_sha256=sha(ROOT/'results/multichrom_sv_v1/manifest.json'),
                      index_size=INDEX.stat().st_size,index_mtime_ns=INDEX.stat().st_mtime_ns),
        strict_scoring_parameters=dict(max_component_calls=48,max_component_truth=8,
                                      max_compound_comparisons=100000,max_subset_states=10000000),
        query_unit='Distinct constructed full-chromosome haplotype; paired into128 diploid samples',
        truth_frequency='Carrier counts fixed before alignment; each recipe occurs at every frequency once across eight distinct sites',
        founder_model='Eight shared SNP backbones,32 privately mutated haplotypes each; carrier counts balanced across founders',
        zero_sv_controls='Last two haplotypes in every founder block are reserved no-SV controls; additional zero-edit haplotypes retained',
        software=dict(minimap2=subprocess.check_output([str(BIN/'minimap2'),'--version'],text=True).strip(),
                      bcftools=subprocess.check_output([str(BIN/'bcftools'),'--version'],text=True).splitlines()[0]))
    _write_json(out/'manifest.json',m)
    print(json.dumps(dict(chromosomes=list(plans),haplotypes=config['n_haplotypes'],events=sum(len(p['events']) for p in plans.values()),
        queries=len(plans)*config['n_haplotypes'],event_carrier_occurrences=sum(e['carrier_count'] for p in plans.values() for e in p['events'])),indent=2))


def run_command(command,log,**kwargs):
    log.write('$ '+' '.join(map(str,command))+'\n');log.flush()
    subprocess.run(list(map(str,command)),check=True,stderr=log,**kwargs)


def ensure_founders(m,chrom,root):
    folder=root/'chromosomes'/chrom/'founders';folder.mkdir(parents=True,exist_ok=True)
    chrom_index=m['config']['chromosomes'].index(chrom)
    with (folder/'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (folder/'complete.json').exists():return
        with pysam.FastaFile(m['reference']) as ref:original=ref.fetch(chrom).upper()
        for i,rate in enumerate(m['config']['founder_array_substitution_rates']):
            target=folder/f'F{i}.fa.gz'
            c=dict(m['config'],outside_array_substitution_rate=m['config']['founder_arm_substitution_rate'])
            mutated,checks=mutate_background(original,chrom,rate,m['config']['base_seed']+chrom_index*100003+i*1009,
                m['plans'][chrom]['active_intervals'],c,folder/f'F{i}.snvs.tsv.gz')
            with tempfile.TemporaryDirectory(dir=os.getenv('TMPDIR','/tmp'),prefix='sv_founder_') as temp:
                fasta=Path(temp)/'founder.fa'
                with fasta.open('w') as f:write_sequence(f,chrom,mutated)
                pysam.tabix_compress(str(fasta),str(target),force=True)
            pysam.faidx(str(target))
            _write_json(folder/f'F{i}.json',dict(sequence_sha256=sequence_sha(mutated),**checks))
        _write_json(folder/'complete.json',dict(completed=True,founders=8))


def materialize(m,chrom,index,template):
    events=[]
    for original in m['plans'][chrom]['events']:
        if index not in original['carriers']:continue
        e={k:v for k,v in original.items() if k!='carriers'}
        a,b=e['start'],e['end'];kind=e['kind']
        if e.get('payload_id'):
            e['alt']=m['payloads'][e['payload_id']]['sequence']
        elif kind in ('INS','TANDEM','DISPERSED'):
            e['alt']=template[e['source_start']:e['source_end']]
        elif kind=='INV':e['alt']=reverse_complement(template[a:b])
        elif kind=='REPLACEMENT':e['alt']=m['payloads']['long_L1_B_alpha']['sequence'][:340]
        else:e['alt']=''
        e.update(deleted_sha256=sequence_sha(template[a:b]),inserted_sha256=sequence_sha(e['alt']),
            deleted_length=b-a,inserted_length=len(e['alt']),net_length_change=len(e['alt'])-(b-a),
            support_start=a,support_end=max(a+1,b))
        assert e['alt']!=template[a:b]
        assert not (set(e['alt'])-set('ACGT'))
        events.append(e)
    return events


def prepare_case(m,chrom,index,root,temp,query_handle,save_individual=False):
    case_dir=root/'chromosomes'/chrom/'haps'/f'H{index:03d}';case_dir.mkdir(parents=True,exist_ok=True)
    founder=index//(m['config']['n_haplotypes']//8)
    path=root/'chromosomes'/chrom/'founders'/f'F{founder}.fa.gz'
    expected=json.loads(path.with_name(f'F{founder}.json').read_text())['sequence_sha256']
    for attempt in range(3):
        try:
            with pysam.FastaFile(str(path)) as fa:base=fa.fetch(chrom).upper()
            if sequence_sha(base)!=expected:
                raise ValueError('Founder FASTA checksum differs from independently replayed source')
            break
        except (OSError,ValueError):
            if attempt==2:raise
            time.sleep(2)
    c=dict(m['config'],outside_array_substitution_rate=m['config']['private_arm_substitution_rate'])
    seed=m['config']['base_seed']+m['config']['chromosomes'].index(chrom)*1000003+index*10007+7919
    template,checks=mutate_background(base,chrom,m['config']['private_array_substitution_rate'],seed,
        m['plans'][chrom]['active_intervals'],c,case_dir/'private_snvs.tsv.gz')
    events=materialize(m,chrom,index,template);query,ledger=apply_truth(template,events)
    sample=sample_name(index);qname=sample+'_'+chrom
    template_path=temp/f'H{index:03d}.template.fa'
    with template_path.open('w') as f:write_sequence(f,chrom,template)
    pysam.faidx(str(template_path))
    write_sequence(query_handle,qname,query)
    individual=temp/f'H{index:03d}.query.fa'
    if save_individual:
        with individual.open('w') as f:write_sequence(f,qname,query)
    truth=dict(id=f'{chrom}_H{index:03d}',chrom=chrom,haplotype_index=index,sample_name=sample,
        query_name=qname,founder_index=founder,truth=ledger,template_sha256=sequence_sha(template),
        query_sha256=sequence_sha(query),query_length=len(query),private_seed=seed,
        negative_windows=windows_for(chrom,m['plans'][chrom]['active_intervals']),**checks)
    _write_json(case_dir/'truth.json',truth)
    oracle=validate_truth(template_path,m['reference'],ledger,truth['negative_windows'])
    assert oracle['passed'];_write_json(case_dir/'oracle.json',oracle)
    return dict(index=index,path=case_dir,template=template_path,query=individual,qname=qname)


def stage_index(m):
    scratch=Path('/scratch/jwein22/sv_cohort_frequency_cache')
    try:scratch.mkdir(parents=True,exist_ok=True)
    except OSError:return Path(m['index'])
    dest=scratch/'reference.asm5.mmi'
    with (scratch/'.index.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        metadata=scratch/'index_identity.json'
        expected=dict(source=m['index'],size=m['identity']['index_size'],mtime_ns=m['identity']['index_mtime_ns'])
        if not dest.exists() or dest.stat().st_size!=expected['size'] or not metadata.exists() or json.loads(metadata.read_text())!=expected:
            partial=scratch/f'index.{os.getpid()}.tmp';shutil.copyfile(m['index'],partial);os.replace(partial,dest)
            _write_json(metadata,expected)
    return dest


def map_batch(m,index,cases,query,temp,log):
    sam=temp/'batch.sam'
    with sam.open('wb') as output:
        run_command([BIN/'minimap2',*m['alignment_arguments'],index,query],log,stdout=output)
    byname={c['qname']:c for c in cases};writers={}
    with pysam.AlignmentFile(str(sam),'r') as alignments:
        assert dict(zip(alignments.references,alignments.lengths))==m['reference_lengths']
        for c in cases:
            c['unsorted']=temp/f"H{c['index']:03d}.unsorted.bam"
            writers[c['qname']]=pysam.AlignmentFile(str(c['unsorted']),'wb',header=alignments.header)
        try:
            for rec in alignments:
                assert rec.query_name in byname
                writers[rec.query_name].write(rec)
        finally:
            for w in writers.values():w.close()
    sam.unlink()


def call_case(m,case,log,*,output=None):
    folder=output or case['path'];folder.mkdir(parents=True,exist_ok=True)
    bam=folder/'alignment.bam'
    run_command([BIN/'samtools','sort','-@','3','-m','1G','-o',bam,case['unsorted']],log)
    run_command([BIN/'samtools','index','-@','3',bam],log)
    env=dict(os.environ,PYTHONPATH=str(CALLER),SVIM_ASM_SOURCE=str(CALLER))
    command=[sys.executable,ROOT/'scripts/trace_svim_asm.py','haploid',folder/'svim',bam,m['reference'],
             '--sample',sample_name(case['index'])]
    for key,value in m['caller_parameters'].items():command.extend(['--'+key,str(value)])
    run_command(command,log,env=env,stdout=log)
    run_command([BIN/'bcftools','view','-T',m['active_bed'],'-Oz','-o',folder/'active.vcf.gz',folder/'svim/variants.vcf'],log)
    run_command([BIN/'bcftools','index','-t',folder/'active.vcf.gz'],log)


def finish_case(m,case):
    folder=case['path'];truth=json.loads((folder/'truth.json').read_text())
    reference=CompositeReference(case['template'],m['reference'])
    try:
        calls=parse_complete_calls(folder/'active.vcf.gz',reference,folder/'svim/candidate_trace.json')
        result=score_complete_truth(reference,truth['truth'],calls,truth['negative_windows'],
            **m['strict_scoring_parameters'])
        result['metadata']=dict(scope='production_active',case_id=truth['id'],target_chrom=truth['chrom'])
        _write_json(folder/'scoring/production_active.json',result)
        _write_json(folder/'scoring/parsed_calls_production_active.json',calls)
        errors=classify_call_errors(reference,truth['truth'],calls,result)
        _write_json(folder/'scoring/call_errors.json',errors)
    finally:reference.close()
    with pysam.VariantFile(str(folder/'active.vcf.gz')) as vcf:
        with pysam.VariantFile(str(folder/'merge_input.vcf.gz'),'wz',header=vcf.header) as output:
            for rec in vcf:
                assert rec.id and ';' not in rec.id
                rec.id=f"H{case['index']:03d}__{rec.id}";output.write(rec)
    pysam.tabix_index(str(folder/'merge_input.vcf.gz'),preset='vcf',force=True)
    _write_json(folder/'complete.json',dict(completed=True,case_id=truth['id'],query_sha256=truth['query_sha256'],
        template_sha256=truth['template_sha256'],records=len(calls),metrics=result['metrics'],
        execution_code_identity_sha256=m.get('_execution_code_identity_sha256'),
        job=os.getenv('SLURM_JOB_ID'),node=os.uname().nodename))


def vcf_records(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else Path(path).open()) as f:
        return Counter(tuple(line.rstrip().split('\t')[:2]+line.rstrip().split('\t')[3:])
                       for line in f if not line.startswith('#'))


def run_chunk(manifest_path,task,pilot=False):
    m=json.loads(Path(manifest_path).read_text());root=Path(manifest_path).resolve().parent
    ensure_assignment_available()
    if pilot:
        chrom=m['config']['chromosomes'][1];indices=[0,30,64,224]
    else:
        chunks=m['config']['n_haplotypes']//m['config']['chunk_size']
        chrom=m['config']['chromosomes'][task//chunks]
        start=(task%chunks)*m['config']['chunk_size'];indices=list(range(start,start+m['config']['chunk_size']))
        gate=json.loads((root/'batch_equivalence.json').read_text());assert gate['passed']
        identity_path=root/'execution_code_identity.json'
        expected=json.loads(identity_path.read_text())
        assert all(sha(path)==value for path,value in expected.items()), 'Execution code changed after freeze'
        m['_execution_code_identity_sha256']=sha(identity_path)
    ensure_founders(m,chrom,root);index=stage_index(m)
    work=root/'logs';work.mkdir(exist_ok=True)
    with (work/('pilot.log' if pilot else f'task_{task:02d}.log')).open('a') as log:
        remaining=[i for i in indices if not (root/'chromosomes'/chrom/'haps'/f'H{i:03d}'/'complete.json').exists()]
        for offset in range(0,len(remaining),m['config']['query_batch_size']):
            batch=remaining[offset:offset+m['config']['query_batch_size']]
            with tempfile.TemporaryDirectory(dir=os.getenv('TMPDIR','/tmp'),prefix='sv_cohort_') as scratch:
                temp=Path(scratch);query=temp/'queries.fa'
                with query.open('w') as f:cases=[prepare_case(m,chrom,i,root,temp,f,pilot) for i in batch]
                print('Mapping',chrom,batch,flush=True);map_batch(m,index,cases,query,temp,log)
                for case in cases:
                    call_case(m,case,log)
                    if pilot:
                        single=temp/f"single_{case['index']}";single.mkdir()
                        sam=single/'individual.sam'
                        with sam.open('wb') as output:
                            run_command([BIN/'minimap2',*m['alignment_arguments'],index,case['query']],log,stdout=output)
                        with pysam.AlignmentFile(str(sam),'r') as src:
                            with pysam.AlignmentFile(str(single/'individual.bam'),'wb',header=src.header) as dst:
                                for rec in src:dst.write(rec)
                        individual=dict(case,unsorted=single/'individual.bam')
                        call_case(m,individual,log,output=single/'called')
                        assert vcf_records(case['path']/'active.vcf.gz')==vcf_records(single/'called/active.vcf.gz'), 'Batch/individual VCF disagreement'
                        assert vcf_records(case['path']/'svim/variants.vcf')==vcf_records(single/'called/svim/variants.vcf'), 'Batch/individual raw VCF disagreement'
                    finish_case(m,case)
                    print('Completed',chrom,case['index'],flush=True)
        if pilot:
            _write_json(root/'batch_equivalence.json',dict(passed=True,chrom=chrom,haplotypes=indices,
                comparison='Identical active-array VCF data fields except record IDs, individual versus batched mapping with same QNAME and flags',
                node=os.uname().nodename,job=os.getenv('SLURM_JOB_ID')))
        else:
            _write_json(root/'chunks'/f'{task:02d}.json',dict(completed=True,chrom=chrom,haplotypes=indices,
                job=os.getenv('SLURM_JOB_ID'),node=os.uname().nodename))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['prepare','pilot','run'])
    p.add_argument('--config',type=Path,default=ROOT/'config/cohort_frequency_v1.json')
    p.add_argument('--manifest',type=Path,default=ROOT/'results/cohort_frequency_v1/manifest.json')
    p.add_argument('--task',type=int,default=0)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.config)
    else:run_chunk(a.manifest,a.task,a.action=='pilot')
