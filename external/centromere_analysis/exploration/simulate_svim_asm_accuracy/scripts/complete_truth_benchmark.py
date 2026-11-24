#!/usr/bin/env python3
"""Complete sequence-truth scenarios for minimap2/svim-asm evaluation."""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import time

import numpy as np
import pysam

from native_truth import reverse_complement
from realistic_sv_simulation import BASE,BIN,REF,BED,INDEX,CALLER,ROOT,fasta_write
from score_sv_recovery import normalize_alleles
from simulate_sv_recovery import digest,run,write_json

CONFIG=ROOT/'config/complete_truth_v1.json'


def panel_pool():
    folder=BASE/'sv_extraction/decompose/mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200/chr1'
    path=folder/'chr1_variants.tsv';pool=[]
    csv.field_size_limit(sys.maxsize)
    with path.open() as handle:
        for row_index,row in enumerate(csv.DictReader(handle,delimiter='\t')):
            if row['svtype'] not in ['DEL','INS']:continue
            for allele_index,alt in enumerate(row['alt'].split(',')):
                if set(alt.upper()+row['ref'].upper())-set('ACGT'):continue
                start,end,allele=normalize_alleles(int(row['pos'])-1,row['ref'].upper(),alt.upper())
                if start==end and allele:kind='INS';size=len(allele)
                elif start<end and not allele:kind='DEL';size=end-start
                else:continue
                pool.append(dict(chrom='chr1',start=start,end=end,alt=allele,kind=kind,size=size,
                                 origin='observed_panel_allele',panel_row=row_index,panel_alt_index=allele_index,
                                 panel_id_sha256=digest(row['id'])))
    return pool,str(path)


def supports(event):
    result=[(event['start'],max(event['start']+1,event['end']))]
    if 'source_start' in event:result.append((event['source_start'],event['source_end']))
    return result


def available(candidate,events,padding=10000):
    return all(b+padding<=c or d+padding<=a for a,b in supports(candidate)
               for old in events for c,d in supports(old))


def choose_panel(pool,kind,low,high,events,rng):
    choices=[x for x in pool if x['kind']==kind and low<=x['size']<=high]
    for idx in rng.permutation(len(choices)):
        event=dict(choices[idx])
        if available(event,events):return event
    raise ValueError(f'No geometrically eligible panel allele: {kind} {low}-{high}')


def build_blueprint(pool,seed,start,end):
    rng=np.random.default_rng(seed);events=[]
    for label,kind,low,high in [('panel_DEL_tail','DEL',50000,200000),
                              ('panel_DEL_340','DEL',339,341),('panel_INS_340','INS',339,341),
                              ('panel_DEL_680','DEL',678,682),('panel_INS_680','INS',678,682),
                              ('panel_DEL_mid','DEL',1000,3000),('panel_INS_long','INS',5000,15000)]:
        event=choose_panel(pool,kind,low,high,events,rng)
        event['label']=label;events.append(event)
    for label,kind,size in [('copy_200kb','TANDEM',200000),('copy_340','TANDEM',340),
                            ('copy_680','TANDEM',680),('inversion_10kb','INV',10000),
                            ('dispersed_680','DISPERSED',680)]:
        for attempt in range(20000):
            p=int(rng.integers(start+size+20000,end-size-20000))
            event=dict(chrom='chr1',start=p,end=p+size if kind=='INV' else p,kind=kind,
                       size=size,label=label,origin='constructed_structural_control')
            if kind=='TANDEM':event.update(source_start=p-size,source_end=p)
            if kind=='DISPERSED':
                source=int(rng.integers(start+10000,end-size-10000))
                if abs(source-p)<500000:continue
                event.update(source_start=source,source_end=source+size)
            if available(event,events):events.append(event);break
        else:raise ValueError(f'Could not allocate {label}')
    standing=[]
    for i in range(6):
        event=choose_panel(pool,'DEL' if i%2==0 else 'INS',50,2000,events+standing,rng)
        event['label']=f'standing_{i}';standing.append(event)
    return dict(focal=events,standing=standing)


def mutate_template(reference,rate,seed,config,ledger_path):
    rng=np.random.default_rng(seed);data=bytearray(reference.encode());array=np.frombuffer(data,dtype=np.uint8)
    alphabet=np.frombuffer(b'ACGT',dtype=np.uint8);lut=np.full(256,-1,dtype=np.int16)
    lut[alphabet]=np.arange(4)
    start,end=config['array_start']-500000,config['array_end']+500000
    bin_size=config['cluster_bin_size'];weights=rng.lognormal(0,config['cluster_lognormal_sigma'],
                                                            (end-start+bin_size-1)//bin_size)
    weights/=weights.mean();n_mutated=0;regional=0
    with gzip.open(ledger_path,'wt') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['chrom','position_0based','ref','alt','zone'])
        blocks=[(a,min(a+bin_size,right),inside)
                for left,right,inside in [(0,start,False),(start,end,True),(end,len(data),False)]
                for a in range(left,right,bin_size)]
        for a,b,in_array in blocks:
            if in_array:
                weight=weights[min((max(a,start)-start)//bin_size,len(weights)-1)]
                local_rate=min(rate*weight,config['max_cluster_rate'])
            else:local_rate=config['outside_array_substitution_rate'] if rate else 0
            n=min(b-a,int(rng.poisson((b-a)*local_rate)))
            if not n:continue
            pos=np.sort(rng.choice(b-a,n,replace=False))+a
            codes=lut[array[pos]];valid=codes>=0;pos=pos[valid];codes=codes[valid]
            old=array[pos].copy();new=alphabet[(codes+rng.integers(1,4,len(pos)))%4]
            array[pos]=new
            for p,x,y in zip(pos,old,new):writer.writerow(['chr1',int(p),chr(int(x)),chr(int(y)),
                                                        'array_and_flanks' if in_array else 'chromosome_arms'])
            n_mutated+=len(pos);regional+=int(in_array)*len(pos)
    replay=bytearray(reference.encode());seen=set()
    with gzip.open(ledger_path,'rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            p=int(row['position_0based'])
            if p in seen or replay[p]!=ord(row['ref']) or row['ref']==row['alt']:
                raise ValueError('Invalid SNV truth ledger')
            replay[p]=ord(row['alt']);seen.add(p)
    if replay!=data:raise ValueError('SNV ledger fails exact reconstruction')
    return data.decode(),dict(snv_count=n_mutated,array_and_flanks_snv_count=regional,
                              reference_length=len(data),ledger_replay_verified=True)


def materialize_events(blueprint,scenario,sequence,pool,seed):
    if scenario=='snv_only':return []
    events=[dict(e) for e in blueprint['focal']]
    if scenario in ['standing','compound']:events.extend(dict(e) for e in blueprint['standing'])
    rng=np.random.default_rng(seed+10000)
    for event in events:
        if event['kind'] in ['TANDEM','DISPERSED']:
            event['alt']=sequence[event['source_start']:event['source_end']]
        elif event['kind']=='INV':event['alt']=reverse_complement(sequence[event['start']:event['end']])
        if scenario=='compound' and event['label'] in ['copy_340','copy_680','inversion_10kb']:
            old_size=event['size']
            donor=choose_panel(pool,'INS',678 if old_size==340 else 339,
                               682 if old_size==340 else 341,[],rng)
            event.update(original_kind=event['kind'],kind='REPLACEMENT',
                         end=event['start']+old_size,alt=donor['alt'],insert_panel_row=donor['panel_row'])
            event.pop('source_start',None);event.pop('source_end',None)
        event.update(id=event['label'],deleted_sha256=digest(sequence[event['start']:event['end']]),
                     inserted_sha256=digest(event['alt']),deleted_length=event['end']-event['start'],
                     inserted_length=len(event['alt']),net_length_change=len(event['alt'])-(event['end']-event['start']))
        event['support_start']=event['start'];event['support_end']=max(event['start']+1,event['end'])
    return sorted(events,key=lambda e:(e['start'],e['end'],e['id']))


def apply_truth(sequence,events):
    pieces=[];cursor=0;offset=0;ledger=[]
    for event in sorted(events,key=lambda e:(e['start'],e['end'])):
        a,b=event['start'],event['end']
        if not 0<=cursor<=a<=b<=len(sequence):raise ValueError('Overlapping or out-of-range truth')
        pieces.extend([sequence[cursor:a],event['alt']])
        entry=dict(event,query_start=a+offset,query_end=a+offset+len(event['alt']))
        ledger.append(entry);offset+=len(event['alt'])-(b-a);cursor=b
    pieces.append(sequence[cursor:]);query=''.join(pieces)
    back=[];cursor=0
    for e in ledger:
        back.extend([query[cursor:e['query_start']],sequence[e['start']:e['end']]])
        cursor=e['query_end']
    back.append(query[cursor:])
    if ''.join(back)!=sequence:raise ValueError('SV ledger cannot reconstruct baseline')
    return query,ledger


def write_input_query(path,query,fragments=None):
    if not fragments:
        fasta_write(path,'chr1',query);return
    cursor=0
    with Path(path).open('w') as handle:
        for fragment in fragments:
            a,b=fragment['query_start'],fragment['query_end']
            if a!=cursor or digest(query[a:b])!=fragment['sha256']:
                raise ValueError('Input fragmentation ledger changed DNA')
            handle.write('>'+fragment['name']+'\n')
            for i in range(a,b,80):handle.write(query[i:min(i+80,b)]+'\n')
            cursor=b
    if cursor!=len(query):raise ValueError('Input fragments omit sequence')


def add_fragmented_scenarios(manifest_path):
    path=Path(manifest_path);m=json.loads(path.read_text())
    if any(c.get('input_fragments') for c in m['cases']):raise ValueError('Fragmented scenarios already added')
    array_start,array_end=m['config']['array_start'],m['config']['array_end']
    additions=[]
    for case in m['cases']:
        case['input_architecture']='full_contig'
        if case['scenario'] not in ['snv_only','standing']:continue
        candidates=[(w['start']+w['end'])//2 for w in case['negative_windows'] if w['stratum']=='active_array']
        desired=[array_start+(array_end-array_start)*i//3 for i in [1,2]]
        cuts=sorted({min(candidates,key=lambda x:abs(x-target)) for target in desired})
        if len(cuts)!=2 or cuts[1]-cuts[0]<m['config'].get('min_middle_contig_length',100000):
            raise ValueError('No valid independent contig cuts')
        qcuts=[cut+sum(len(e['alt'])-(e['end']-e['start']) for e in case['truth'] if e['end']<=cut) for cut in cuts]
        bg=m['backgrounds'][case['background']]
        with pysam.FastaFile(bg['template']) as fa:sequence=fa.fetch('chr1').upper()
        query,_=apply_truth(sequence,case['truth'])
        if digest(query)!=case['query_sha256']:raise ValueError('Fragment control differs from paired full contig')
        new=copy.deepcopy(case);new['id']+='_fragmented';new['scenario']+='_fragmented'
        new.update(input_architecture='three_contigs',paired_full_contig_case=case['id'],reference_contig_cuts=cuts)
        bounds=[0]+qcuts+[len(query)]
        new['input_fragments']=[dict(name=f'chr1_part{i+1}',query_start=a,query_end=b,sha256=digest(query[a:b]))
                                for i,(a,b) in enumerate(zip(bounds,bounds[1:]))]
        for window in new['negative_windows']:
            window['contains_input_contig_cut']=any(window['start']<=cut<window['end'] for cut in cuts)
        additions.append(new)
    m['cases'].extend(additions)
    m['config']['assembly_scenarios']=['snv_only_fragmented','standing_fragmented']
    m['contiguity_convention']='Same known genomic sequence; input breaks add no biological SV truth; architecture and name/length hash effects are joint'
    (path.parent/'manifest_before_fragmentation.json').write_bytes(path.read_bytes())
    write_json(path,m);print('Added',len(additions),'paired fragmentation cases; total',len(m['cases']),flush=True)


def negative_windows(events,length,config):
    start,end=config['array_start'],config['array_end'];pad=config['negative_support_padding'];windows=[]
    intervals=[(max(0,a-pad),min(length,b+pad)) for e in events for a,b in supports(e)]
    for label,left,right,width in [('active_array',start,end,50000),('left_flank',start-500000,start,50000),
                                   ('right_flank',end,end+500000,50000),('left_arm',0,start-500000,1000000),
                                   ('right_arm',end+500000,length,1000000)]:
        for a in range(left,right,width):
            b=min(a+width,right)
            if any(a<d and b>c for c,d in intervals):continue
            windows.append(dict(id=f'{label}_{a}',chrom='chr1',start=a,end=b,stratum=label))
    return windows


def prepare(config_path):
    config=json.loads(Path(config_path).read_text());out=ROOT/config['output_dir'];out.mkdir(parents=True,exist_ok=True)
    if (out/'manifest.json').exists():raise FileExistsError(out/'manifest.json')
    pool,panel_path=panel_pool();blueprints={seed:build_blueprint(pool,seed,config['array_start'],config['array_end']) for seed in config['seeds']}
    with pysam.FastaFile(str(REF)) as fa:reference=fa.fetch('chr1').upper()
    cases=[];backgrounds={}
    for rate in config['array_substitution_rates']:
        for seed in config['seeds']:
            name=f'rate{int(rate*1000):03d}_seed{seed}';folder=out/'backgrounds'/name;folder.mkdir(parents=True,exist_ok=True)
            seq,metrics=mutate_template(reference,rate,seed,config,folder/'snvs.tsv.gz')
            fasta_write(folder/'template.fa','chr1',seq);pysam.faidx(str(folder/'template.fa'))
            bg=dict(id=name,template=str(folder/'template.fa'),nominal_array_rate=rate,seed=seed,
                    template_sha256=digest(seq),**metrics)
            backgrounds[name]=bg;write_json(folder/'metadata.json',bg)
            for scenario in config['scenarios']:
                events=materialize_events(blueprints[seed],scenario,seq,pool,seed)
                query,ledger=apply_truth(seq,events)
                case=dict(id=name+'_'+scenario,background=name,scenario=scenario,truth=ledger,
                          negative_windows=negative_windows(ledger,len(seq),config),
                          query_length=len(query),query_sha256=digest(query),
                          sv_ledger_replay_verified=True)
                cases.append(case)
                print(case['id'],len(ledger),case['query_length'],flush=True)
            del seq
    manifest=dict(config=config,reference=str(REF),reference_chr1_sha256=digest(reference),
                  panel_source=panel_path,backgrounds=backgrounds,cases=cases,
                  truth_definition='All SNVs and simulation-defined SV transformations are logged; no baseline-call subtraction',
                  scope='Class-balanced conditional panel-allele challenge plus constructed controls, not population sensitivity',
                  independence='Three seeds; paired divergence/scenarios; rate0 SNV-only controls have identical DNA')
    write_json(out/'manifest.json',manifest)


def execute(manifest_path,task):
    m=json.loads(Path(manifest_path).read_text());case=m['cases'][task];bg=m['backgrounds'][case['background']]
    out=Path(manifest_path).parent/'cases'/case['id'];out.mkdir(parents=True,exist_ok=True)
    if (out/'complete.json').exists():return
    started=time.time()
    with pysam.FastaFile(bg['template']) as fa:template=fa.fetch('chr1').upper()
    if digest(template)!=bg['template_sha256']:raise ValueError('Template changed')
    query,ledger=apply_truth(template,case['truth'])
    if digest(query)!=case['query_sha256']:raise ValueError('Query differs from recorded truth')
    write_json(out/'truth.json',case)
    with (out/'commands.log').open('w') as log,tempfile.TemporaryDirectory(prefix='complete_sv_',dir=os.getenv('TMPDIR','/tmp')) as tmp:
        tmp=Path(tmp);write_input_query(tmp/'query.fa',query,case.get('input_fragments'));del query,template
        with (tmp/'alignment.sam').open('wb') as handle:
            run([BIN/'minimap2','-a','-x','asm5','--cs','-r2k','-t','3',INDEX,tmp/'query.fa'],log,stdout=handle)
        run([BIN/'samtools','sort','-@','3','-m','1G','-T',tmp/'sort','-o',out/'alignment.bam',tmp/'alignment.sam'],log)
        run([BIN/'samtools','index','-@','3',out/'alignment.bam'],log)
        env=os.environ.copy();env['PYTHONPATH']=str(CALLER);env['SVIM_ASM_SOURCE']=str(CALLER)
        command=[sys.executable,ROOT/'scripts/trace_svim_asm.py','haploid',out/'svim',out/'alignment.bam',REF,'--sample',case['id']]
        parameters=json.loads((ROOT/'config/chm13_chr1_pilot.json').read_text())['caller_parameters']
        for k,v in parameters.items():command.extend(['--'+k,str(v)])
        run(command,log,env=env,stdout=log)
        run([BIN/'bcftools','view','-T',BED,'-Oz','-o',out/'active.vcf.gz',out/'svim/variants.vcf'],log)
        run([BIN/'bcftools','index','-t',out/'active.vcf.gz'],log)
    write_json(out/'complete.json',dict(elapsed_seconds=time.time()-started,node=os.uname().nodename,
                                       job=os.getenv('SLURM_JOB_ID'),task=os.getenv('SLURM_ARRAY_TASK_ID')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','augment','run'])
    p.add_argument('--config',default=str(CONFIG));p.add_argument('--manifest',default='results/complete_truth_v1/manifest.json')
    p.add_argument('--task',type=int);a=p.parse_args()
    if a.action=='prepare':prepare(a.config)
    elif a.action=='augment':add_fragmented_scenarios(a.manifest)
    else:execute(a.manifest,a.task)
