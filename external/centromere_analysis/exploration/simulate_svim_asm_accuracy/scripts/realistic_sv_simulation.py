#!/usr/bin/env python3
"""Native HPRC single-edit benchmark with paired full-alignment context."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import pysam

from simulate_sv_recovery import mutate, run, write_json, digest
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BASE = Path(CENSV_ROOT + '')
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/native_hprc_v1'
BIN = BASE / '.pixi/envs/default/bin'
REF = BASE / 'resources/chm13v2.0.fa.gz'
BED = BASE / 'resources/chm13v2.0_censat_v2.1.active.bed'
INDEX = ROOT / 'results/chm13_chr1_pilot/reference.asm5.mmi'
START, END = 121796048, 126300487
BACKGROUNDS = ['CHM13', 'HG00099_hap2', 'HG002_hap1', 'HG002_hap2', 'HG00099_hap1']
CALLER = BASE / 'software/svim-asm/src'


def read_fasta(path):
    with pysam.FastaFile(str(path)) as fa:
        name = fa.references[0]
        return name, fa.fetch(name).upper()


def fasta_write(path, name, sequence):
    with Path(path).open('w') as handle:
        handle.write('>' + name + '\n')
        for i in range(0, len(sequence), 80):
            handle.write(sequence[i:i+80] + '\n')


def cigar_blocks(rec, length):
    q, r = 0, rec.reference_start
    for op, n in rec.cigartuples or []:
        if op in (0, 7, 8):
            a, b = max(r, START), min(r+n, END)
            if a < b:
                qa, qb = q+a-r, q+b-r
                if rec.is_reverse:
                    qa, qb = length-qb, length-qa
                yield dict(qstart=qa, qend=qb, rstart=a, rend=b, reverse=rec.is_reverse)
            q += n
            r += n
        elif op in (1, 4, 5):
            q += n
        elif op in (2, 3):
            r += n


def select_sites(blocks):
    stable = max(blocks, key=lambda b: b['qend']-b['qstart'])
    stable_pos = (stable['qstart']+stable['qend'])//2
    edge = min(blocks, key=lambda b: abs(b['rstart']-(START+10000)))
    edge_r = min(max(START+10000, edge['rstart']), edge['rend']-1)
    edge_q = (edge['qend']-(edge_r-edge['rstart']) if edge['reverse'] else
              edge['qstart']+edge_r-edge['rstart'])
    merged = []
    for b in sorted(blocks, key=lambda b:b['qstart']):
        if merged and b['qstart'] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b['qend'])
        else:
            merged.append([b['qstart'], b['qend']])
    gaps = [(b[0]-a[1], a[1], b[0]) for a,b in zip(merged, merged[1:])]
    gap_details = {}
    if gaps:
        _, a, b = max(gaps)
        complex_q = (a+b)//2
        lb=max((x for x in blocks if x['qend']<=a),key=lambda x:x['qend'])
        rb=min((x for x in blocks if x['qstart']>=b),key=lambda x:x['qstart'])
        gap_details=dict(native_gap_start=a,native_gap_end=b,native_gap_length=b-a,
                         left_block=lb,right_block=rb,
                         orientation_consistent=lb['reverse']==rb['reverse'],
                         physical_array_membership='unverified_candidate')
    else:
        complex_q = (merged[0][0]+merged[-1][1])//2
    near = min(blocks, key=lambda b:min(abs(b['qstart']-complex_q), abs(b['qend']-complex_q)))
    return [dict(context='longest_aligned_block', position=stable_pos,
                 reverse=stable['reverse'],rough_reference_position=(stable['rstart']+stable['rend'])//2),
            dict(context='largest_unaligned_gap', position=complex_q,
                 reverse=near['reverse'],gap_details=gap_details,
                 rough_reference_position=(near['rstart']+near['rend'])//2),
            dict(context='reference_array_edge_candidate', position=edge_q,
                 reverse=edge['reverse'],rough_reference_position=edge_r)]


def prepare_background(task):
    name = BACKGROUNDS[task]
    folder = OUT / 'backgrounds' / name
    folder.mkdir(parents=True, exist_ok=True)
    if (folder/'metadata.json').exists():
        return
    if name == 'CHM13':
        with pysam.FastaFile(str(REF)) as fa:
            sequence = fa.fetch('chr1').upper()
        contig = 'chr1'
        sites = [dict(context=c, position=p, rough_reference_position=p) for c,p in
                 [('longest_aligned_block',124048267), ('largest_unaligned_gap',125000000),
                  ('reference_array_edge_candidate',START+10000)]]
        with pysam.AlignmentFile(str(ROOT/'results/chm13_chr1_pilot/cases/control/alignment.bam')) as bam:
            with pysam.AlignmentFile(str(folder/'unchanged.bam'),'wb',header=bam.header):
                pass
        metadata = dict(background=name, contig=contig, source_fasta=str(REF), sites=sites)
    else:
        profile = json.loads((ROOT/'results/hprc_background_profiles'/f'{name}.json').read_text())
        scores = Counter()
        for r in profile['records']:
            if r['mapq'] >= 20:
                scores[r['contig']] += min(r['ref_end'],END)-max(r['ref_start'],START)
        contig = scores.most_common(1)[0][0]
        fasta = BASE/'input/HPRC_bgzip'/f'{name}.fa.gz'
        with pysam.FastaFile(str(fasta)) as fa:
            sequence = fa.fetch(contig).upper()
        blocks = []
        bam_path = BASE/'alignments/svim_asm'/f'{name}_to_chm13.bam'
        with pysam.AlignmentFile(str(bam_path)) as bam:
            header = bam.header.to_dict()
            pg = [x for x in header.get('PG',[]) if x.get('PN')=='minimap2']
            if not pg or pg[0]['VN'] != '2.30-r1287' or '-r2k' not in pg[0]['CL']:
                raise ValueError(f'Unexpected upstream aligner provenance: {pg}')
            with pysam.AlignmentFile(str(folder/'unchanged.bam'),'wb',header=bam.header) as dest:
                for rec in bam.fetch(until_eof=True):
                    if rec.query_name != contig:
                        dest.write(rec)
                    elif (not rec.is_secondary and not rec.is_unmapped and rec.mapping_quality>=20
                          and rec.reference_name=='chr1' and rec.reference_start<END and rec.reference_end>START):
                        blocks.extend(cigar_blocks(rec,len(sequence)))
        sites = select_sites(blocks)
        metadata = dict(background=name, contig=contig, source_fasta=str(fasta),
                        source_bam=str(bam_path), upstream_pg=pg, sites=sites,
                        profile=profile['local_cs_counts'],
                        selection='Dominant array-overlapping contig; sites selected before injection outcomes',
                        observed_mapq20_reference_coverage=profile['mapq20_reference_coverage'])
    fasta_write(folder/'contig.fa',contig,sequence)
    pysam.faidx(str(folder/'contig.fa'))
    metadata.update(contig_length=len(sequence), contig_sha256=digest(sequence),
                    coordinate_convention='Native forward query coordinates, zero-based half-open')
    write_json(folder/'metadata.json',metadata)
    print(name,contig,len(sequence),sites,flush=True)


def plan():
    manifest = OUT/'manifest.json'
    if manifest.exists():
        raise FileExistsError(manifest)
    for name in BACKGROUNDS:
        path=OUT/'backgrounds'/name/'metadata.json'
        bg=json.loads(path.read_text())
        if name=='CHM13':
            for site in bg['sites']:
                site['reverse']=False
                if site['context']=='largest_unaligned_gap':site['context']='interior_control'
        else:
            blocks=[]
            with pysam.AlignmentFile(bg['source_bam']) as bam:
                for rec in bam.fetch(until_eof=True):
                    if (rec.query_name==bg['contig'] and not rec.is_secondary and not rec.is_unmapped
                        and rec.mapping_quality>=20 and rec.reference_name=='chr1'
                        and rec.reference_start<END and rec.reference_end>START):
                        blocks.extend(cigar_blocks(rec,bg['contig_length']))
            bg['sites']=select_sites(blocks)
        write_json(path,bg)
    donor = json.loads((OUT/'backgrounds/HG00099_hap1/metadata.json').read_text())
    with pysam.FastaFile(str(OUT/'backgrounds/HG00099_hap1/contig.fa')) as fa:
        center = donor['sites'][0]['position']
        donor_seq = fa.fetch(donor['contig'],center,center+10000).upper()
    backgrounds, events = {}, []
    for name in BACKGROUNDS[:4]:
        metadata = json.loads((OUT/'backgrounds'/name/'metadata.json').read_text())
        backgrounds[name] = metadata
        events.append(dict(id=name+'_control', background=name, kind='CONTROL', size=0,
                           context='control', position=metadata['sites'][0]['position'],
                           rough_reference_position=metadata['sites'][0]['rough_reference_position'], seed=20261029))
        for site in metadata['sites']:
            for kind in ['DEL','TANDEM','DONOR_INS']:
                for size in [340,680,10000]:
                    event = dict(site, id=f"{name}_{site['context']}_{kind.lower()}_{size}",
                                 background=name, kind=kind, size=size, seed=20261029)
                    if kind=='DONOR_INS':
                        from native_truth import reverse_complement
                        flip=donor['sites'][0]['reverse'] != site['reverse']
                        allele=reverse_complement(donor_seq[:size]) if flip else donor_seq[:size]
                        event.update(allele=allele, donor_background='HG00099_hap1',
                                     donor_contig=donor['contig'],donor_start=center,
                                     donor_end=center+size,donor_reverse_complemented=flip,
                                     donor_local_reverse=donor['sites'][0]['reverse'])
                    events.append(event)
    write_json(manifest,dict(backgrounds=backgrounds, events=events, reference=str(REF),
                             index=str(INDEX), donor=donor,
                             sequence_truth='Native query edit; reference correspondence independently assessed per site',
                             alignment_context='Full input alignment context; unchanged query-contig alignments reused identically',
                             limitations=['Background selection from four indexed haplotypes, not cohort-wide strata',
                                          'Local correspondence uniqueness is bounded, not genome-wide']))
    print(f'Planned {len(events)} queries',flush=True)


def edit_query(sequence,event):
    if event['kind']=='DONOR_INS':
        p=event['position']; allele=event['allele']
        return sequence[:p]+allele+sequence[p:],allele
    return mutate(sequence,event)


def execute(task,replicate=False):
    manifest=json.loads((OUT/'manifest.json').read_text())
    event=dict(manifest['events'][task])
    if replicate:
        if event['kind']!='CONTROL':raise ValueError('Only controls may be replicated')
        event['id']+='_replicate'
    bg=manifest['backgrounds'][event['background']]
    case=OUT/'cases'/event['id']; case.mkdir(parents=True,exist_ok=True)
    if (case/'complete.json').exists():
        return
    started=time.time()
    contig,sequence=read_fasta(OUT/'backgrounds'/event['background']/'contig.fa')
    if digest(sequence)!=bg['contig_sha256']:
        raise ValueError('Baseline input changed')
    changed,allele=edit_query(sequence,event)
    truth=dict(event,contig=contig,allele=allele,baseline_sha256=digest(sequence),
               edited_sha256=digest(changed),allele_sha256=digest(allele))
    if event['kind']=='DONOR_INS':
        truth['donor_equals_adjacent_native_copy']=(allele==sequence[event['position']-event['size']:event['position']])
    write_json(case/'truth.json',truth)
    with (case/'commands.log').open('w') as log, tempfile.TemporaryDirectory(prefix='native_sv_',dir=os.getenv('TMPDIR','/tmp')) as tmp:
        tmp=Path(tmp)
        fasta_write(tmp/'query.fa',contig,changed)
        del sequence,changed
        with (tmp/'alignment.sam').open('wb') as handle:
            run([BIN/'minimap2','-a','-x','asm5','--cs','-r2k','-t','3',INDEX,tmp/'query.fa'],log,stdout=handle)
        run([BIN/'samtools','sort','-@','3','-m','1G','-T',tmp/'sort','-o',case/'edited_contig.bam',tmp/'alignment.sam'],log)
        run([BIN/'samtools','index','-@','3',case/'edited_contig.bam'],log)
        run([BIN/'samtools','merge','-@','3','-o',tmp/'full_context.bam',
             OUT/'backgrounds'/event['background']/'unchanged.bam',case/'edited_contig.bam'],log)
        run([BIN/'samtools','index','-@','3',tmp/'full_context.bam'],log)
        params=json.loads((ROOT/'config/chm13_chr1_pilot.json').read_text())['caller_parameters']
        env=os.environ.copy();env['PYTHONPATH']=str(CALLER)
        command=[sys.executable,CALLER/'svim_asm/svim-asm','haploid',case/'svim',tmp/'full_context.bam',REF,'--sample',event['id']]
        for key,value in params.items():
            command.extend(['--'+key,str(value)])
        run(command,log,env=env,stdout=log)
        run([BIN/'bcftools','view','-T',BED,'-Oz','-o',case/'active.vcf.gz',case/'svim/variants.vcf'],log)
        run([BIN/'bcftools','index','-t',case/'active.vcf.gz'],log)
    write_json(case/'complete.json',dict(elapsed_seconds=time.time()-started,
                                        node=os.uname().nodename,slurm_job=os.getenv('SLURM_JOB_ID'),
                                        array_task=os.getenv('SLURM_ARRAY_TASK_ID'),completed=time.time()))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prep','plan','run','replicate'])
    parser.add_argument('--task',type=int)
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.action=='prep': prepare_background(args.task)
    elif args.action=='plan': plan()
    else: execute(args.task,replicate=args.action=='replicate')
