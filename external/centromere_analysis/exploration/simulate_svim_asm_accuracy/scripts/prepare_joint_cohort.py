#!/usr/bin/env python3
"""Combine frozen per-query outputs into four-chromosome haplotype callsets."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import gzip
import json
import os
import socket
from pathlib import Path
import subprocess
import tempfile

import pysam

from cohort_frequency_joint_reference import JointSNVOverlayReference, joint_truth, verify_joint_code_identity
from cohort_frequency_metrics import sample_name
from evaluate_multichrom_sv import _write_json
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BCFTOOLS = CENSV_ROOT + '/.pixi/envs/default/bin/bcftools'


def rewrite_record_ids(source, output, prefix, include_header):
    for line in source:
        if line.startswith('#'):
            if include_header:
                output.write(line)
            continue
        columns = line.rstrip('\n').split('\t')
        if len(columns) < 8 or not columns[2] or columns[2] == '.':
            raise ValueError('VCF record requires a source ID and at least eight fields')
        columns[2] = prefix + columns[2]
        output.write('\t'.join(columns) + '\n')


def prepare_one(args):
    from score_cohort_joint_truth import score_joint_truth
    root, h = Path(args[0]), args[1]
    manifest = json.loads((root / 'manifest.json').read_text())
    identity_sha256 = verify_joint_code_identity(root)
    folder = root / 'joint_cohort/haps' / f'H{h:03d}'
    marker = folder / 'complete.json'
    if marker.exists():
        if json.loads(marker.read_text()).get('joint_execution_code_identity_sha256') != identity_sha256:
            raise ValueError(f'Joint completion belongs to another code version: {marker}')
        return h, 'already_complete'
    (folder / 'scoring').mkdir(parents=True, exist_ok=True)
    truth = joint_truth(root, manifest, h)
    parsed, source_hashes = [], {}
    with tempfile.TemporaryDirectory(prefix=f'joint_H{h:03d}_') as temp:
        unsorted = Path(temp) / 'union.vcf'
        first = root / 'chromosomes' / manifest['config']['chromosomes'][0] / 'haps' / f'H{h:03d}/active.vcf.gz'
        with pysam.VariantFile(str(first)) as source:
            header = source.header.copy()
        if list(header.samples) != [sample_name(h)]:
            raise ValueError(f'Incorrect first source sample H{h:03d}')
        observed = set()
        with unsorted.open('w') as output:
            output.write(str(header))
            for chrom in manifest['config']['chromosomes']:
                source_folder = root / 'chromosomes' / chrom / 'haps' / f'H{h:03d}'
                source_vcf = source_folder / 'active.vcf.gz'
                source_parsed = source_folder / 'scoring/parsed_calls_production_active.json'
                calls = json.loads(source_parsed.read_text())
                call_ids = {c['id'] for c in calls}
                input_ids = set()
                with pysam.VariantFile(str(source_vcf)) as source:
                    if list(source.header.samples) != [sample_name(h)]:
                        raise ValueError(f'Incorrect source sample: {source_vcf}')
                    for record in source:
                        original = record.id
                        if original not in call_ids or original in input_ids:
                            raise ValueError(f'Invalid source record identity: {source_vcf}: {original}')
                        input_ids.add(original)
                        identity = f'{chrom}__{original}'
                        if identity in observed:
                            raise ValueError('Duplicate joint record ID')
                        observed.add(identity)
                with gzip.open(source_vcf, 'rt') as serialized:
                    rewrite_record_ids(serialized, output, prefix=f'{chrom}__', include_header=False)
                if input_ids != call_ids or len(call_ids) != len(calls):
                    raise ValueError(f'Source VCF/parsed coverage differs: {source_vcf}')
                for call in calls:
                    parsed.append(dict(call, id=f"{chrom}__{call['id']}", source_query_chrom=chrom,
                                       source_record_id=call['id']))
                source_hashes[chrom] = dict(vcf_sha256=hashlib.sha256(source_vcf.read_bytes()).hexdigest(),
                    parsed_sha256=hashlib.sha256(source_parsed.read_bytes()).hexdigest())
        subprocess.run([BCFTOOLS, 'sort', '-Oz', '-o', str(folder / 'active.vcf.gz'), str(unsorted)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run([BCFTOOLS, 'index', '-f', '-t', str(folder / 'active.vcf.gz')], check=True)
    reference = JointSNVOverlayReference(root, manifest, h)
    try:
        identities = reference.verify_templates()
        score = score_joint_truth(reference, truth['truth'], parsed, truth['negative_windows'],
                                  **manifest['strict_scoring_parameters'])
    finally:
        reference.close()
    score['metadata'] = dict(scope='joint_production_active', case_id=truth['id'],
        query_chromosomes=manifest['config']['chromosomes'], source_union=True,
        joint_execution_code_identity_sha256=identity_sha256,
        template_sha256_by_chromosome=identities)
    _write_json(folder / 'truth.json', truth)
    _write_json(folder / 'scoring/production_active.json', score)
    _write_json(folder / 'scoring/parsed_calls_production_active.json', parsed)
    with tempfile.TemporaryDirectory(prefix=f'merge_H{h:03d}_') as temp:
        text_input = Path(temp) / 'merge_input.vcf'
        with gzip.open(folder / 'active.vcf.gz', 'rt') as source, text_input.open('w') as output:
            rewrite_record_ids(source, output, prefix=f'H{h:03d}__', include_header=True)
        pysam.tabix_compress(str(text_input), str(folder / 'merge_input.vcf.gz'), force=True)
    pysam.tabix_index(str(folder / 'merge_input.vcf.gz'), preset='vcf', force=True)
    _write_json(marker, dict(completed=True, haplotype_index=h, records=len(parsed), metrics=score['metrics'],
        source_hashes=source_hashes, template_sha256_by_chromosome=identities,
        manifest_sha256=hashlib.sha256((root / 'manifest.json').read_bytes()).hexdigest(),
        joint_execution_code_identity_sha256=identity_sha256,
        job=os.environ.get('SLURM_JOB_ID'), node=socket.gethostname(),
        requested_cpus=os.environ.get('SLURM_CPUS_PER_TASK')))
    return h, 'completed'


def prepare(root, workers=4):
    from ensure_scipy_assignment import ensure_assignment_available
    ensure_assignment_available()
    root = Path(root).resolve()
    manifest = json.loads((root / 'manifest.json').read_text())
    if not 1 <= workers <= 4:
        raise ValueError('Joint preparation requires one to four worker CPUs')
    missing = [f'{chrom}/H{h:03d}' for chrom in manifest['config']['chromosomes']
               for h in range(manifest['config']['n_haplotypes'])
               if not (root / 'chromosomes' / chrom / 'haps' / f'H{h:03d}/complete.json').exists()]
    if missing:
        raise ValueError(f'Incomplete source cohort: {len(missing)} cases; examples={missing[:5]}')
    tasks = [(str(root), h) for h in range(manifest['config']['n_haplotypes'])]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for h, status in pool.map(prepare_one, tasks):
            print(f'Joint source H{h:03d} {status}', flush=True)
    _write_json(root / 'joint_cohort/complete.json', dict(completed=True, haplotypes=len(tasks),
        source_queries=len(tasks)*len(manifest['config']['chromosomes']),
        method='Substantive-record-preserving union of independent haploid query callsets',
        job=os.environ.get('SLURM_JOB_ID'), node=socket.gethostname(),
        requested_cpus=os.environ.get('SLURM_CPUS_PER_TASK')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default='results/cohort_frequency_v1')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    prepare(args.root, args.workers)
