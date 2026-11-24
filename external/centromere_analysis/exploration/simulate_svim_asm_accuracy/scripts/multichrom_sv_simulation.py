#!/usr/bin/env python3
"""Chromosome-wide known-truth simulations with paired TE/alpha payloads."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np
import pysam

from complete_truth_benchmark import apply_truth, supports
from native_truth import reverse_complement
from realistic_sv_simulation import BASE, BED, BIN, CALLER, INDEX, REF, ROOT, fasta_write
from simulate_sv_recovery import digest, run, write_json


def merge_intervals(intervals):
    merged = []
    for a, b in sorted(intervals):
        if a >= b:
            raise ValueError('Empty active interval')
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return merged


def read_active_bed(path):
    intervals = {}
    for line in Path(path).read_text().splitlines():
        if line and not line.startswith('#'):
            chrom, a, b, *_ = line.split('\t')
            intervals.setdefault(chrom, []).append((int(a), int(b)))
    return {chrom: merge_intervals(rows) for chrom, rows in intervals.items()}


def blueprint(chrom, intervals, seed, config):
    rng = np.random.default_rng(seed)
    free = [[a + 15000, b - 15000] for a, b in intervals if b - a > 30000]
    pad = config['event_separation_bp']

    def reserve(size, avoid=None, min_distance=0):
        candidates = []
        for i, (a, b) in enumerate(free):
            ranges = [(a, b)] if avoid is None else [(a, min(b, avoid[0]-min_distance)),
                                                     (max(a, avoid[1]+min_distance), b)]
            candidates.extend((i, left, right) for left, right in ranges if right-left >= size)
        if not candidates:
            raise ValueError(f'{chrom}: insufficient predeclared space for {size} bp')
        weights = np.array([b-a-size+1 for _, a, b in candidates], dtype=float)
        i, a, b = candidates[int(rng.choice(len(candidates), p=weights/weights.sum()))]
        start = a if rng.integers(2) == 0 else b-size
        original_a, original_b = free.pop(i)
        free.extend([[original_a, start-pad], [start+size+pad, original_b]])
        free[:] = [x for x in free if x[0] < x[1]]
        return start

    large = max(b-a for a, b in intervals) >= config['large_pack_min_interval_bp']
    specifications = [
        ('paired_299', 'INS', 299, 'paired_payload'),
        ('paired_6087', 'INS', 6087, 'paired_payload'),
        ('donor_340', 'INS', 340, 'array_derived'),
        ('donor_680', 'INS', 680, 'array_derived'),
        ('deletion_340', 'DEL', 340, 'array_derived'),
        ('deletion_680', 'DEL', 680, 'array_derived'),
        ('deletion_6087', 'DEL', 6087, 'array_derived'),
        ('tandem_340', 'TANDEM', 340, 'array_derived'),
        ('tandem_680', 'TANDEM', 680, 'array_derived'),
        ('dispersed_680', 'DISPERSED', 680, 'array_derived'),
        ('inversion_10000', 'INV', 10000, 'array_derived'),
        ('replacement_680_340', 'REPLACEMENT', 680, 'array_derived')]
    if large:
        specifications += [('deletion_100000', 'DEL', 100000, 'array_derived'),
                           ('tandem_100000', 'TANDEM', 100000, 'array_derived'),
                           ('donor_50000', 'INS', 50000, 'array_derived')]
    result = []
    for label, kind, size, sequence_class in sorted(specifications, key=lambda x: -x[2]):
        event = dict(id=label, label=label, chrom=chrom, kind=kind, size=size,
                     sequence_class=sequence_class, origin='constructed_size_control')
        if kind == 'TANDEM':
            a = reserve(size + 1)
            event.update(start=a+size, end=a+size, source_start=a, source_end=a+size)
        elif kind in ['DEL', 'INV', 'REPLACEMENT']:
            a = reserve(size)
            event.update(start=a, end=a+size)
        elif label.startswith('paired_'):
            a = reserve(1)
            event.update(start=a, end=a, pair_id='Alu_A' if size == 299 else 'long_L1_B',
                         origin='observed_matched_payload')
        else:
            source = reserve(size)
            a = reserve(1, (source, source+size), 30000) if kind == 'DISPERSED' else reserve(1)
            event.update(start=a, end=a, source_start=source, source_end=source+size)
        result.append(event)
    return dict(events=sorted(result, key=lambda e: e['start']), large_pack_eligible=large,
                large_pack_omission_reason=None if large else 'largest_active_interval_below_900000_bp')


def windows_for(chrom, intervals):
    return [dict(id=f'{chrom}_{a}', chrom=chrom, start=a, end=min(a+50000, right),
                 stratum='active_array', full_50kb=min(a+50000, right)-a == 50000)
            for left, right in intervals for a in range(left, right, 50000)]


def prepare(config_path):
    config = json.loads(Path(config_path).read_text())
    out = ROOT / config['output_dir']; out.mkdir(parents=True, exist_ok=True)
    if (out/'manifest.json').exists():
        raise FileExistsError(out/'manifest.json')
    intervals = read_active_bed(BED)
    with pysam.FastaFile(str(REF)) as reference:
        lengths = dict(zip(reference.references, reference.lengths))
    donors_path = ROOT/'results/native_te_v2/manifest.json'
    donors = json.loads(donors_path.read_text())['payloads']
    payloads = {p['id']: p for p in donors if p['pair_id'] in config['payload_pairs']}
    assert len(payloads) == 4
    for p in payloads.values():
        assert digest(p['sequence']) == p['sequence_sha256'] and p['rm_fraction'] >= .8
    plans = {}
    bundles = []
    for i, chrom in enumerate(config['chromosomes']):
        seed = config['base_seed'] + i*1009
        plans[chrom] = blueprint(chrom, intervals[chrom], seed, config)
        for rate in config['array_substitution_rates']:
            bundles.append(dict(id=f'{chrom}_rate{int(rate*1000):03d}', chrom=chrom,
                                seed=seed, array_rate=rate, active_intervals=intervals[chrom],
                                chromosome_length=lengths[chrom], simulations=['control', 'alpha_mixed', 'TE_mixed']))
    caller = json.loads((ROOT/'config/chm13_chr1_pilot.json').read_text())['caller_parameters']
    manifest = dict(config=config, reference=str(REF), reference_lengths=lengths,
                    active_bed=str(BED), index=str(INDEX), caller_source=str(CALLER),
                    caller_parameters=caller, alignment_arguments=['-a', '-x', 'asm5', '--cs', '-r2k', '-t', '3'],
                    blueprint=plans, payloads=payloads, bundles=bundles,
                    truth_definition='All constructed structural edits and background substitutions are known; no baseline-call subtraction',
                    primary_type_comparison='alpha_mixed plus control; shared events are not pooled again from TE_mixed',
                    paired_payload_comparison='Only paired_299 and paired_6087 at same positions, sizes and backgrounds',
                    selection='Geometry and fixed seed only; no alignment outcome used',
                    provenance=dict(active_bed_sha256=digest(BED.read_text()),
                                    payload_manifest_sha256=digest(donors_path.read_text()),
                                    config_sha256=digest(Path(config_path).read_text())))
    write_json(out/'manifest.json', manifest)
    print(json.dumps(dict(bundles=len(bundles), simulations=len(bundles)*3,
                          large_omitted=[c for c, v in plans.items() if not v['large_pack_eligible']]), indent=2))


def mutate_background(sequence, chrom, rate, seed, intervals, config, ledger):
    rng = np.random.default_rng(seed)
    data = bytearray(sequence.encode()); array = np.frombuffer(data, dtype=np.uint8)
    alphabet = np.frombuffer(b'ACGT', dtype=np.uint8)
    lookup = np.full(256, -1, dtype=np.int16); lookup[alphabet] = np.arange(4)
    regions = merge_intervals([(max(0, a-500000), min(len(data), b+500000)) for a, b in intervals])
    edges = sorted({0, len(data), *(x for pair in regions for x in pair)})
    width = config['cluster_bin_size']
    bins = [(a, min(a+width, right), any(left >= x and right <= y for x, y in regions))
            for left, right in zip(edges, edges[1:]) for a in range(left, right, width)]
    weights = rng.lognormal(0, config['cluster_lognormal_sigma'], sum(r for _, _, r in bins))
    weights /= weights.mean()
    count = inside = weight_index = 0
    with gzip.open(ledger, 'wt') as handle:
        writer = csv.writer(handle, delimiter='\t'); writer.writerow(['chrom', 'position_0based', 'ref', 'alt', 'zone'])
        for a, b, regional in bins:
            local = min(rate*weights[weight_index], config['max_cluster_rate']) if regional else (config['outside_array_substitution_rate'] if rate else 0)
            if regional:
                weight_index += 1
            n = min(b-a, int(rng.poisson((b-a)*local)))
            if not n:
                continue
            positions = np.sort(rng.choice(b-a, n, replace=False)) + a
            codes = lookup[array[positions]]; valid = codes >= 0; positions, codes = positions[valid], codes[valid]
            old = array[positions].copy(); new = alphabet[(codes+rng.integers(1, 4, len(positions))) % 4]
            array[positions] = new
            for pos, x, y in zip(positions, old, new):
                writer.writerow([chrom, int(pos), chr(int(x)), chr(int(y)), 'array_and_flanks' if regional else 'arms'])
            count += len(positions); inside += int(regional)*len(positions)
    replay = bytearray(sequence.encode()); seen = set()
    with gzip.open(ledger, 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            pos = int(row['position_0based'])
            assert row['chrom'] == chrom and pos not in seen and replay[pos] == ord(row['ref']) and row['ref'] != row['alt']
            replay[pos] = ord(row['alt']); seen.add(pos)
    assert replay == data
    return data.decode(), dict(snv_count=count, array_and_flanks_snv_count=inside, ledger_replay_verified=True)


def materialize(plan, sequence, scenario, payloads):
    if scenario == 'control':
        return []
    events = []
    for original in plan['events']:
        event = dict(original); a, b, kind = event['start'], event['end'], event['kind']
        if 'pair_id' in event:
            p = payloads[event['pair_id'] + ('_TE' if scenario == 'TE_mixed' else '_alpha')]
            event.update(alt=p['sequence'], sequence_class='TE_rich' if scenario == 'TE_mixed' else 'alpha_rich',
                         payload_id=p['id'], repeat_fraction=p['rm_fraction'], payload_sha256=p['sequence_sha256'])
        elif kind in ['TANDEM', 'DISPERSED', 'INS']:
            event['alt'] = sequence[event['source_start']:event['source_end']]
        elif kind == 'INV':
            event['alt'] = reverse_complement(sequence[a:b])
        elif kind == 'REPLACEMENT':
            event['alt'] = payloads['long_L1_B_alpha']['sequence'][:340]
        else:
            event['alt'] = ''
        assert not (set(event['alt']) - set('ACGT'))
        assert event['alt'] != sequence[a:b]
        event.update(deleted_sha256=digest(sequence[a:b]), inserted_sha256=digest(event['alt']),
                     deleted_length=b-a, inserted_length=len(event['alt']), net_length_change=len(event['alt'])-(b-a),
                     support_start=a, support_end=max(a+1, b))
        events.append(event)
    return events


def execute(manifest_path, task, *, align=True, diagnostics=False):
    m = json.loads(Path(manifest_path).read_text()); bundle = m['bundles'][task]
    out = Path(manifest_path).resolve().parent/'bundles'/bundle['id']; out.mkdir(parents=True, exist_ok=True)
    chrom, config = bundle['chrom'], m['config']
    if (out/'complete.json').exists() and align:
        return
    from evaluate_multichrom_sv import score_simulation, validate_truth
    started = time.time()
    background = out/'template.fa.gz'
    if not background.exists():
        with pysam.FastaFile(m['reference']) as fa:
            original = fa.fetch(chrom).upper()
        mutated, metrics = mutate_background(original, chrom, bundle['array_rate'], bundle['seed'],
                                             bundle['active_intervals'], config, out/'snvs.tsv.gz')
        with tempfile.TemporaryDirectory(prefix='sv_bg_', dir=os.getenv('TMPDIR', '/tmp')) as temp:
            source = Path(temp)/'template.fa'; fasta_write(source, chrom, mutated)
            pysam.tabix_compress(str(source), str(background), force=True)
        pysam.faidx(str(background))
        write_json(out/'background.json', dict(chrom=chrom, template_sha256=digest(mutated),
                   reference_chromosome_sha256=digest(original), **metrics))
        del original, mutated
    metadata = json.loads((out/'background.json').read_text())
    for scenario in bundle['simulations']:
        case_dir = out/scenario; case_dir.mkdir(exist_ok=True)
        if align and (case_dir/'complete.json').exists():
            continue
        with pysam.FastaFile(str(background)) as fa:
            template = fa.fetch(chrom).upper()
        assert digest(template) == metadata['template_sha256']
        events = materialize(m['blueprint'][chrom], template, scenario, m['payloads'])
        query, ledger = apply_truth(template, events)
        case = dict(bundle, id=bundle['id']+'_'+scenario, scenario=scenario, truth=ledger,
                    query_sha256=digest(query), query_length=len(query),
                    template_sha256=metadata['template_sha256'], negative_windows=windows_for(chrom, bundle['active_intervals']))
        write_json(case_dir/'truth.json', case)
        oracle = validate_truth(background, m['reference'], ledger, case['negative_windows'])
        assert oracle['passed']
        write_json(case_dir/'oracle.json', oracle)
        if align:
            with (case_dir/'commands.log').open('a') as log, tempfile.TemporaryDirectory(prefix='sv_mc_', dir=os.getenv('TMPDIR', '/tmp')) as temp:
                temp = Path(temp); fasta_write(temp/'query.fa', chrom, query); del query, template
                with (temp/'alignment.sam').open('wb') as handle:
                    run([BIN/'minimap2', *m['alignment_arguments'], m['index'], temp/'query.fa'], log, stdout=handle)
                run([BIN/'samtools', 'sort', '-@', '3', '-m', '1G', '-T', temp/'sort', '-o', case_dir/'alignment.bam', temp/'alignment.sam'], log)
                run([BIN/'samtools', 'index', '-@', '3', case_dir/'alignment.bam'], log)
                with pysam.AlignmentFile(str(case_dir/'alignment.bam')) as bam:
                    assert dict(zip(bam.references, bam.lengths)) == m['reference_lengths'], 'Index does not match the complete reference'
                env = os.environ.copy(); env['PYTHONPATH'] = str(CALLER); env['SVIM_ASM_SOURCE'] = str(CALLER)
                command = [sys.executable, ROOT/'scripts/trace_svim_asm.py', 'haploid', case_dir/'svim', case_dir/'alignment.bam', m['reference'], '--sample', case['id']]
                for key, value in m['caller_parameters'].items():
                    command += ['--'+key, str(value)]
                run(command, log, env=env, stdout=log)
                run([BIN/'bcftools', 'view', '-T', m['active_bed'], '-Oz', '-o', case_dir/'active.vcf.gz', case_dir/'svim/variants.vcf'], log)
                run([BIN/'bcftools', 'index', '-t', case_dir/'active.vcf.gz'], log)
        else:
            del query, template
        score_simulation(case_dir, background, m['reference'], chrom, diagnostics=diagnostics)
        marker = 'complete.json' if align else 'score_complete.json'
        write_json(case_dir/marker, dict(completed=True, node=os.uname().nodename,
                   job=os.getenv('SLURM_JOB_ID'), elapsed_bundle_seconds=time.time()-started))
        print('Completed', case['id'], flush=True)
    write_json(out/('complete.json' if align else 'score_complete.json'), dict(completed=True, node=os.uname().nodename,
               job=os.getenv('SLURM_JOB_ID'), elapsed_seconds=time.time()-started))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run', 'score'])
    parser.add_argument('--config', type=Path, default=ROOT/'config/multichrom_sv_v1.json')
    parser.add_argument('--manifest', type=Path, default=ROOT/'results/multichrom_sv_v1/manifest.json')
    parser.add_argument('--task', type=int)
    parser.add_argument('--diagnostics', action='store_true')
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare(args.config)
    else:
        execute(args.manifest, args.task, align=args.action == 'run', diagnostics=args.diagnostics)


if __name__ == '__main__':
    main()
