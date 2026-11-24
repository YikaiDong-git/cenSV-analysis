#!/usr/bin/env python3
"""Annotate changed sequences, keeping direct evidence separate from reference context."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UPSTREAM = Path(CENSV_ROOT + '')
MEI = UPSTREAM / 'discovery/2026-09-14_mei_annotation_centromeric_sv/output/T1_mei'
ATLAS = UPSTREAM / 'discovery/2026-09-23_monomer_atlas/output/hormon_chm13v2/chr1'
TOOLS = UPSTREAM.parent
RM_ENV = TOOLS / 'tools/CDR-Finder/env_rm/.pixi/envs/default'
SD_ENV = TOOLS / '_shared/tools/hor_tools/.pixi/envs/default'


def digest(sequence):
    return hashlib.sha256(sequence.upper().encode()).hexdigest()


def tool_aliases(hashes):
    hashes = list(hashes)
    aliases = {sha[:32]: sha for sha in hashes}
    if len(aliases) != len(hashes):
        raise ValueError('Tool sequence alias collision')
    return aliases


def union_bp(intervals):
    total, last = 0, -1
    for start, end in sorted(intervals):
        total += max(0, end - max(start, last))
        last = max(last, end)
    return total


def repeat_group(name, family):
    if name.replace('\\', '/').lower() == 'alr/alpha':
        return 'alpha_satellite'
    if family.startswith('Satellite'):
        return 'other_satellite'
    if family.split('/')[0] in {'LINE', 'SINE', 'LTR', 'DNA', 'RC', 'Retroposon'}:
        return 'transposable_element'
    if family in {'Simple_repeat', 'Low_complexity'}:
        return family.lower()
    return 'other_repeat'


def parse_rm(path):
    with Path(path).open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 15 or not fields[0].isdigit():
                continue
            yield {'query': fields[4], 'start': int(fields[5]) - 1,
                   'end': int(fields[6]), 'repeat': fields[9], 'family': fields[10],
                   'repeat_id': fields[14], 'group': repeat_group(fields[9], fields[10])}


def coverage(hits, length):
    by_group, by_family = defaultdict(list), defaultdict(list)
    for hit in hits:
        start, end = max(0, hit['start']), min(length, hit['end'])
        if end <= start:
            continue
        by_group[hit['group']].append((start, end))
        by_family[hit['family']].append((start, end))
    return {'length': length, 'annotated_bp': union_bp([iv for x in by_group.values() for iv in x]),
            'group_bp': {k: union_bp(v) for k, v in sorted(by_group.items())},
            'family_bp': {k: union_bp(v) for k, v in sorted(by_family.items())}}


def classify_component(length, rm, sd):
    if rm is not None:
        groups = rm['group_bp']
        alpha = groups.get('alpha_satellite', 0) / length
        te = groups.get('transposable_element', 0)
        if te:
            return 'mixed_TE_satellite' if groups.get('alpha_satellite', 0) or groups.get('other_satellite', 0) else 'TE_supported'
        if alpha >= .8:
            return 'alpha_satellite_RM'
        if rm['annotated_bp']:
            return 'other_or_partial_RM'
        return 'unmasked_RM'
    if sd and sd['high_identity_bp'] / length >= .8:
        return 'alpha_like_SD_only'
    return 'unresolved'


def parse_sd(path, lengths, aliases=None):
    blocks = defaultdict(list)
    with Path(path).open() as handle:
        for line in handle:
            f = line.rstrip().split('\t')
            if len(f) < 5 or not f[2].isdigit():
                continue
            name, monomer = f[:2]
            name = (aliases or {}).get(name, name)
            start, end, identity = int(f[2]), int(f[3]) + 1, float(f[4])
            if name not in lengths or not (0 <= start < end <= lengths[name]):
                raise ValueError('StringDecomposer interval outside its actual allele')
            blocks[name].append((start, end, identity, not monomer.startswith('?'), monomer))
    result = {}
    for name, rows in blocks.items():
        good = [(s, e) for s, e, ident, reliable, _ in rows if reliable and ident >= 80]
        denom = sum(e-s for s, e, _, _, _ in rows)
        result[name] = {'blocks': len(rows), 'high_identity_bp': union_bp(good),
                        'weighted_identity': sum((e-s)*ident for s, e, ident, _, _ in rows) / denom,
                        'monomer_labels': sorted({m for _, _, _, _, m in rows})}
    return result


def read_fasta_selected(path, wanted):
    name, pieces = None, []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name in wanted:
                    yield name, ''.join(pieces).upper()
                name, pieces = line[1:].split()[0], []
            elif name in wanted:
                pieces.append(line.strip())
        if name in wanted:
            yield name, ''.join(pieces).upper()


def atomic_json(path, obj):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, indent=2) + '\n')
    temp.replace(path)


def file_digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def prepare(out):
    import pysam
    manifest_path = out / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    work = out / 'repeat_annotation'; work.mkdir(exist_ok=True)
    sequences, events, literal_hashes = {}, [], set()
    handles = {key: pysam.FastaFile(value['template']) for key, value in manifest['backgrounds'].items()}
    for case in manifest['cases']:
        for event in case['truth']:
            row = {'case': case['id'], 'id': event['id'], 'kind': event['kind'],
                   'background': case['background'], 'components': [], 'reference_intervals': []}
            components = [('inserted', event['alt']),
                          ('removed', handles[case['background']].fetch(event['chrom'], event['start'], event['end']).upper())]
            for role, sequence in components:
                if not sequence:
                    continue
                sha = digest(sequence)
                sequences.setdefault(sha, sequence)
                row['components'].append({'role': role, 'sha256': sha, 'length': len(sequence)})
                if event['kind'] == 'INS' and role == 'inserted':
                    literal_hashes.add(sha)
            row['reference_intervals'].append({'role': 'reference_site', 'chrom': event['chrom'],
                'start': event['start'], 'end': max(event['start'] + 1, event['end'])})
            if 'source_start' in event:
                row['reference_intervals'].append({'role': 'copy_source', 'chrom': event['chrom'],
                    'start': event['source_start'], 'end': event['source_end']})
            events.append(row)
    for handle in handles.values():
        handle.close()
    oldhash = defaultdict(list)
    for sha, sequence in sequences.items():
        oldhash[hashlib.sha1(sequence.encode()).hexdigest()[:16]].append(sha)
    verified = set()
    for path in sorted((MEI / 'ins_arm').glob('shard_*/query.fa')):
        for key, sequence in read_fasta_selected(path, oldhash):
            sha = digest(sequence)
            if sha not in oldhash[key]:
                raise ValueError('Truncated sequence hash collision')
            verified.add(sha)
    rm_hits = defaultdict(list)
    for path in sorted((MEI / 'ins_arm').glob('shard_*/query.fa.out')):
        for hit in parse_rm(path):
            for sha in oldhash.get(hit['query'], []):
                if sha in verified:
                    rm_hits[sha].append(hit)
    missing_literal = sorted(literal_hashes - verified)
    aliases = tool_aliases(sequences)
    for name, names in [('actual_changed_sequences.fa', sequences), ('missing_literal_INS.fa', missing_literal)]:
        with (work / name).open('w') as handle:
            for sha in names:
                handle.write(f'>{sha[:32]}\n{sequences[sha]}\n')
    metadata = {'manifest_sha256': hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                'sequence_lengths': {k: len(v) for k, v in sequences.items()},
                'tool_id_to_sha256': aliases,
                'literal_INS_sequences': sorted(literal_hashes), 'prior_verified': sorted(verified),
                'missing_literal_INS': missing_literal, 'prior_rm_hits': dict(rm_hits), 'events': events,
                'sources': {'prior_rm': str(MEI / 'ins_arm/shard_*/query.fa.out'),
                            'reference_rm': str(MEI / 'chm13_hor_rm/shards/shard_*/q.fa.out'),
                            'monomers': str(ATLAS / 'monomers.fa'),
                            'provenance': str(MEI / 'step0/MANIFEST.txt')},
                'rules': {'alpha_predominance_fraction': .8, 'SD_identity_floor': 80,
                          'TE_support': 'Any RepeatMasker LINE/SINE/LTR/DNA/RC/Retroposon homology; not proof of transposition',
                          'negative_claim': 'No TE annotation is not a validated absence; SD-only evidence does not test TE content',
                          'reference_context': 'Inherited CHM13 composition, separate from actual SNV-mutated sequence',
                          'reference_scope': '41 active-HOR intervals, not a whole-genome annotation',
                          'mixed': 'Group coverages are interval unions within class and may overlap across classes'}}
    atomic_json(work / 'prepared.json', metadata)
    print(json.dumps({'unique_sequences': len(sequences), 'total_bp': sum(map(len, sequences.values())),
                      'literal_INS_sequences': len(literal_hashes), 'prior_verified': len(verified),
                      'missing_literal_INS': len(missing_literal),
                      'missing_literal_INS_bp': sum(len(sequences[k]) for k in missing_literal)}, indent=2))


def run(out):
    import os
    work = out / 'repeat_annotation'
    meta = json.loads((work / 'prepared.json').read_text())
    commands = []
    if meta['missing_literal_INS']:
        rm_dir = work / 'rm_missing'; rm_dir.mkdir(exist_ok=True)
        command = [str(RM_ENV / 'bin/RepeatMasker'), '-s', '-species', 'human', '-pa', '1',
                   '-dir', str(rm_dir), str(work / 'missing_literal_INS.fa')]
        commands.append(command)
        env = dict(os.environ, PATH=str(RM_ENV / 'bin') + ':' + os.environ.get('PATH', ''), OMP_NUM_THREADS='4')
        with (work / 'rm_missing.log').open('w') as log:
            subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    command = [str(SD_ENV / 'bin/stringdecomposer'), str(work / 'actual_changed_sequences.fa'),
               str(ATLAS / 'monomers.fa'), '-t', '4', '-o', str(work / 'sd'), '--out-file', 'decomp']
    commands.append(command)
    env = dict(os.environ, PATH=str(SD_ENV / 'bin') + ':' + os.environ.get('PATH', ''), OMP_NUM_THREADS='4')
    with (work / 'sd.log').open('w') as log:
        subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    atomic_json(work / 'commands.json', {'commands': commands})
    collect(out)


def reference_hits():
    hits, intervals = defaultdict(list), defaultdict(set)
    for path in sorted((MEI / 'chm13_hor_rm/shards').glob('shard_*/q.fa.out')):
        for hit in parse_rm(path):
            chrom, interval = hit['query'].split(':')
            start, end = map(int, interval.split('-'))
            intervals[chrom].add((start, end))
            hit = dict(hit, start=start + hit['start'], end=start + hit['end'])
            hits[chrom].append(hit)
    return hits, intervals


def collect(out):
    work = out / 'repeat_annotation'
    meta = json.loads((work / 'prepared.json').read_text())
    if meta['manifest_sha256'] != hashlib.sha256((out / 'manifest.json').read_bytes()).hexdigest():
        raise ValueError('Manifest changed after annotation preparation')
    rm_hits = meta['prior_rm_hits']
    newly_tested = set(meta['missing_literal_INS'])
    if newly_tested:
        raw = work / 'rm_missing/missing_literal_INS.fa.out'
        if not raw.exists():
            raise ValueError('Missing RepeatMasker output; do not infer zero hits from a failed run')
        for hit in parse_rm(raw):
            sha = meta['tool_id_to_sha256'][hit['query']]
            if sha not in newly_tested:
                raise ValueError('Unexpected sequence in RepeatMasker output')
            rm_hits.setdefault(sha, []).append(dict(hit, query=sha))
    tested = set(meta['prior_verified']) | newly_tested
    sd = parse_sd(work / 'sd/decomp.tsv', meta['sequence_lengths'], meta['tool_id_to_sha256'])
    component_results = {}
    for sha, length in meta['sequence_lengths'].items():
        rm = coverage(rm_hits.get(sha, []), length) if sha in tested else None
        component_results[sha] = {'length': length, 'rm': rm, 'sd': sd.get(sha),
                                 'rm_source': 'exact_prior_sequence' if sha in meta['prior_verified'] else
                                              'new_literal_INS_mask' if sha in newly_tested else 'not_directly_masked',
                                 'class': classify_component(length, rm, sd.get(sha))}
    ref_hits, ref_intervals = reference_hits()
    rows = []
    for event in meta['events']:
        components = [component_results[c['sha256']] for c in event['components']]
        classes = {c['class'] for c in components}
        if classes & {'TE_supported', 'mixed_TE_satellite'}:
            label = 'TE_or_mixed_supported'
        elif classes and classes <= {'alpha_satellite_RM', 'alpha_like_SD_only'}:
            label = 'alpha_satellite_supported'
        else:
            label = 'other_or_unresolved'
        contexts = []
        for interval in event['reference_intervals']:
            start, end, chrom = interval['start'], interval['end'], interval['chrom']
            hits = [dict(h, start=max(h['start'], start)-start, end=min(h['end'], end)-start)
                    for h in ref_hits.get(chrom, []) if h['start'] < end and h['end'] > start]
            context = coverage(hits, end-start)
            context['assayed_bp'] = union_bp([(max(s,start)-start, min(e,end)-start)
                for s,e in ref_intervals.get(chrom,[]) if s < end and e > start])
            contexts.append(dict(interval, composition=context))
        event['reference_context'] = contexts
        event['composition_class'] = label
        rows.append({'case': event['case'], 'id': event['id'], 'kind': event['kind'],
                     'composition_class': label, 'component_classes': ';'.join(sorted(classes)),
                     'n_components': len(components), 'n_direct_RM_components': sum(c['rm'] is not None for c in components),
                     'n_SD_components': sum(c['sd'] is not None for c in components),
                     'component_sha256': ';'.join(c['sha256'] for c in event['components'])})
    result = dict(meta, sequence_annotations=component_results,
                  summary={'n_events': len(rows), 'n_unique_sequences': len(component_results),
                           'n_direct_RM_sequences': len(tested), 'n_SD_sequences': len(sd),
                           'event_classes': dict(Counter(r['composition_class'] for r in rows)),
                           'unique_sequence_classes': dict(Counter(c['class'] for c in component_results.values()))})
    sources = sorted((MEI / 'ins_arm').glob('shard_*/query.fa.out'))
    sources += sorted((MEI / 'chm13_hor_rm/shards').glob('shard_*/q.fa.out'))
    sources += [ATLAS / 'monomers.fa', ATLAS / 'placement.json', MEI / 'step0/MANIFEST.txt',
                Path(__file__).resolve(), work / 'prepared.json', work / 'sd/decomp.tsv']
    if newly_tested:
        sources.append(work / 'rm_missing/missing_literal_INS.fa.out')
    result['source_sha256'] = {str(p): file_digest(p) for p in sources}
    result['instruments'] = {'repeatmasker': '4.1.9', 'dfam': '3.9, partitions 0 and 7',
                             'rmblast': '2.17.1+', 'stringdecomposer': '1.1.2',
                             'instrument_provenance': str(MEI / 'step0/MANIFEST.txt'),
                             'SD_method': str(UPSTREAM / 'discovery/2026-09-23_monomer_atlas/code/06_persv/table/sd_alleles.py')}
    atomic_json(out / 'repeat_annotation.json', result)
    receipt = work / 'submission_receipt.json'
    if receipt.exists():
        record = json.loads(receipt.read_text())
        record.update(status='completed', annotation_summary=result['summary'],
                      annotation_sha256=file_digest(out / 'repeat_annotation.json'))
        atomic_json(receipt, record)
    temp = out / 'repeat_annotations.tsv.tmp'
    with temp.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    temp.replace(out / 'repeat_annotations.tsv')
    print(json.dumps(result['summary'], indent=2))


def report_repeat_section(path, assets, events, tex, table, figure):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    path, assets = Path(path), Path(assets)
    data = json.loads(path.with_name('repeat_annotation.json').read_text())
    if data['manifest_sha256'] != file_digest(path.with_name('manifest.json')):
        raise ValueError('Repeat annotations are stale relative to the report manifest')
    with path.open() as handle:
        ann = {(r['case'], r['id']): r for r in csv.DictReader(handle, delimiter='\t')}
    expected = {(c['id'], e['id']) for c in json.loads(path.with_name('manifest.json').read_text())['cases'] for e in c['truth']}
    if set(ann) != expected or len(ann) != data['summary']['n_events']:
        raise ValueError('Repeat annotation event identities do not match the manifest')
    summary, seqs = data['summary'], data['sequence_annotations']
    n_literal = len(data['literal_INS_sequences'])
    literal_te = sum(seqs[s]['rm']['group_bp'].get('transposable_element', 0) > 0
                     for s in data['literal_INS_sequences'] if seqs[s]['rm'] is not None)
    text = [r'\subsection{Changed-allele repeat composition}', tex(
        f"The {summary['n_events']} injected event instances contain {summary['n_unique_sequences']} distinct changed sequences. "
        f"{summary['n_direct_RM_sequences']} sequences have direct RepeatMasker evidence and {summary['n_SD_sequences']} have StringDecomposer output. "
        "Annotations describe the actual inserted and removed sequences; CHM13 site/source annotations are retained separately and are not substituted for inserted-allele composition. "
        "Repeat content does not establish the mutation mechanism or prove that complete monomers were gained or lost."),
        tex(f"All {n_literal} distinct panel insertion payloads were tested with RepeatMasker; {literal_te} have annotated TE homology. "
            "TE-rich insertion sensitivity remains unmeasured in the absence of such payloads. A negative RepeatMasker annotation is instrument-dependent and is not proof of biological absence."),
        tex("RepeatMasker 4.1.9 with Dfam 3.9 and RMBLAST is reused only after exact sequence verification; missing panel insertion alleles are masked with the same instrument. "
            "The alpha-satellite class requires at least 80% ALR/Alpha coverage, accepting either slash spelling in the annotation. "
            "For sequences without direct masking, alpha-like support requires at least 80% coverage by reliable StringDecomposer blocks with identity at least 80% to the chr1 monomer library. "
            "These descriptive annotation thresholds do not alter SV scoring. SD-only support does not test TE absence. TE homology and satellite homology may coexist; interval unions avoid double-counting within a class.")]
    unknown = [a for a in seqs.values() if a['class'] == 'unresolved']
    direct_alpha = [a for a in seqs.values() if a['class'] == 'alpha_satellite_RM']
    sd_below = sum(a['sd'] is None or a['sd']['high_identity_bp'] / a['length'] < .8
                   for a in direct_alpha)
    if unknown:
        text.append(tex(
            f"The {len(unknown)} unclassified distinct sequences were {min(a['length'] for a in unknown)}--{max(a['length'] for a in unknown)} bp long. "
            f"All received decomposition output, but failed the conservative monomer-coverage/identity rule. In this selected set, {sd_below}/{len(direct_alpha)} "
            "directly RepeatMasker-supported alpha-satellite sequences also failed the SD-only rule. Therefore, 'other or unresolved' means repeat-annotation uncertainty, "
            "not evidence for a distinct biological repeat class or absence of alpha-satellite. It is separate from unresolved SV recovery. "
            "Recovery differences between annotation groups should not be interpreted as effects of repeat class; the groups also differ in sequence length and available annotation evidence."))
    counts = Counter(c['class'] for c in seqs.values())
    text.append(table(['Evidence class', 'Distinct sequences'], sorted(counts.items()),
        'Unique actual changed sequences, not independent natural variants. The same sequence may recur in many queries.', 'tab:repeat_evidence', 'lr'))
    if literal_te == 0:
        text.append(tex("No tested panel insertion payload had direct TE-positive homology. Consequently, this experiment cannot compare TE-insertion sensitivity with satellite-insertion sensitivity. A TE-rich allele challenge would be required for that comparison; good performance on alpha-satellite payloads alone does not establish TE recovery."))
    stats = defaultdict(lambda: [0, 0])
    for event in events:
        if event.get('scope') != 'raw_full':
            continue
        row = ann.get((event['case'], event['id']))
        if row is None:
            raise ValueError('Scored event missing repeat annotation')
        key = (row['composition_class'], event['kind'])
        stats[key][0] += 1
        stats[key][1] += event.get('recovered') in (True, 'True')
    text.append(table(['Composition', 'SV type', 'Strict recovered / injected'],
        [(c.replace('_', ' '), k, f'{v[1]}/{v[0]}') for (c,k),v in sorted(stats.items())],
        'Truth-composition-stratified strict recall using unfiltered caller output from the chr1 query. Unexplained calls were not composition-annotated, so class-specific precision is not estimated.',
        'tab:repeat_recovery', 'p{3.1in}lr'))
    groups = sorted(stats)
    if groups:
        fig, ax = plt.subplots(figsize=(8, max(3, .38*len(groups))))
        fractions = [stats[g][1]/stats[g][0] for g in groups]
        ax.barh(range(len(groups)), fractions, color='#247b86')
        ax.set_yticks(range(len(groups)), [f'{k}: {c.replace("_", " ")}' for c,k in groups], fontsize=8)
        for y,g in enumerate(groups):
            ax.text(min(.98, fractions[y] + .02), y, f'{stats[g][1]}/{stats[g][0]}', va='center', fontsize=8)
        ax.set(xlim=(0,1.15), xlabel='Strict event recovery fraction')
        ax.spines[['top','right']].set_visible(False)
        fig.savefig(assets / 'repeat_composition_recovery.pdf', bbox_inches='tight')
        fig.savefig(assets / 'repeat_composition_recovery.png', dpi=160, bbox_inches='tight'); plt.close(fig)
        text.append(figure('repeat_composition_recovery.pdf',
            'Sequence-composition evidence and event type using unfiltered caller output from the chr1 query. Fractions pool paired mutation recipes across divergence and contiguity conditions and are descriptive, not independent-sample confidence intervals.',
            'fig:repeat_recovery'))
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('results/complete_truth_v1'))
    parser.add_argument('--stage', choices=['prepare', 'run', 'collect'], required=True)
    args = parser.parse_args()
    globals()[args.stage](args.output_dir.resolve())


if __name__ == '__main__':
    main()
