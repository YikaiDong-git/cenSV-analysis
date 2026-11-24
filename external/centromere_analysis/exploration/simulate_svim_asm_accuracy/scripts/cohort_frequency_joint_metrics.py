"""Safety primitives for joint cohort analysis across source chromosomes."""
from collections import Counter

from cohort_frequency_metrics import frequency_bin, row_class


def joint_call_id(source_chrom, haplotype_index, original_id):
    if not source_chrom or '__' in source_chrom or not original_id or ';' in original_id:
        raise ValueError('Invalid source chromosome or call ID')
    if not 0 <= haplotype_index < 256:
        raise ValueError('Invalid haplotype index')
    return f'{source_chrom}__H{haplotype_index:03d}__{original_id}'


def joint_event_id(truth_chrom, event_id):
    if not truth_chrom or '__' in truth_chrom or not event_id:
        raise ValueError('Invalid truth chromosome or event ID')
    return f'{truth_chrom}__{event_id}'


def namespace_truth(source_chrom, truth):
    result = []
    for event in truth:
        if event['chrom'] != source_chrom:
            raise ValueError('Truth event is not on its source chromosome')
        entry = dict(event)
        entry['original_event_id'] = event['id']
        entry['id'] = joint_event_id(source_chrom, event['id'])
        result.append(entry)
    if len({e['id'] for e in result}) != len(result):
        raise ValueError('Duplicated truth event ID')
    return result


def joint_filter_and_deduplicate(rows, n):
    best = {}
    for index, row in enumerate(rows):
        if row['row_index'] != index:
            raise ValueError('Row indices must be contiguous and ordered')
        ac = row['ac']
        if not 0 <= ac <= n:
            raise ValueError('Invalid carrier count')
        row['af'] = ac / n
        row['maf'] = min(row['af'], 1 - row['af'])
        row['frequency_bin'] = frequency_bin(row['maf'])
        row['nonreference_major'] = row['af'] > .5
        row['passes_ac_af'] = ac >= 2 and ac <= int(n * .99)
        row['dedup_retained'] = False
        if row['passes_ac_af']:
            key = (row['chrom'], row['pos'], row['end'], row['svtype'], row['svlen'])
            if key not in best or ac > rows[best[key]]['ac']:
                best[key] = index
    for index in best.values():
        rows[index]['dedup_retained'] = True
    return rows


def classify_joint_carriers(carrier_members):
    if not carrier_members:
        raise ValueError('A merged row must have a carrier')
    carrier_statuses = []
    ambiguous = []
    for haplotype, statuses in carrier_members.items():
        if not statuses or set(statuses) - {'tp', 'fp', 'unresolved'}:
            raise ValueError('Missing or invalid source-call status')
        if len(statuses) != 1:
            ambiguous.append(haplotype)
            carrier_statuses.append('unresolved')
        else:
            carrier_statuses.append(statuses[0])
    counts = Counter(carrier_statuses)
    return dict(ac=len(carrier_members), classification=row_class(carrier_statuses),
                n_tp=counts['tp'], n_fp=counts['fp'], n_unresolved=counts['unresolved'],
                ambiguous_multi_source_carriers=sorted(ambiguous))


def load_joint_sources(joint_root, n):
    import json
    from pathlib import Path
    from cohort_frequency_metrics import sample_name

    root = Path(joint_root)
    sources, score_events, cases = {}, {}, {}
    for h in range(n):
        folder = root / 'haps' / f'H{h:03d}'
        if not (folder / 'complete.json').exists():
            raise FileNotFoundError(f'Incomplete joint sample {folder}')
        case = json.loads((folder / 'truth.json').read_text())
        score = json.loads((folder / 'scoring/production_active.json').read_text())
        calls = json.loads((folder / 'scoring/parsed_calls_production_active.json').read_text())
        if case['haplotype_index'] != h or case['sample_name'] != sample_name(h):
            raise ValueError(f'Joint sample identity mismatch: {folder}')
        if len(score['calls']) != len(calls):
            raise ValueError(f'Joint call count mismatch: {folder}')
        events = {e['id']: e for e in score['events']}
        if len(events) != len(case['truth']) or set(events) != {e['id'] for e in case['truth']}:
            raise ValueError(f'Joint truth IDs do not match score: {folder}')
        for call, parsed in zip(score['calls'], calls):
            if call['id'] != parsed['id'] or call['status'] not in ('tp', 'fp', 'unresolved'):
                raise ValueError(f'Joint call score mismatch: {folder}')
            identity = f'H{h:03d}__{call["id"]}'
            if identity in sources:
                raise ValueError(f'Duplicate joint source ID: {identity}')
            sources[identity] = dict(h=h, sample=sample_name(h), status=call['status'],
                matched_truth_ids=call.get('matched_truth_ids', []), raw_id=call['id'], full=parsed)
        score_events[h] = events
        cases[h] = case
    return sources, score_events, cases


def validate_joint_vcfs(joint_root, n, sources):
    import pysam
    from pathlib import Path
    from cohort_frequency_metrics import presence, sample_name

    paths = []
    for h in range(n):
        path = Path(joint_root) / 'haps' / f'H{h:03d}' / 'merge_input.vcf.gz'
        paths.append(path)
        with pysam.VariantFile(str(path)) as vcf:
            if list(vcf.header.samples) != [sample_name(h)]:
                raise ValueError(f'Wrong joint VCF sample: {path}')
            found = set()
            for rec in vcf:
                identity = rec.id
                if not identity or identity not in sources or sources[identity]['h'] != h or identity in found:
                    raise ValueError(f'Invalid joint VCF source ID: {identity}')
                gt = rec.samples[sample_name(h)].get('GT')
                if not presence(gt):
                    raise ValueError(f'Joint source call is not present: {identity}')
                alleles = {x for x in gt if x > 0}
                if len(alleles) != 1 or not rec.alts or max(alleles) > len(rec.alts):
                    raise ValueError(f'Ambiguous joint source allele: {identity}')
                sources[identity].update(source_alt=rec.alts[next(iter(alleles))-1],
                                         source_ref=rec.ref, source_start=rec.start)
                found.add(identity)
            expected = {k for k, value in sources.items() if value['h'] == h}
            if found != expected:
                raise ValueError(f'Joint VCF does not match saved call rows: {path}')
    return paths


def collect_all_reported_chromosomes(merged_vcf, n, sources):
    from cohort_frequency_metrics import _collect_merged

    reported = sorted({s['full'].get('reported_chrom', s['full'].get('chrom'))
                       for s in sources.values()})
    if any(not chrom for chrom in reported):
        raise ValueError('Joint source is missing its reported chromosome')
    rows, memberships = [], []
    for chrom in reported:
        local_rows, local_members = _collect_merged(merged_vcf, chrom, n, sources)
        offset = len(rows)
        for row in local_rows:
            row['row_index'] += offset
        for member in local_members:
            member['row_index'] += offset
        rows.extend(local_rows)
        memberships.extend(local_members)
    if {member['source_id'] for member in memberships} != set(sources):
        raise ValueError('Joint merged output does not account for every source record')
    return rows, memberships, reported


def joint_truth_recovery(manifest, score_events, memberships, rows, filtered_scores):
    from cohort_frequency_metrics import _event_recovery

    n = int(manifest['config']['n_haplotypes'])
    result = []
    for chrom in manifest['config']['chromosomes']:
        events = []
        for entry in manifest['plans'][chrom]['events']:
            namespaced = dict(entry)
            namespaced['id'] = joint_event_id(chrom, entry['id'])
            events.append(namespaced)
        part = _event_recovery(chrom, {'events': events}, n, score_events,
                               memberships, rows, filtered_scores)
        for row in part:
            row['original_event_id'] = row['event_id'].split('__', 1)[1]
        result.extend(part)
    return result


def aggregate_joint(manifest, rows, memberships, event_carriers, output):
    import json
    import shutil
    import tempfile
    from pathlib import Path
    from cohort_frequency_metrics import aggregate

    output = Path(output)
    if (output / 'merged').exists():
        raise FileExistsError(f'Joint aggregate already exists: {output / "merged"}')
    output.mkdir(parents=True, exist_ok=True)
    chroms = manifest['config']['chromosomes']
    extras = [row for row in rows if row['chrom'] not in chroms]
    with tempfile.TemporaryDirectory(prefix='joint_cohort_aggregate_', dir=output) as temporary:
        temp = Path(temporary)
        (temp / 'manifest.json').write_text(json.dumps(manifest) + '\n')
        for chrom in chroms:
            folder = temp / 'analysis' / chrom
            folder.mkdir(parents=True)
            grouped_rows = [r for r in rows if r['chrom'] == chrom]
            if chrom == chroms[0]:
                grouped_rows.extend(extras)
            entry = dict(chrom=chrom, rows=grouped_rows,
                         event_carriers=[r for r in event_carriers if r['chrom'] == chrom])
            (folder / 'analysis.json').write_text(json.dumps(entry) + '\n')
        summary = aggregate(temp / 'manifest.json')
        dest = output / 'merged'
        if dest.exists():
            raise FileExistsError(f'Joint aggregate already exists: {dest}')
        shutil.copytree(temp / 'analysis' / 'merged', dest)
    reported_chromosomes = sorted({row['chrom'] for row in rows})
    for chrom in sorted(set(chroms) | set(reported_chromosomes)):
        folder = output / chrom
        folder.mkdir(exist_ok=True)
        indices = {row['row_index'] for row in rows if row['chrom'] == chrom}
        entry = dict(chrom=chrom, rows=[r for r in rows if r['chrom'] == chrom],
                     memberships=[m for m in memberships if m['row_index'] in indices],
                     event_carriers=[r for r in event_carriers if r['chrom'] == chrom])
        (folder / 'analysis.json').write_text(json.dumps(entry, indent=2) + '\n')
    summary['reported_chromosomes'] = reported_chromosomes
    summary['truth_chromosomes'] = list(chroms)
    summary['joint_cohort'] = True
    (output / 'merged' / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


def analyze_joint(manifest_path, rescore_fn, *, bcftools=None):
    import json
    from pathlib import Path
    import subprocess
    from cohort_frequency_metrics import (BCFTOOLS, _apply_filtered_scores,
                                           _write_tsv, validate_manifest)

    manifest_path = Path(manifest_path)
    root = manifest_path.parent
    manifest = json.loads(manifest_path.read_text())
    validate_manifest(manifest)
    n = int(manifest['config']['n_haplotypes'])
    joint_root = root / 'joint_cohort'
    sources, score_events, _ = load_joint_sources(joint_root, n)
    paths = validate_joint_vcfs(joint_root, n, sources)
    output = root / 'joint_analysis'
    output.mkdir(parents=True, exist_ok=True)
    listing = output / 'merge_inputs.txt'
    listing.write_text(''.join(str(path.resolve()) + '\n' for path in paths))
    merged = output / 'merged.vcf.gz'
    command = [str(bcftools or BCFTOOLS), 'merge', '-l', str(listing), '-Oz', '-o', str(merged)]
    subprocess.run(command, check=True)
    subprocess.run([str(bcftools or BCFTOOLS), 'index', '-f', '-t', str(merged)], check=True)
    rows, members, reported = collect_all_reported_chromosomes(merged, n, sources)
    rescored = rescore_fn(root, joint_root, manifest, rows, members, sources)
    rows, members = _apply_filtered_scores(rows, members, rescored)
    event_carriers = joint_truth_recovery(manifest, score_events, members, rows, rescored)
    _write_tsv(output / 'merged_rows.tsv', rows, ('row_index', 'chrom'))
    _write_tsv(output / 'source_membership.tsv', members, ('chrom', 'row_index', 'source_id'))
    _write_tsv(output / 'event_carriers.tsv', event_carriers, ('chrom', 'event_id', 'haplotype_index'))
    global_result = dict(rows=rows, memberships=members, event_carriers=event_carriers,
                         reported_chromosomes=reported, merged_vcf=str(merged), command=command,
                         n_haplotypes=n)
    (output / 'analysis.json').write_text(json.dumps(global_result, indent=2, allow_nan=False) + '\n')
    summary = aggregate_joint(manifest, rows, members, event_carriers, output)
    return dict(reported_chromosomes=reported, n_rows=len(rows), n_memberships=len(members),
                summary=summary, merged_vcf=str(merged), command=command)
