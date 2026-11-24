#!/usr/bin/env python3
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
ALLOWED_NODES = {'node03', 'node04', 'node08', 'node12'}
CHROMS = {'chr1', 'chr7', 'chr13', 'chr19'}
THRESHOLDS = (0, 1, 5, 10)
BINS = ('<1%', '1-<5%', '5-<10%', '>=10%')


class AuditFailure(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise AuditFailure(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_json(path, evidence=None):
    path = Path(path)
    data = path.read_bytes()
    if evidence is not None:
        evidence[str(path)] = hashlib.sha256(data).hexdigest()
    return json.loads(data)


def same(actual, expected, label):
    if isinstance(expected, float):
        require(isinstance(actual, (int, float)) and math.isfinite(actual)
                and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12), label)
    else:
        require(actual == expected, label)


def fields(actual, expected, label):
    for key, value in expected.items():
        require(key in actual, f'{label}: missing {key}')
        same(actual[key], value, f'{label}: {key}: {actual[key]!r} != {value!r}')


def ratio(a, b):
    return a / b if b else None


def classify(statuses):
    c = Counter(statuses)
    require(c and not (set(c) - {'tp', 'fp', 'unresolved'}), 'Invalid carrier statuses')
    kind = ('any_unresolved' if c['unresolved'] else 'mixed_tp_fp' if c['tp'] and c['fp']
            else 'all_tp' if c['tp'] else 'all_fp')
    return dict(classification=kind, n_tp=c['tp'], n_fp=c['fp'], n_unresolved=c['unresolved'])


def frequency_bin(value):
    return BINS[0] if value < .01 else BINS[1] if value < .05 else BINS[2] if value < .1 else BINS[3]


def mixed(row):
    t, f, u = row['n_tp'], row['n_fp'], row['n_unresolved']
    return t > 0 and f > 0, t + f + u >= 2 and t + u >= 1 and f + u >= 1


def audit_score(score, truth_ids, call_ids, params):
    events, calls, metrics = score['events'], score['calls'], score['metrics']
    require(len(events) == len(truth_ids) and {e['id'] for e in events} == set(truth_ids),
            'Score lost or duplicated truth events')
    require(len(calls) == len(call_ids) and {c['id'] for c in calls} == set(call_ids),
            'Score lost or duplicated call records')
    c = Counter(row['status'] for row in calls)
    require(not (set(c) - {'tp', 'fp', 'unresolved'}), 'Unknown call status')
    t = sum(e['recovered'] is True for e in events)
    u = sum(e['status'] == 'unresolved' for e in events)
    require(all(not e['recovered'] for e in events if e['status'] == 'unresolved'),
            'Recovered event also unresolved')
    fields(metrics, dict(n_truth=len(events), n_calls=len(calls), tp_events=t,
        unresolved_events=u, fn_events=len(events)-t-u, fn_events_strict=len(events)-t,
        tp_calls=c['tp'], fp_calls=c['fp'], unresolved_calls=c['unresolved'],
        sensitivity_lower=ratio(t, len(events)), sensitivity_upper=ratio(t+u, len(events)),
        precision_lower=ratio(c['tp'], len(calls)),
        precision_upper=ratio(c['tp']+c['unresolved'], len(calls))), 'Score metrics')
    fields(score['conventions'], params, 'Canonical scoring parameters')
    for event in events:
        if event['recovered']:
            require(event['matched_call_ids'] and set(event['matched_call_ids']) <= set(call_ids),
                    'Recovered event lacks surviving supporting records')
    return ({r['id']: r['status'] for r in calls}, {r['id']: r for r in events})


def _option(command, name):
    for i, token in enumerate(command):
        if token.startswith(name + '='):
            return token.split('=', 1)[1]
        if token == name and i + 1 < len(command):
            return command[i + 1]
    return None


def audit_resources(ledger, config, observed_nodes):
    allowed = set(ledger['node_allowlist'])
    require(allowed and allowed <= ALLOWED_NODES and len(allowed) <= 4, 'Invalid resource node allowlist')
    require(set(config['node_allowlist']) == allowed, 'Config/ledger node disagreement')
    require(observed_nodes <= allowed, 'Completed case ran outside allowed nodes')
    require(0 < ledger['cpus_per_job'] <= 4 and 0 < config['cpus_per_job'] <= 4,
            'CPU allocation exceeds four')
    count = 0
    for entry in ledger['submissions']:
        tasks, command = entry['tasks'], entry['command']
        require(isinstance(tasks, int) and tasks > 0, 'Invalid allocation count')
        count += tasks
        if 'indices' in entry:
            require(len(entry['indices']) == len(set(entry['indices'])) == tasks, 'Array task count mismatch')
        if command[0] == 'sbatch':
            require(_option(command, '--cpus-per-task') is not None
                    and 0 < int(_option(command, '--cpus-per-task')) <= 4, 'Unbounded job CPU request')
            require(_option(command, '--nodes') == '1' and _option(command, '--ntasks') == '1',
                    'Unexpected distributed allocation')
            excluded = set((_option(command, '--exclude') or '').split(','))
            require(not (allowed & excluded), 'Submitted job excludes declared allowed node')
        else:
            require(command[:2] == ['scontrol', 'requeue'], 'Unknown allocation ledger operation')
    fields(ledger, dict(submitted_new_tasks=count,
                        total_submitted_tasks=ledger['previous_tasks'] + count), 'Resource ledger')
    require(0 <= ledger['previous_tasks'] and ledger['total_submitted_tasks'] <= ledger['user_ceiling'] <= 150,
            'Resource allocation ceiling exceeded')
    return dict(total_allocations=ledger['total_submitted_tasks'], observed_nodes=sorted(observed_nodes),
                allowed_nodes=sorted(allowed), cpus_per_allocation=ledger['cpus_per_job'],
                scope='Submission ledger and completed-case node receipts; no new scheduler queries.')


def identity_catalog(root, manifest, evidence):
    archives = list((root / 'code_snapshots').rglob('*')) + list(root.glob('*.json'))
    archive_hashes = {sha(p) for p in archives if p.is_file()}
    candidates = [root / 'execution_code_identity.json', root / 'analysis_code_identity.json']
    candidates += list((root / 'code_snapshots').rglob('*identity.json'))
    candidates += list((root / 'code_snapshots').rglob('previous_identity.json'))
    catalog = {}
    for path in set(candidates):
        mapping = read_json(path, evidence)
        require(mapping and all(isinstance(v, str) and re.fullmatch('[0-9a-f]{64}', v)
                                for v in mapping.values()), f'Invalid identity map: {path}')
        for name, expected in mapping.items():
            source = Path(name)
            available = expected in archive_hashes or (source.is_file() and sha(source) == expected)
            require(available, f'Unpreserved code epoch source: {path}: {name}')
            if source.name == 'manifest.json':
                copies = [p for p in root.glob('*.json') if sha(p) == expected]
                if source.is_file() and sha(source) == expected:
                    copies.append(source)
                require(copies, f'Epoch manifest not archived: {path}')
                prior = read_json(copies[0])
                for key in ('plans', 'config', 'reference', 'active_bed', 'index', 'payloads',
                            'caller_parameters', 'alignment_arguments'):
                    same(prior[key], manifest[key], f'Biological/calling design changed across epochs: {key}')
        catalog[sha(path)] = str(path)
    return catalog


def table_matches(path, expected, evidence):
    evidence[str(path)] = sha(path)
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        actual = list(reader)
    require(len(actual) == len(expected), f'Table row count differs: {path}')
    def cell(value):
        if isinstance(value, (dict, list, tuple)):
            return json.dumps(value, separators=(',', ':'), sort_keys=True)
        return '' if value is None else str(value)
    for i, (row, wanted) in enumerate(zip(actual, expected)):
        fields(row, {k: cell(v) for k, v in wanted.items()}, f'{path}: row {i}')


def group_metrics(rows):
    c = Counter(r['classification'] for r in rows)
    t, f, u = (sum(r[key] for r in rows) for key in ('n_tp', 'n_fp', 'n_unresolved'))
    n = len(rows)
    return dict(n_rows=n, all_tp=c['all_tp'], all_fp=c['all_fp'], mixed_tp_fp=c['mixed_tp_fp'],
        any_unresolved=c['any_unresolved'], proven_mixed=sum(mixed(r)[0] for r in rows),
        possible_mixed=sum(mixed(r)[1] for r in rows),
        all_fp_lower=ratio(c['all_fp'], n), all_fp_upper=ratio(sum(r['n_tp'] == 0 for r in rows), n),
        any_fp_lower=ratio(sum(r['n_fp'] > 0 for r in rows), n),
        any_fp_upper=ratio(sum(r['n_fp'] > 0 or r['n_unresolved'] > 0 for r in rows), n),
        carrier_tp=t, carrier_fp=f, carrier_unresolved=u,
        carrier_fp_lower=ratio(f, t+f+u), carrier_fp_upper=ratio(f+u, t+f+u))


def audit_chromosome(root, chrom, manifest, epochs, evidence):
    n = manifest['config']['n_haplotypes']
    plans = {e['id']: e for e in manifest['plans'][chrom]['events']}
    founder_dir = root / 'chromosomes' / chrom / 'founders'
    fields(read_json(founder_dir / 'complete.json', evidence), dict(completed=True, founders=8), 'Founder completion')
    for i in range(8):
        founder = read_json(founder_dir / f'F{i}.json', evidence)
        require(founder.get('ledger_replay_verified') is True, 'Founder SNP replay proof missing')
        require(re.fullmatch('[0-9a-f]{64}', founder['sequence_sha256']), 'Invalid founder sequence identity')
    analysis = read_json(root / 'analysis' / chrom / 'analysis.json', evidence)
    fields(analysis, dict(chrom=chrom, n_haplotypes=n), 'Chromosome analysis identity')
    rows, members = analysis['rows'], analysis['memberships']
    require([r['row_index'] for r in rows] == list(range(len(rows))), 'Noncontiguous row indices')
    members_by_row, member_ids, retained = defaultdict(list), set(), defaultdict(set)
    for member in members:
        identity = member['source_id']
        h = member['haplotype_index']
        require(identity not in member_ids and identity.startswith(f'H{h:03d}__') and 0 <= h < n,
                'Duplicated or misassigned source membership')
        require(0 <= member['row_index'] < len(rows), 'Invalid membership row index')
        member_ids.add(identity)
        members_by_row[member['row_index']].append(member)
        row = rows[member['row_index']]
        for p in THRESHOLDS:
            if row['dedup_retained'] and row['maf'] >= p / 100:
                retained[h, p].add(identity.split('__', 1)[1])
    nodes, templates, source_ids, used_epochs = set(), set(), set(), Counter()
    call_status, event_status = {}, {}
    for h in range(n):
        folder = root / 'chromosomes' / chrom / 'haps' / f'H{h:03d}'
        done, truth, oracle = (read_json(folder / name, evidence)
                               for name in ('complete.json', 'truth.json', 'oracle.json'))
        require(done['completed'] is True and oracle['passed'] is True, f'Incomplete case/oracle: {folder}')
        required_checks = ('all_truth_recovered', 'no_false_or_unresolved_calls', 'no_unresolved_events',
                           'whole_sequence_reconstructed', 'no_false_negative_windows')
        require(all(oracle['checks'].get(k) is True for k in required_checks), f'Failed oracle checks: {folder}')
        require(truth.get('ledger_replay_verified') is True, f'Missing SNP replay proof: {folder}')
        fields(truth, dict(id=f'{chrom}_H{h:03d}', chrom=chrom, haplotype_index=h,
                          sample_name=f'SIM{h//2:03d}_hap{h%2+1}'), 'Truth identity')
        fields(done, dict(case_id=truth['id'], query_sha256=truth['query_sha256'],
                          template_sha256=truth['template_sha256']), 'Completion identity')
        require(truth['template_sha256'] not in templates, f'Duplicate haplotype template on {chrom}')
        require(all(re.fullmatch('[0-9a-f]{64}', truth[k]) for k in ('template_sha256', 'query_sha256')),
                'Invalid sequence identity digest')
        templates.add(truth['template_sha256'])
        nodes.add(done['node'].split('.')[0])
        epoch = done.get('execution_code_identity_sha256') or done.get('pilot_generation_code_identity_sha256')
        require(epoch in epochs, f'Unarchived execution identity: {folder}: {epoch}')
        used_epochs[epoch] += 1
        wanted = {name for name, e in plans.items() if h in e['carriers']}
        require(len(truth['truth']) == len(wanted) and {e['id'] for e in truth['truth']} == wanted,
                f'Truth carrier assignment differs: {folder}')
        for event in truth['truth']:
            fields(event, {k: v for k, v in plans[event['id']].items() if k != 'carriers'}, 'Frozen event design')
            same(hashlib.sha256(event['alt'].encode()).hexdigest(), event['inserted_sha256'], 'Inserted allele digest')
        fields(oracle['metrics'], dict(n_truth=len(wanted), tp_events=len(wanted), fp_calls=0,
                                      unresolved_calls=0, full_sequence_reconstruction=True), 'Oracle metrics')
        parsed = read_json(folder / 'scoring/parsed_calls_production_active.json', evidence)
        ids = [r['id'] for r in parsed]
        require(len(ids) == len(set(ids)), 'Duplicated parsed source call ID')
        source_ids.update(f'H{h:03d}__{i}' for i in ids)
        for p in ('raw', *THRESHOLDS):
            name = 'production_active.json' if p == 'raw' else f'filtered_maf_{p}.json'
            score = read_json(folder / 'scoring' / name, evidence)
            calls = ids if p == 'raw' else retained[h, p]
            call_status[h, p], event_status[h, p] = audit_score(score, wanted, calls, manifest['strict_scoring_parameters'])
            if p == 'raw':
                same(done['metrics'], score['metrics'], 'Completion/raw-score metrics')
                same(done['records'], len(ids), 'Completion/raw record count')
                if score.get('metadata', {}).get('capacity_refinement'):
                    archive = folder / 'scoring_before_capacity_canonicalization'
                    same(score['metadata']['original_score_sha256'], sha(archive / 'production_active.json'), 'Archived raw-score hash')
                    same(done['original_completion_sha256'], sha(archive / 'complete.json'), 'Archived completion hash')
                    same(done['canonical_raw_scoring_parameters'], manifest['strict_scoring_parameters'], 'Canonicalization receipt')
                    canonical = read_json(root / 'raw_capacity_canonicalization.json', evidence)
                    same(canonical['cases'][truth['id']]['new_metrics'], score['metrics'], 'Canonicalization/root receipt')
    require(len(templates) == 256 and member_ids == source_ids, 'Incomplete template/source membership set')
    winners = {}
    for row in rows:
        group = members_by_row[row['row_index']]
        require(set(row['source_ids']) == {x['source_id'] for x in group}, 'Row/source membership mismatch')
        require(len({x['haplotype_index'] for x in group}) == len(group) == row['ac'], 'Row AC mismatch')
        fields(row, dict(an=n, af=row['ac']/n, maf=min(row['ac']/n, 1-row['ac']/n),
                         frequency_bin=frequency_bin(min(row['ac']/n, 1-row['ac']/n))), 'Row frequency')
        raw = [call_status[x['haplotype_index'], 'raw'][x['source_id'].split('__', 1)[1]] for x in group]
        fields(row, {'raw_' + k: v for k, v in classify(raw).items()}, 'Raw row classifications')
        passes = row['ac'] >= 2 and row['ac'] <= int(.99*n)
        same(row['passes_ac_af'], passes, 'Production AC filter')
        if passes:
            key = (row['pos'], row['end'], row['svtype'], row['svlen'])
            if key not in winners or row['ac'] > rows[winners[key]]['ac']:
                winners[key] = row['row_index']
        for p in THRESHOLDS:
            survives = row['dedup_retained'] and row['maf'] >= p/100
            require((str(p) in row['rowstatus_by_maf']) == survives, 'Missing/extraneous filtered row status')
            if survives:
                statuses = [call_status[x['haplotype_index'], p][x['source_id'].split('__', 1)[1]] for x in group]
                fields(row['rowstatus_by_maf'][str(p)], classify(statuses), 'Filtered row classification')
                if p == 0:
                    fields(row, classify(statuses), 'Baseline row classification')
            for member in group:
                expected = (call_status[member['haplotype_index'], p][member['source_id'].split('__', 1)[1]]
                            if survives else None)
                same(member[f'status_maf_{p}'], expected, 'Filtered membership status')
    require({r['row_index'] for r in rows if r['dedup_retained']} == set(winners.values()), 'Deduplication differs')
    expected_pairs = {(name, h) for name, e in plans.items() for h in e['carriers']}
    observed_pairs = {(e['event_id'], e['haplotype_index']) for e in analysis['event_carriers']}
    require(len(analysis['event_carriers']) == len(expected_pairs) and observed_pairs == expected_pairs,
            'Truth-event carrier denominator incomplete')
    for row in analysis['event_carriers']:
        name, h = row['event_id'], row['haplotype_index']
        e = plans[name]
        fields(row, dict(true_ac=len(e['carriers']), true_af=len(e['carriers'])/n,
                         recipe=e['recipe'], true_type=e['kind'], size=e['size']), 'True event frequency')
        fields(row, dict(status=event_status[h, 'raw'][name]['status'],
                         raw_recovered=event_status[h, 'raw'][name]['recovered']), 'Raw event recovery')
        for p in THRESHOLDS:
            fields(row, {f'status_maf_{p}': event_status[h, p][name]['status'],
                         f'recovered_maf_{p}': event_status[h, p][name]['recovered']}, 'Filtered event recovery')
    for name, values in (('merged_rows.tsv', rows), ('source_membership.tsv', members),
                         ('event_carriers.tsv', analysis['event_carriers'])):
        table_matches(root / 'analysis' / chrom / name, values, evidence)
    return analysis, nodes, used_epochs


def audit_aggregate(root, manifest, analyses, evidence):
    folder = root / 'analysis/merged'
    summary = read_json(folder / 'summary.json', evidence)
    rows = [r for a in analyses for r in a['rows'] if r['dedup_retained']]
    events = [r for a in analyses for r in a['event_carriers']]
    fields(summary, dict(chromosomes=manifest['config']['chromosomes'], n_haplotypes=256,
                         baseline_retained_rows=len(rows), truth_event_count=128), 'Aggregate summary')
    expected_keys = set(itertools.product(('all', 'af_le_0.5'), ['ALL', *manifest['config']['chromosomes']],
                                         (0., .01, .05, .1), BINS))
    observed = set()
    for group in summary['frequency_metrics']:
        key = (group['population'], group['chrom'], group['maf_min'], group['frequency_bin'])
        require(key not in observed, 'Duplicated aggregate frequency group')
        observed.add(key)
        population, chrom, threshold, bin_name = key
        selected = [dict(r, **r['rowstatus_by_maf'][str(int(threshold*100))]) for r in rows
                    if (chrom == 'ALL' or r['chrom'] == chrom) and r['maf'] >= threshold
                    and (population == 'all' or r['af'] <= .5) and r['frequency_bin'] == bin_name]
        fields(group, group_metrics(selected), 'Aggregate frequency metrics')
    require(observed == expected_keys, 'Incomplete frequency/threshold group coverage')
    expected_enrichment = set(itertools.product(('all', 'af_le_0.5'),
        ['ALL', *manifest['config']['chromosomes']], (.01, .05), ('all_false', 'any_false', 'mixed')))
    seen_enrichment = set()
    for entry in summary['rare_vs_common']:
        key = entry['population'], entry['chrom'], entry['rare_maf_below'], entry['endpoint']
        require(key not in seen_enrichment, 'Duplicate enrichment comparison')
        seen_enrichment.add(key)
        population, chrom, cutoff, endpoint = key
        universe = [r for r in rows if (chrom == 'ALL' or r['chrom'] == chrom)
                    and (population == 'all' or r['af'] <= .5)]
        rare = [r for r in universe if r['maf'] < cutoff]
        common = [r for r in universe if r['maf'] >= .05]
        def bounds(group):
            if endpoint == 'all_false':
                return sum(r['n_tp'] == 0 and r['n_unresolved'] == 0 for r in group), sum(r['n_tp'] == 0 for r in group)
            if endpoint == 'any_false':
                return sum(r['n_fp'] > 0 for r in group), sum(r['n_fp'] + r['n_unresolved'] > 0 for r in group)
            require(endpoint == 'mixed', 'Unknown enrichment endpoint')
            return sum(mixed(r)[0] for r in group), sum(mixed(r)[1] for r in group)
        rl, ru = bounds(rare)
        cl, cu = bounds(common)
        rlo, rhi, clo, chi = ratio(rl, len(rare)), ratio(ru, len(rare)), ratio(cl, len(common)), ratio(cu, len(common))
        fields(entry, dict(rare_rows=len(rare), common_rows=len(common), common_maf_at_least=.05,
            rare_lower_count=rl, rare_upper_count=ru, common_lower_count=cl, common_upper_count=cu,
            rare_lower=rlo, rare_upper=rhi, common_lower=clo, common_upper=chi,
            difference_lower=rlo-chi if rlo is not None and chi is not None else None,
            difference_upper=rhi-clo if rhi is not None and clo is not None else None,
            ratio_lower=rlo/chi if rlo is not None and chi else None,
            ratio_upper=rhi/clo if rhi is not None and clo else None), 'Rare/common enrichment')
    require(seen_enrichment == expected_enrichment, 'Incomplete enrichment comparison coverage')
    recovery = defaultdict(list)
    for row in events:
        for p in ('raw', *THRESHOLDS):
            key = (row['chrom'], row['recipe'], row['true_type'], row['size'], row['true_ac'],
                   'raw' if p == 'raw' else p/100)
            recovery[key].append((row['raw_recovered'], row['status']) if p == 'raw' else
                                 (row[f'recovered_maf_{p}'], row[f'status_maf_{p}']))
    seen = set()
    for row in summary['truth_recovery']:
        key = tuple(row[k] for k in ('chrom', 'recipe', 'true_type', 'size', 'true_ac', 'filter'))
        require(key in recovery and key not in seen, 'Invalid aggregate recovery group')
        seen.add(key)
        values = recovery[key]
        t, u = sum(x[0] for x in values), sum(x[1] == 'unresolved' for x in values)
        fields(row, dict(carriers=len(values), recovered=t, unresolved=u, lower=t/len(values),
                         upper=(t+u)/len(values)), 'Aggregate recovery')
    require(seen == set(recovery), 'Missing recovery group')
    require(len(events) == 4080, 'Final true-carrier denominator is not 4080')
    for name, expected in (('frequency_metrics.tsv', summary['frequency_metrics']),
                           ('truth_carrier_recovery.tsv', summary['truth_recovery']),
                           ('rare_vs_common.tsv', summary['rare_vs_common']),
                           ('retained_rows.tsv', rows), ('event_carriers.tsv', events)):
        table_matches(folder / name, expected, evidence)
    with (folder / 'truth_event_retention.tsv').open() as handle:
        retention = list(csv.DictReader(handle, delimiter='\t'))
    by_event = defaultdict(list)
    for row in events:
        by_event[row['chrom'], row['event_id']].append(row)
    require(len(retention) == len(by_event) == 128, 'Event retention denominator differs')
    seen = set()
    for row in retention:
        key = row['chrom'], row['event_id']
        require(key in by_event and key not in seen, 'Invalid retention event')
        seen.add(key)
        values = by_event[key]
        same(int(row['carriers']), len(values), 'Retention true carriers')
        same(int(row['true_ac']), values[0]['true_ac'], 'Retention true AC')
        same(float(row['true_af']), values[0]['true_af'], 'Retention true AF')
        raw_count = sum(e['raw_recovered'] for e in values)
        same(int(row['raw_recovered_carriers']), raw_count, 'Raw recovered carriers')
        same(row['raw_event_discovered'], str(raw_count > 0), 'Raw discovered event')
        same(int(row['unresolved_carriers']), sum(e['status'] == 'unresolved' for e in values), 'Raw unresolved carriers')
        for p in THRESHOLDS:
            count = sum(e[f'recovered_maf_{p}'] for e in values)
            same(int(row[f'recovered_carriers_maf_{p}']), count, 'Retained true carriers')
            same(row[f'event_discovered_maf_{p}'], str(count > 0), 'Event discovery')
    evidence[str(folder / 'truth_event_retention.tsv')] = sha(folder / 'truth_event_retention.tsv')
    return dict(baseline_rows=len(rows), truth_events=len(by_event), truth_carrier_occurrences=len(events))


def run_audit(root, output, semantic_receipt=None):
    root, output = Path(root).resolve(), Path(output).resolve()
    receipt = dict(passed=False, started_utc=datetime.now(timezone.utc).isoformat(),
                   root=str(root), status='audit_in_progress', source_sha256={})
    output.parent.mkdir(parents=True, exist_ok=True)
    def save():
        output.write_text(json.dumps(receipt, indent=2, allow_nan=False) + '\n')
    save()
    try:
        evidence = receipt['source_sha256']
        manifest = read_json(root / 'manifest.json', evidence)
        require(manifest['config']['n_haplotypes'] == 256 and set(manifest['config']['chromosomes']) == CHROMS,
                'Expected the frozen four-chromosome, 256-haplotype design')
        recipe_counts = defaultdict(Counter)
        for plan in manifest['plans'].values():
            require(len(plan['events']) == 32, 'Expected 32 truth sites per chromosome')
            for event in plan['events']:
                require(len(event['carriers']) == len(set(event['carriers'])) == event['carrier_count']
                        and all(isinstance(h, int) and 0 <= h < 256 for h in event['carriers']), 'Invalid fixed carrier set')
                recipe_counts[event['recipe']][event['carrier_count']] += 1
        require(len(recipe_counts) == 16 and all(c == Counter([1, 2, 4, 8, 16, 32, 64, 128])
                                                for c in recipe_counts.values()), 'Unbalanced recipe/frequency design')
        missing = [f'{chrom}/H{h:03d}' for chrom in manifest['config']['chromosomes'] for h in range(256)
                   if not (root / 'chromosomes' / chrom / 'haps' / f'H{h:03d}' / 'complete.json').exists()]
        require(not missing, f'Incomplete cohort: {len(missing)} of 1024 cases missing; examples={missing[:8]}')
        gate = read_json(root / 'batch_equivalence.json', evidence)
        require(gate.get('passed') is True, 'Batch/individual mapping gate did not pass')
        epochs = identity_catalog(root, manifest, evidence)
        analyses, nodes, used_epochs = [], set(), Counter()
        for chrom in manifest['config']['chromosomes']:
            analysis, found_nodes, found_epochs = audit_chromosome(root, chrom, manifest, epochs, evidence)
            analyses.append(analysis)
            nodes.update(found_nodes)
            used_epochs.update(found_epochs)
        receipt['resources'] = audit_resources(read_json(root / 'submissions.json', evidence), manifest['config'], nodes)
        receipt['reconciliation'] = audit_aggregate(root, manifest, analyses, evidence)
        gate_path = Path(semantic_receipt) if semantic_receipt else root / 'analysis/merged_semantics_validation.json'
        semantics = read_json(gate_path, evidence)
        require(semantics.get('passed') is True, 'Merged symbolic/BND semantic gate failed')
        same(semantics['manifest_sha256'], sha(root / 'manifest.json'), 'Semantic receipt manifest identity')
        same(semantics['validator_sha256'], sha(Path(__file__).with_name('validate_cohort_merged_semantics.py')),
             'Semantic validator code identity')
        require(len(semantics['chromosomes']) == 4 and {r['chrom'] for r in semantics['chromosomes']} == CHROMS,
                'Semantic gate does not cover all four chromosomes')
        for item in semantics['chromosomes']:
            require(item['passed'] is True and item['n_haplotypes'] == 256, 'Incomplete semantic chromosome gate')
            matching_analysis = next(a for a in analyses if a['chrom'] == item['chrom'])
            same(Path(item['merged_vcf']).resolve(), Path(matching_analysis['merged_vcf']).resolve(),
                 'Semantic gate used a different merged VCF')
            same(sha(item['merged_vcf']), item['merged_vcf_sha256'], 'Stale merged semantic receipt')
            fields(item['counts'], dict(merged_rows=len(matching_analysis['rows']),
                                        source_carrier_assertions=len(matching_analysis['memberships'])),
                   'Semantic gate record counts')
            require(len(item['source_vcf_sha256']) == 256, 'Incomplete semantic source coverage')
            expected_sources = {str((root / 'chromosomes' / item['chrom'] / 'haps' / f'H{h:03d}' / 'merge_input.vcf.gz').resolve())
                                for h in range(256)}
            require({str(Path(p).resolve()) for p in item['source_vcf_sha256']} == expected_sources,
                    'Semantic gate used a different source callset')
            for path, expected in item['source_vcf_sha256'].items():
                same(sha(path), expected, 'Stale source semantic receipt')
        receipt.update(passed=True, status='complete', completed_cases=1024,
                       distinct_templates_per_chromosome=256, execution_epochs=dict(used_epochs),
                       finished_utc=datetime.now(timezone.utc).isoformat(), auditor_sha256=sha(__file__))
        save()
        return receipt
    except Exception as error:
        receipt.update(passed=False, status='failed', error=f'{type(error).__name__}: {error}',
                       finished_utc=datetime.now(timezone.utc).isoformat())
        save()
        raise AuditFailure(receipt['error']) from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT / 'results/cohort_frequency_v1')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--semantic-receipt', type=Path)
    args = parser.parse_args()
    output = args.output or args.root / 'analysis/final_run_audit.json'
    result = run_audit(args.root, output, args.semantic_receipt)
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}))


if __name__ == '__main__':
    main()
