#!/usr/bin/env python3
"""Type-stratified metrics from frozen complete-truth SV benchmark evidence."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

from simulate_sv_recovery import write_json
from evaluate_complete_truth import table

ROOT = Path(__file__).resolve().parents[1]
TRUTH_TYPES = ['INS', 'DEL', 'TANDEM', 'DISPERSED', 'INV', 'REPLACEMENT']
CALL_TYPES = ['INS', 'DEL', 'DUP:TANDEM', 'DUP:INT', 'INV', 'BND']


def read_table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def yes(value):
    return value is True or str(value).lower() == 'true'


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def size_bin(size):
    if size < 50:
        return '<50'
    if size <= 1000:
        return '50-1000'
    if size <= 10000:
        return '1001-10000'
    if size <= 50000:
        return '10001-50000'
    return '>50000'


def sensitivity(rows):
    total = len(rows)
    tp = sum(yes(row['recovered']) for row in rows)
    unknown = sum(row['status'] == 'unresolved' for row in rows)
    single = [row for row in rows if yes(row['recovered']) and len(json.loads(row['matched_call_ids'])) == 1]
    multiple = sum(yes(row['recovered']) and len(json.loads(row['matched_call_ids'])) > 1 for row in rows)
    return dict(n_truth=total, tp_events=tp, fn_events=total-tp-unknown, unresolved_events=unknown,
                sensitivity=ratio(tp, total) if not unknown else None,
                sensitivity_lower=ratio(tp, total), sensitivity_upper=ratio(tp+unknown, total),
                exact_representation_events=sum(yes(row['exact']) for row in rows),
                single_call_recovered=len(single), multiple_call_recovered=multiple,
                single_call_literal_type_correct=sum(yes(row['literal_type_match']) for row in single),
                single_call_type_compatible=sum(yes(row['type_compatible']) for row in single),
                n_positive_queries=len({row['case'] for row in rows}),
                n_backgrounds=len({row['background'] for row in rows}),
                n_seeds=len({row['seed'] for row in rows}))


def precision(rows):
    total = len(rows)
    tp = sum(row['status'] == 'tp' for row in rows)
    fp = sum(row['status'] == 'fp' for row in rows)
    unknown = total-tp-fp
    return dict(n_calls=total, tp_calls=tp, fp_calls=fp, unresolved_calls=unknown,
                precision=ratio(tp, total) if not unknown else None,
                precision_lower=ratio(tp, total), precision_upper=ratio(tp+unknown, total))


def false_ids(window, calls, kind, minimum_size):
    identifiers = set()
    for ident in json.loads(window['fp_call_ids']):
        call = calls[(window['case'], window['scope'], ident)]
        if (kind == 'ALL' or call['kind'] == kind) and int(call['called_edit_size_bp']) >= minimum_size:
            identifiers.add(ident)
    return identifiers


def control_specificity(windows, calls, kind='ALL', minimum_size=10):
    eligible = [row for row in windows if yes(row['eligible_negative']) and row['stratum'] == 'active_array']
    per_query = defaultdict(set)
    for row in eligible:
        per_query[row['case']].update(false_ids(row, calls, kind, minimum_size))
    full_windows = [row for row in eligible if int(row['end'])-int(row['start']) == 50000]
    fp_windows = sum(bool(false_ids(row, calls, kind, minimum_size)) for row in full_windows)
    uncertain = sum(yes(row['unresolved_window']) for row in full_windows)
    n_fail = sum(bool(ids) for ids in per_query.values())
    n_fp_records = sum(len(ids) for ids in per_query.values())
    exposure = sum(int(row['end'])-int(row['start']) for row in eligible)
    return dict(n_null_queries=len(per_query), null_queries_with_fp=n_fail,
                null_query_specificity=ratio(len(per_query)-n_fail, len(per_query)),
                negative_windows_50kb=len(full_windows), fp_windows_50kb=fp_windows,
                unresolved_windows_50kb=uncertain,
                window_specificity_50kb=ratio(len(full_windows)-fp_windows, len(full_windows)) if not uncertain else None,
                fp_records_in_active_arrays=n_fp_records,
                fp_records_per_null_query=ratio(n_fp_records, len(per_query)),
                active_exposure_bp=exposure, fp_records_per_Mb=ratio(n_fp_records*1e6, exposure))


def safe_window_specificity(windows, calls, kind='ALL', minimum_size=10):
    eligible = [row for row in windows if yes(row['eligible_negative']) and row['stratum'] == 'active_array'
                and int(row['end'])-int(row['start']) == 50000]
    fp = sum(bool(false_ids(row, calls, kind, minimum_size)) for row in eligible)
    unresolved = sum(yes(row['unresolved_window']) for row in eligible)
    return dict(negative_windows_50kb=len(eligible), fp_windows_50kb=fp, unresolved_windows_50kb=unresolved,
                window_specificity_50kb=ratio(len(eligible)-fp, len(eligible)) if not unresolved else None)


def grouped(rows, dimensions):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[key] for key in dimensions)].append(row)
    for key, values in sorted(groups.items(), key=lambda item: str(item[0])):
        yield dict(zip(dimensions, key)), values


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_evidence(source):
    manifest = json.loads((source/'manifest.json').read_text())
    case_map = {case['id']: case for case in manifest['cases']}
    assert len(case_map) == 54
    annotation = json.loads((source/'repeat_annotation.json').read_text())
    repeats = {(row['case'], row['id']): row['composition_class'] for row in annotation['events']}
    events, calls, windows = [read_table(source/name) for name in ['events.tsv', 'calls.tsv', 'negative_windows.tsv']]
    for rows in [events, calls, windows]:
        for row in rows:
            case = case_map[row['case']]
            architecture = case['input_architecture']
            if architecture not in ['full_contig','three_contigs']:
                raise ValueError('Unknown input architecture: '+architecture)
            row['architecture'] = 'full_contig' if architecture == 'full_contig' else 'fragmented'
            row['is_null'] = not case['truth']
            row['array_rate'] = float(row['array_rate'])
            row['seed'] = int(row['seed'])
    for row in events:
        row['truth_type'] = row['kind']
        row['size_class'] = size_bin(int(row['size_bp']))
        row['repeat_class'] = repeats[(row['case'], row['id'])]
        truth = next(event for event in case_map[row['case']]['truth'] if event['id'] == row['id'])
        assert truth['kind'] == row['kind'] and max(len(truth['alt']), truth['end']-truth['start']) == int(row['size_bp'])
    lookup = {(row['case'], row['scope'], row['id']): row for row in calls}
    assert len(lookup) == len(calls)
    for row in calls:
        row['called_type'] = row['kind']
        row['size_class'] = size_bin(int(row['called_edit_size_bp']))
    for row in events:
        matched = json.loads(row['matched_call_ids'])
        if not yes(row['recovered']):
            row['recovery_output'] = 'NO_CERTIFIED_MATCH'
        elif len(matched) > 1:
            row['recovery_output'] = 'MULTIPLE_RECORDS'
        else:
            row['recovery_output'] = lookup[(row['case'], row['scope'], matched[0])]['kind']
    metrics = read_table(source/'metrics.tsv')
    for metric in metrics:
        chosen = [r for r in events if (r['case'], r['scope']) == (metric['case'], metric['scope'])]
        emitted = [r for r in calls if (r['case'], r['scope']) == (metric['case'], metric['scope'])]
        assert len(chosen) == int(metric['n_truth'])
        assert sum(yes(r['recovered']) for r in chosen) == int(metric['tp_events'])
        assert len(emitted) == int(metric['n_calls'])
        assert sum(r['status'] == 'tp' for r in emitted) == int(metric['tp_calls'])
    components = read_table(source/'components.tsv')
    assert not any(yes(row['search_bounded']) for row in components)
    assert not any(row['status'] == 'unresolved' for row in events+calls)
    assert not any(yes(row['unresolved_window']) for row in windows)
    assert not any(row['status'] == 'tp' and int(row['called_edit_size_bp']) < 50 for row in calls)
    assert len(events) == 1188
    return manifest, events, calls, windows, lookup


def run(source, output):
    if (output/'complete.json').exists():
        raise ValueError('Completed metrics already exist; choose a new output directory')
    output.mkdir(parents=True, exist_ok=True)
    manifest, events, calls, windows, lookup = load_evidence(source)
    outputs = {}
    for name, dimensions in [
        ('sensitivity_by_type', ['architecture','scope','truth_type']),
        ('sensitivity_all_query_observations', ['scope','truth_type']),
        ('sensitivity_by_type_divergence', ['architecture','scope','truth_type','array_rate']),
        ('sensitivity_by_type_size_divergence', ['architecture','scope','truth_type','size_class','array_rate']),
        ('sensitivity_by_scenario_type_divergence', ['architecture','scope','scenario','truth_type','array_rate']),
        ('sensitivity_by_repeat_type', ['architecture','scope','truth_type','repeat_class']),
        ('sensitivity_by_seed', ['architecture','scope','scenario','truth_type','array_rate','seed']),
        ('sensitivity_pooled_scenarios_by_seed', ['architecture','scope','truth_type','array_rate','seed'])]:
        outputs[name] = [dict(meta, **sensitivity(rows)) for meta, rows in grouped(events, dimensions)]
    for name, dimensions in [
        ('precision_by_called_type', ['architecture','scope']),
        ('precision_by_called_type_divergence', ['architecture','scope','array_rate']),
        ('precision_by_scenario_type_divergence', ['architecture','scope','scenario','array_rate']),
        ('precision_by_seed', ['architecture','scope','array_rate','seed'])]:
        result = []
        for meta, group_windows in grouped(windows, dimensions):
            selected = [row for row in calls if all(row[key] == value for key,value in meta.items())]
            for minimum in [10,50]:
                for kind in CALL_TYPES+['ALL']:
                    rows = [row for row in selected if (kind == 'ALL' or row['kind'] == kind)
                            and int(row['called_edit_size_bp']) >= minimum]
                    result.append(dict(meta, called_type=kind, reported_record_min_size=minimum,
                                       n_queries=len({r['case'] for r in group_windows}), **precision(rows)))
        outputs[name] = result
    specificity_rows, controls = [], []
    for dimensions in [['architecture','scope'], ['architecture','scope','array_rate'],
                       ['architecture','scope','array_rate','seed']]:
        for meta, rows in grouped(windows, dimensions):
            for minimum in [10,50]:
                for kind in CALL_TYPES+['ALL']:
                    base = dict(meta, called_type=kind, reported_record_min_size=minimum)
                    specificity_rows.append(dict(base, **safe_window_specificity(rows,lookup,kind,minimum)))
                    null = [r for r in rows if r['is_null']]
                    controls.append(dict(base, **control_specificity(null,lookup,kind,minimum)))
    outputs['specificity_truth_free_50kb'] = specificity_rows
    outputs['specificity_no_SV_controls'] = controls
    outputs['type_confusion'] = [dict(meta,n_events=len(rows)) for meta,rows in grouped(events,
                                 ['architecture','scope','truth_type','recovery_output'])]
    outputs['event_level'] = events
    outputs['call_level'] = calls
    for name, rows in outputs.items():
        table(output/(name+'.tsv'), rows)
    primary = lambda row: row.get('architecture')=='full_contig' and row['scope']=='production_active'
    primary_cases=[c for c in manifest['cases'] if c['input_architecture']=='full_contig']
    summary = dict(dataset='CHM13-derived complete-truth benchmark',
                   primary_scope='full_contig / production_active; all divergence settings and scenarios',
                   truth_matching='background-SNP-restored exact sequence consequence; repeat/split equivalence allowed',
                   n_queries=len(primary_cases),n_positive_queries=sum(bool(c['truth']) for c in primary_cases),
                   n_null_queries=sum(not c['truth'] for c in primary_cases),
                   n_seeds=len({r['seed'] for r in events if r['architecture']=='full_contig'}),
                   n_backgrounds=len({c['background'] for c in primary_cases}),
                   sensitivity=[r for r in outputs['sensitivity_by_type'] if primary(r)],
                   precision=[r for r in outputs['precision_by_called_type'] if primary(r) and r['reported_record_min_size']==10],
                   null_specificity=[r for r in controls if primary(r) and 'array_rate' not in r and r['reported_record_min_size']==10],
                   window_specificity=[r for r in specificity_rows if primary(r) and 'array_rate' not in r and r['reported_record_min_size']==10],
                   n_unresolved_event_rows=sum(r['status']=='unresolved' for r in events),
                   TE_rich_insertion_truth_count=sum(r['kind']=='INS' and 'TE' in r['repeat_class'] for r in events if r['scope']=='production_active'),
                   biological_independence='Reused genomes, sites and payloads; descriptive assay fractions, not population estimates or independent-event CIs')
    write_json(output/'summary.json',summary)
    inputs = ['manifest.json','events.tsv','calls.tsv','metrics.tsv','negative_windows.tsv','components.tsv','repeat_annotation.json']
    write_json(output/'provenance.json',dict(source=str(source),inputs_sha256={name:digest(source/name) for name in inputs},
               script_sha256=digest(Path(__file__)),production_rerun=False,slurm_jobs_submitted=0,
               record_size_note='Size strata condition on emitted records; they do not claim an independent caller run at a new size threshold'))
    print(json.dumps(dict(sensitivity=[(r['truth_type'],r['tp_events'],r['n_truth']) for r in summary['sensitivity']],
                         precision=[(r['called_type'],r['tp_calls'],r['n_calls']) for r in summary['precision']],
                         control_specificity=summary['null_specificity']),indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'results/complete_truth_v1')
    parser.add_argument('--output',type=Path,default=ROOT/'results/sv_type_metrics_v1')
    args=parser.parse_args()
    run(args.source.resolve(),args.output.resolve())
