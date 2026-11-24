#!/usr/bin/env python3
"""Describe natural catalog strata beside frozen, truth-scored simulations."""
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

from analyze_array_position_accuracy import position, read_bed, read_tsv, write_tsv
from summarize_panel_frequency import maf_from_export
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(__file__).resolve().parents[1]
FREQUENCIES = ('<1%', '1-<5%', '>=5%')
SIZES = ('<=1 kb', '1-7 kb', '>7 kb', 'Unsized')
DEPTHS = ('0-20%', '20-40%', '40-60%', '60-80%', '80-100%', 'Outside active arrays')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frequency(value):
    return '<1%' if value < Fraction(1, 100) else '1-<5%' if value < Fraction(5, 100) else '>=5%'


def size_band(size):
    size = abs(int(size))
    return 'Unsized' if not size else '<=1 kb' if size <= 1000 else '1-7 kb' if size <= 7000 else '>7 kb'


def boolean(x):
    return str(x).lower() in ('true', '1')


def metric(group, endpoint):
    n = len(group)
    correct = sum(r['status'] == 'tp' if endpoint == 'precision' else boolean(r['recovered']) for r in group)
    unresolved = sum(r['status'] == 'unresolved' for r in group)
    return dict(n=n, correct=correct, unresolved=unresolved,
                lower=correct / n if n else None,
                upper=(correct + unresolved) / n if n else None,
                site_descriptors=len({r['site_id'] for r in group}) if endpoint == 'sensitivity' else None,
                query_cases=len({r['case_id'] for r in group}))


def aggregate(rows, endpoint):
    out = []
    dimensions = ['frequency_bin', 'size_bin', 'relative_bin']
    combinations = [()] + [ks for n in (1, 2, 3) for ks in itertools.combinations(dimensions, n)]
    combinations += [('svtype',), ('rate',), ('chrom',), ('frequency_bin', 'svtype'),
                     ('frequency_bin', 'size_bin', 'svtype'), ('relative_bin', 'size_bin', 'rate')]
    for keys in combinations:
        groups = defaultdict(list)
        for r in rows:
            groups[(r['dataset'], *(r.get(k, 'Not simulated') for k in keys))].append(r)
        for values, group in sorted(groups.items()):
            out.append(dict(dataset=values[0], stratification='+'.join(keys) or 'All',
                            **dict(zip(keys, values[1:])), **metric(group, endpoint)))
    return out


def panel_groups(rows):
    out = []
    for keys in [(), ('frequency_bin',), ('size_bin',), ('relative_bin',),
                 ('frequency_bin', 'size_bin'), ('frequency_bin', 'relative_bin'),
                 ('size_bin', 'relative_bin'), ('frequency_bin', 'size_bin', 'relative_bin'),
                 ('frequency_bin', 'svtype'), ('chrom',)]:
        groups = defaultdict(list)
        for r in rows:
            groups[tuple(r[k] for k in keys)].append(r)
        for values, group in sorted(groups.items()):
            out.append(dict(stratification='+'.join(keys) or 'All', **dict(zip(keys, values)),
                            rows=len(group), all_alpha_rich=sum(boolean(r['all_alleles_alpha_rich']) for r in group),
                            any_TE_rich=sum(int(r['n_TE_rich_alleles']) > 0 for r in group),
                            any_TE_homology=sum(int(r['n_TE_positive_alleles']) > 0 for r in group),
                            annotated=sum(r['annotation_status'] == 'annotated' for r in group)))
    return out


def cohort_rows(analysis, annotations):
    index = {(r['case_id'], r['record_id']): r for r in annotations if r['dataset'] == 'Cohort production'}
    assert len(index) == 7659
    retained = {r['row_index']: r for r in analysis['rows'] if r['dedup_retained']}
    output = []
    counts = Counter()
    for member in analysis['memberships']:
        if member['status_maf_0'] is None:
            continue
        row = retained[member['row_index']]
        h = f"H{member['haplotype_index']:03d}"
        source_id = member['source_id'].split('__', 1)[1]
        r = dict(index[h, source_id])
        assert r['status'] == member['status_maf_0']
        r.update(frequency_bin=frequency(Fraction(min(row['ac'], row['an']-row['ac']), row['an'])),
                 called_maf=row['maf'], merged_row_index=row['row_index'])
        output.append(r)
        counts[row['row_index'], r['status']] += 1
    for row in retained.values():
        for status in ('tp', 'fp', 'unresolved'):
            assert counts[row['row_index'], status] == row['n_'+status]
    assert len(output) == len(index)
    return output, retained


def inspect_natural_examples(rows, source_manifest, output_dir):
    chosen = {}
    for row in rows:
        if row['svtype'] not in ('INS', 'DEL') or row['size_bin'] == 'Unsized':
            continue
        key = (row['frequency_bin'], row['size_bin'], row['relative_bin'])
        chosen.setdefault(key, row)
    by_chrom = defaultdict(dict)
    for r in chosen.values():
        by_chrom[r['chrom']][r['catalog_row_id']] = r
    inspected = []
    for source in source_manifest['sources']:
        path = Path(source['path'])
        chrom = path.parent.name
        todo = by_chrom.get(chrom)
        if not todo:
            continue
        with path.open() as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if int(row['ac']) < 2 or float(row['af']) > .99:
                    continue
                key = ':'.join(row[k] for k in ('chrom', 'pos', 'end', 'svtype', 'svlen'))
                if key not in todo:
                    continue
                annotation = todo.pop(key)
                assert int(row['ac']) == int(annotation['ac'])
                assert maf_from_export(row) == maf_from_export(annotation)
                ref = row['ref'].upper()
                alts = row['alt'].upper().split(',')
                literals = [a for a in alts if a and set(a) <= set('ACGTN')]
                normalized = []
                for alt in literals:
                    left, right = 0, 0
                    while right < min(len(ref), len(alt)) and ref[-1-right] == alt[-1-right]:
                        right += 1
                    r, a = (ref[:-right], alt[:-right]) if right else (ref, alt)
                    while left < min(len(r), len(a)) and r[left] == a[left]:
                        left += 1
                    r, a = r[left:], a[left:]
                    normalized.append(dict(operation='INS' if not r and a else 'DEL' if r and not a else 'REPLACEMENT',
                                           removed_bp=len(r), added_bp=len(a)))
                inspected.append(dict(catalog_row_id=key, source=str(path), **{k:annotation[k] for k in
                    ('frequency_bin', 'size_bin', 'relative_bin', 'composition')},
                    n_alts=len(alts), n_literal_alts=len(literals),
                    ref_sha256=hashlib.sha256(ref.encode()).hexdigest(),
                    alt_field_sha256=hashlib.sha256(row['alt'].encode()).hexdigest(),
                    normalized_alleles=json.dumps(normalized, separators=(',', ':')),
                    af=row['af'], ac=row['ac'], an=row['an']))
                if not todo:
                    break
        if todo:
            raise ValueError(f'Missing source examples: {list(todo)}')
    assert len(inspected) == len(chosen)
    write_tsv(output_dir/'actual_allele_examples.tsv', inspected)
    return len(inspected)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--output-dir', type=Path, default=ROOT/'results/panel_accuracy_assessment_v1')
    args = p.parse_args()
    root, out = args.root.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    csv.field_size_limit(10_000_000)
    bed = Path(CENSV_ROOT + '/resources/chm13v2.0_censat_v2.1.active.bed')
    arrays = read_bed(bed)
    catalog = read_tsv(root/'results/panel_repeat_composition.tsv')
    source_manifest = json.loads((root/'results/panel_sv_size_distribution.json').read_text())
    source_checks = []
    for entry in source_manifest['sources']:
        path = Path(entry['path'])
        source_checks.append(dict(path=str(path), recorded_bytes=entry['bytes'],
                                  current_bytes=path.stat().st_size, byte_count_matches=entry['bytes']==path.stat().st_size))
    assert all(r['byte_count_matches'] for r in source_checks)
    for r in catalog:
        r['frequency_bin'] = frequency(maf_from_export(r))
        r['size_bin'] = size_band(r['svlen']) if r['svtype'] != 'BND' else 'Unsized'
        loc = position(r['chrom'], int(r['pos'])-1, int(r['pos'])-1, arrays, {})
        r.update({k:loc.get(k, '') for k in ('relative_bin', 'relative_depth', 'distance_to_end_bp')})
    assert len(catalog) == 188591
    assert Counter(r['frequency_bin'] for r in catalog) == {'<1%':118220, '1-<5%':62634, '>=5%':7737}
    write_tsv(out/'actual_calls_annotated.tsv', catalog)
    write_tsv(out/'actual_strata.tsv', panel_groups(catalog))
    n_examples = inspect_natural_examples(catalog, source_manifest, out)
    spatial = root/'results/array_position_v1'
    calls, events = read_tsv(spatial/'calls_annotated.tsv'), read_tsv(spatial/'events_annotated.tsv')
    analysis = json.loads((root/'results/cohort_frequency_v1/joint_analysis/analysis.json').read_text())
    cohort, retained = cohort_rows(analysis, calls)
    source_freq = {(r['case_id'], r['record_id']):r['frequency_bin'] for r in cohort}
    calls = [r for r in calls if r['dataset'] != 'Cohort production'] + cohort
    for r in calls:
        if r['dataset'] == 'Cohort MAF>=5%':
            r['frequency_bin'] = source_freq[r['case_id'], r['record_id']]
    event_index = {(f"H{r['haplotype_index']:03d}", r['event_id']):r for r in analysis['event_carriers']}
    for e in events:
        if e['dataset'].startswith('Cohort'):
            truth = event_index[e['case_id'], e['site_id']]
            e['frequency_bin'] = frequency(Fraction(truth['true_ac'], 256))
    write_tsv(out/'simulation_calls.tsv', calls)
    write_tsv(out/'simulation_events.tsv', events)
    write_tsv(out/'simulation_precision.tsv', aggregate(calls, 'precision'))
    write_tsv(out/'simulation_sensitivity.tsv', aggregate(events, 'sensitivity'))
    row_annotations = []
    for row in retained.values():
        loc = position(row['chrom'], row['pos']-1, row['pos']-1, arrays, {})
        row_annotations.append(dict(row_index=row['row_index'], chrom=row['chrom'], svtype=row['svtype'],
            frequency_bin=frequency(Fraction(min(row['ac'],row['an']-row['ac']),row['an'])),
            size_bin=size_band(row['svlen']), relative_bin=loc['relative_bin'],
            n_tp=row['n_tp'], n_fp=row['n_fp'], n_unresolved=row['n_unresolved']))
    row_summary=[]
    for keys in [('frequency_bin',),('frequency_bin','size_bin'),('frequency_bin','size_bin','relative_bin')]:
        groups=defaultdict(list)
        for r in row_annotations:groups[tuple(r[k] for k in keys)].append(r)
        for key, group in sorted(groups.items()):
            row_summary.append(dict(stratification='+'.join(keys),**dict(zip(keys,key)),rows=len(group),
                entirely_false=sum(r['n_tp']==0 and r['n_unresolved']==0 for r in group),
                all_carriers_correct=sum(r['n_fp']==0 and r['n_unresolved']==0 for r in group),
                mixed=sum(r['n_tp']>0 and r['n_fp']>0 for r in group)))
    write_tsv(out/'simulation_row_correctness.tsv',row_summary)
    small = lambda r: r['size_bin'] in ('<=1 kb', '1-7 kb')
    common = [r for r in catalog if r['frequency_bin'] == '>=5%']
    common_small = [r for r in common if small(r)]
    summary = dict(catalog_rows=len(catalog), common_rows=len(common),
        common_small_rows=len(common_small),
        common_small_outer_fifth=sum(r['relative_bin']=='0-20%' for r in common_small),
        common_small_outer_two_fifths=sum(r['relative_bin'] in ('0-20%', '20-40%') for r in common_small),
        all_source_allele_examples=n_examples,
        literal_ALT_examples=sum(int(r['n_literal_alts']) for r in read_tsv(out/'actual_allele_examples.tsv')),
        new_slurm_jobs=0, simulations_are_not_native_truth_calibration=True)
    for dataset in ('Cohort production', 'Cohort MAF>=5%', '23-chromosome', 'Chr1 intact'):
        c = [r for r in calls if r['dataset']==dataset and small(r)]
        e = [r for r in events if r['dataset']==dataset and small(r)]
        summary[dataset] = dict(small_record_correct=sum(r['status']=='tp' for r in c),
            small_record_total=len(c), small_truth_recovered=sum(boolean(r['recovered']) for r in e),
            small_truth_total=len(e))
    (out/'evidence_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    inputs=[root/'results/panel_repeat_composition.tsv', root/'results/panel_repeat_composition.json',
            root/'results/panel_sv_size_distribution.json', spatial/'calls_annotated.tsv', spatial/'events_annotated.tsv',
            root/'results/cohort_frequency_v1/joint_analysis/analysis.json',
            root/'results/cohort_frequency_v1/joint_execution_code_identity.json',bed]
    provenance=dict(inputs=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in inputs],
        upstream_sources=source_checks,script_sha256=sha(__file__),n_examples=n_examples,
        coordinate_conventions='Catalog rows and merged rows: reported POS-1 as point. Simulation calls: normalized edit midpoint. Truth: actual edit midpoint.',
        projection='No natural call receives a simulated truth label or a calibrated probability.',
        n_new_alignments=0,n_new_slurm_jobs=0)
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(dict(actual_rows=len(catalog),allele_examples=n_examples,simulated_calls=len(calls),simulated_events=len(events))))


if __name__ == '__main__':
    main()
