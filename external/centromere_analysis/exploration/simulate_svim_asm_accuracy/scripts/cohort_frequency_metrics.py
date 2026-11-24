#!/usr/bin/env python3
"""Merge called haplotypes and measure production-row errors by frequency."""
import argparse
from collections import Counter, defaultdict
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from vcf_literal_equivalence import exact_padding_equivalence
from evaluate_multichrom_sv import _write_json
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BCFTOOLS = Path(CENSV_ROOT + '/.pixi/envs/default/bin/bcftools')
THRESHOLDS = (0.0, .01, .05, .10)


def sample_name(index):
    return f'SIM{index // 2:03d}_hap{index % 2 + 1}'


def presence(gt):
    return bool(gt) and all(x is not None for x in gt) and any(x > 0 for x in gt)


def row_class(statuses):
    counts = Counter(statuses)
    if not counts or set(counts) - {'tp', 'fp', 'unresolved'}:
        raise ValueError('Missing or invalid strict status')
    if counts['unresolved']:
        return 'any_unresolved'
    if counts['tp'] and counts['fp']:
        return 'mixed_tp_fp'
    return 'all_tp' if counts['tp'] else 'all_fp'


def mixed_bounds(row):
    tp, fp, unknown = row['n_tp'], row['n_fp'], row['n_unresolved']
    lower = tp > 0 and fp > 0
    upper = tp + fp + unknown >= 2 and tp + unknown >= 1 and fp + unknown >= 1
    return bool(lower), bool(upper)


def filter_and_deduplicate(rows, n):
    best = {}
    for row in rows:
        row['passes_ac_af'] = row['ac'] >= 2 and row['ac'] <= int(n * .99)
        row['dedup_retained'] = False
        if not row['passes_ac_af']:
            continue
        key = (row['pos'], row['end'], row['svtype'], row['svlen'])
        if key not in best or row['ac'] > rows[best[key]]['ac']:
            best[key] = row['row_index']
    for index in best.values():
        rows[index]['dedup_retained'] = True
    return rows


def frequency_bin(maf):
    return '<1%' if maf < .01 else '1-<5%' if maf < .05 else '5-<10%' if maf < .10 else '>=10%'


def _write_tsv(path, rows, fields=()):
    cols = list(fields)
    for row in rows:
        for key in row:
            if key not in cols:
                cols.append(key)
    with path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, separators=(',', ':'), sort_keys=True) if isinstance(v, (dict, list, tuple)) else v for k, v in row.items()})


def _read(path):
    return json.loads(path.read_text())


def validate_manifest(m):
    n = int(m['config']['n_haplotypes'])
    chromosomes = m['config']['chromosomes']
    if len(chromosomes) != len(set(chromosomes)) or set(chromosomes) != set(m['plans']):
        raise ValueError('Chromosome plans are incomplete or duplicated')
    recipe_counts = defaultdict(Counter)
    chromosome_recipe_counts = defaultdict(Counter)
    for chrom in chromosomes:
        events = m['plans'][chrom]['events']
        ids = [e['id'] for e in events]
        if len(ids) != len(set(ids)):
            raise ValueError(f'Duplicate truth event ID on {chrom}')
        for e in events:
            carriers = e['carriers']
            if len(set(carriers)) != len(carriers) or len(carriers) != e['carrier_count'] or any(not isinstance(h, int) or not 0 <= h < n for h in carriers):
                raise ValueError(f'Invalid carrier set: {chrom} {e["id"]}')
            recipe_counts[e['recipe']][e['carrier_count']] += 1
            chromosome_recipe_counts[chrom][e['recipe']] += 1
    if n == 256 and len(chromosomes) == 4:
        expected = {1, 2, 4, 8, 16, 32, 64, 128}
        if len(recipe_counts) != 16 or any(set(counts) != expected or any(v != 1 for v in counts.values()) for counts in recipe_counts.values()):
            raise ValueError('Recipe carrier-frequency allocation is not balanced over eight tiers')
        if any(len(counts) != 16 or any(v != 2 for v in counts.values()) for counts in chromosome_recipe_counts.values()):
            raise ValueError('Each chromosome must have two sites per recipe')
    return dict(chromosomes=len(chromosomes), events=sum(sum(c.values()) for c in chromosome_recipe_counts.values()),
                recipes=len(recipe_counts), recipe_frequency_counts={k: dict(v) for k, v in recipe_counts.items()})


def _source_calls(root, chrom, n):
    sources, score_events = {}, {}
    for h in range(n):
        folder = root / 'chromosomes' / chrom / 'haps' / f'H{h:03d}'
        if not (folder / 'complete.json').exists():
            raise FileNotFoundError(f'Incomplete haplotype {folder}')
        truth, score = _read(folder / 'truth.json'), _read(folder / 'scoring/production_active.json')
        parsed = _read(folder / 'scoring/parsed_calls_production_active.json')
        if truth['chrom'] != chrom or truth['haplotype_index'] != h or truth['sample_name'] != sample_name(h):
            raise ValueError(f'Incorrect haplotype identity: {folder}')
        if len(score['calls']) != len(parsed):
            raise ValueError(f'Call table length mismatch: {folder}')
        for call, full in zip(score['calls'], parsed):
            if call['id'] != full['id'] or call['status'] not in ('tp', 'fp', 'unresolved'):
                raise ValueError(f'Call table disagreement: {folder}')
            identity = f'H{h:03d}__{call["id"]}'
            if identity in sources:
                raise ValueError(f'Duplicate source record {identity}')
            sources[identity] = dict(h=h, sample=sample_name(h), status=call['status'],
                matched_truth_ids=call.get('matched_truth_ids', []), raw_id=call['id'], full=full)
        score_events[h] = {e['id']: e for e in score['events']}
        truth_ids = {e['id'] for e in truth['truth']}
        if truth_ids != set(score_events[h]):
            raise ValueError(f'Event table disagreement: {folder}')
    return sources, score_events


def _validate_target_sources(chrom, sources):
    off_target = [identity for identity, item in sources.items()
                  if item['full'].get('reported_chrom', item['full'].get('chrom')) != chrom]
    if off_target:
        raise ValueError(f'{chrom}: {len(off_target)} production-active source call(s) assigned '
                         f'to another chromosome (examples: {off_target[:3]}). '
                         'Joint all-chromosome analysis is required; refusing target-only merge.')
    return 0


def _source_vcfs(root, chrom, n, sources):
    paths = []
    for h in range(n):
        path = root / 'chromosomes' / chrom / 'haps' / f'H{h:03d}' / 'merge_input.vcf.gz'
        paths.append(path)
        with pysam.VariantFile(str(path)) as vcf:
            if list(vcf.header.samples) != [sample_name(h)]:
                raise ValueError(f'Wrong VCF sample in {path}')
            seen = set()
            for rec in vcf:
                identity = rec.id
                if not identity or identity not in sources or sources[identity]['h'] != h or identity in seen:
                    raise ValueError(f'Invalid merge-input ID {identity} in {path}')
                gt = rec.samples[sample_name(h)].get('GT')
                if not presence(gt):
                    raise ValueError(f'Non-carrier source call {identity}')
                alleles = {x for x in gt if x > 0}
                if len(alleles) != 1:
                    raise ValueError(f'Ambiguous source genotype for {identity}')
                allele = next(iter(alleles))
                if not rec.alts or allele > len(rec.alts):
                    raise ValueError(f'Invalid source ALT index for {identity}')
                sources[identity].update(source_alt=rec.alts[allele-1], source_ref=rec.ref, source_start=rec.start)
                seen.add(identity)
            expected = {k for k, v in sources.items() if v['h'] == h}
            if seen != expected:
                raise ValueError(f'Merge-input and scored calls differ in {path}')
    return paths


def _scalar(value, field):
    if isinstance(value, (tuple, list)):
        unique = set(value)
        if len(unique) != 1:
            raise ValueError(f'Multivalued {field} is unsupported by production row parser: {value}')
        value = value[0]
    return value


def _collect_merged(path, chrom, n, sources):
    rows, memberships, observed_ids = [], [], set()
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        expected_samples = [sample_name(h) for h in range(n)]
        if len(samples) != n or set(samples) != set(expected_samples):
            raise ValueError('Merged sample names do not match the cohort')
        for rec in vcf:
            if rec.chrom != chrom:
                continue
            ids = rec.id.split(';') if rec.id else []
            if not ids or len(set(ids)) != len(ids) or any(i not in sources for i in ids):
                raise ValueError(f'Untraceable merged ID: {rec.chrom}:{rec.pos} {rec.id}')
            if observed_ids.intersection(ids):
                raise ValueError(f'Source ID repeated across merged rows: {rec.id}')
            observed_ids.update(ids)
            carriers = {s for s in samples if presence(rec.samples[s].get('GT'))}
            mapped = [sources[i] for i in ids]
            if len({s['sample'] for s in mapped}) != len(mapped) or {s['sample'] for s in mapped} != carriers:
                raise ValueError(f'Merged carrier/source identity mismatch: {rec.chrom}:{rec.pos} {rec.id}')
            for source in mapped:
                if 'source_alt' not in source:
                    continue
                gt = rec.samples[source['sample']].get('GT')
                alleles = {x for x in gt if x is not None and x > 0}
                if len(alleles) != 1 or not rec.alts or max(alleles) > len(rec.alts):
                    raise ValueError(f'Ambiguous merged ALT at {rec.chrom}:{rec.pos}')
                merged_alt = rec.alts[next(iter(alleles))-1]
                original_alt = source['source_alt']
                if any(x in original_alt + merged_alt for x in '<>[]'):
                    same = original_alt == merged_alt and source['source_start'] == rec.start
                else:
                    same = exact_padding_equivalence(source['source_start'], source['source_ref'], original_alt,
                                                     rec.start, rec.ref, merged_alt)
                if not same:
                    raise ValueError(f'Merged genotype points to wrong source allele at {rec.chrom}:{rec.pos} sample {source["sample"]}')
            ac = len(carriers)
            if not ac:
                raise ValueError('Merged row contains no nonreference genotype')
            types = [s['status'] for s in mapped]
            svlen = _scalar(rec.info.get('SVLEN', 0) or 0, 'SVLEN')
            svtype = _scalar(rec.info.get('SVTYPE', 'UNK'), 'SVTYPE')
            info_text = str(rec).split('\t')[7]
            end_field = next((x[4:] for x in info_text.split(';') if x.startswith('END=')), None)
            end = int(end_field) if end_field is not None else rec.pos
            row = dict(row_index=len(rows), chrom=rec.chrom, pos=rec.pos,
                end=end, svtype=svtype, svlen=abs(svlen),
                ref=rec.ref, alt=list(rec.alts or ()), source_ids=ids,
                ac=ac, an=n, af=ac/n, maf=min(ac/n, 1-ac/n),
                nonreference_major=ac/n > .5, frequency_bin=frequency_bin(min(ac/n, 1-ac/n)),
                classification=row_class(types), n_tp=types.count('tp'), n_fp=types.count('fp'),
                n_unresolved=types.count('unresolved'))
            rows.append(row)
            for identity in ids:
                s = sources[identity]
                memberships.append(dict(chrom=chrom, row_index=row['row_index'], source_id=identity,
                    haplotype_index=s['h'], sample_name=s['sample'], status=s['status'],
                    matched_truth_ids=s['matched_truth_ids']))
    if all('full' in item for item in sources.values()):
        expected = {identity for identity, item in sources.items()
                    if item['full'].get('reported_chrom', item['full'].get('chrom')) == chrom}
        if observed_ids != expected:
            raise ValueError(f'Merged target source membership incomplete: missing={len(expected-observed_ids)}, unexpected={len(observed_ids-expected)}')
    return filter_and_deduplicate(rows, n), memberships


def _apply_filtered_scores(rows, memberships, filtered_scores):
    for row in rows:
        row['raw_classification'] = row['classification']
        row['raw_n_tp'], row['raw_n_fp'], row['raw_n_unresolved'] = row['n_tp'], row['n_fp'], row['n_unresolved']
        row['rowstatus_by_maf'] = {}
    grouped = defaultdict(list)
    for member in memberships:
        grouped[member['row_index']].append(member)
    for threshold in (0, 1, 5, 10):
        per_hap = filtered_scores[threshold]
        lookups = {int(h): {c['id']: c for c in result['calls']} for h, result in per_hap.items()}
        for row in rows:
            survives = row['dedup_retained'] and row['maf'] >= threshold / 100
            statuses = []
            for member in grouped[row['row_index']]:
                key = f'status_maf_{threshold}'
                member[key] = None
                if not survives:
                    continue
                raw_id = member['source_id'].split('__', 1)[1]
                status = lookups[member['haplotype_index']].get(raw_id)
                if status is None:
                    raise ValueError(f'Filtered score lacks surviving call: {member["source_id"]}')
                member[key] = status['status']
                statuses.append(status['status'])
            if survives:
                if len(statuses) != row['ac']:
                    raise ValueError('Filtered carrier/status count mismatch')
                counts = Counter(statuses)
                current = dict(classification=row_class(statuses), n_tp=counts['tp'], n_fp=counts['fp'],
                               n_unresolved=counts['unresolved'])
                row['rowstatus_by_maf'][str(threshold)] = current
                if threshold == 0:
                    row.update(current)
    return rows, memberships


def _event_recovery(chrom, plans, n, score_events, memberships, rows, filtered_scores=None):
    source_row = {m['source_id']: rows[m['row_index']] for m in memberships}
    out = []
    for event in plans['events']:
        carriers = set(event['carriers'])
        if len(carriers) != event['carrier_count'] or any(not 0 <= x < n for x in carriers):
            raise ValueError(f'Invalid truth carrier set for {chrom} {event["id"]}')
        for h in range(n):
            e = score_events[h].get(event['id'])
            if h in carriers and e is None:
                raise ValueError(f'Missing carrier truth: {chrom} {event["id"]} H{h:03d}')
            if h not in carriers and e is not None:
                raise ValueError(f'Unexpected carrier truth: {chrom} {event["id"]} H{h:03d}')
            if h not in carriers:
                continue
            matched = [f'H{h:03d}__{x}' for x in e.get('matched_call_ids', [])]
            if e.get('recovered') and not matched:
                raise ValueError('Recovered event lacks matched call IDs')
            if any(x not in source_row for x in matched):
                raise ValueError(f'Matched call missing from merged target rows: {matched}')
            row = dict(chrom=chrom, event_id=event['id'], recipe=event['recipe'], true_type=event['kind'],
                size=event['size'], true_ac=len(carriers), true_af=len(carriers)/n,
                haplotype_index=h, status=e['status'], raw_recovered=bool(e['recovered']), matched_source_ids=matched)
            for threshold in THRESHOLDS:
                key = int(threshold*100)
                original_retained = bool(e['recovered']) and all(
                    source_row[x]['dedup_retained'] and source_row[x]['maf'] >= threshold for x in matched)
                row[f'original_representation_retained_maf_{key}'] = original_retained
                if filtered_scores is None:
                    row[f'recovered_maf_{key}'] = original_retained
                    row[f'status_maf_{key}'] = e['status']
                else:
                    filtered_events = {x['id']: x for x in filtered_scores[key][h]['events']}
                    if event['id'] not in filtered_events:
                        raise ValueError(f'Filtered score lacks truth event: {chrom} {event["id"]} H{h:03d}')
                    filtered = filtered_events[event['id']]
                    row[f'recovered_maf_{key}'] = bool(filtered['recovered'])
                    row[f'status_maf_{key}'] = filtered['status']
            row['matched_called_afs'] = [source_row[x]['af'] for x in matched]
            out.append(row)
    return out


def active_chunk_tasks(root):
    ledger_path=Path(root)/'submissions.json'
    if not ledger_path.exists():
        raise RuntimeError('Submission ledger missing; refusing recovery without scheduler provenance')
    ledger=_read(ledger_path)
    roots={str(s['job_id']).split('_')[0] for s in ledger['submissions'] if s.get('indices') is not None}
    user=os.environ.get('SLURM_JOB_USER') or os.environ.get('USER') or 'jwein22'
    output=subprocess.check_output(['squeue','--array','-h','-u',user,'-o','%F|%K|%T'],text=True)
    records=[line.split('|') for line in output.splitlines() if line.strip()]
    own=os.environ.get('SLURM_JOB_ID')
    if own and not any(r[0]==own for r in records):
        raise RuntimeError('Scheduler did not report the active analysis allocation; refusing recovery')
    result=set()
    for parent,task,state in records:
        if parent in roots:
            if not task.isdigit():raise RuntimeError('Cannot identify an active cohort array index')
            result.add(int(task))
    return result


def ensure_complete_cohort(manifest_path):
    manifest_path=Path(manifest_path);root=manifest_path.parent;m=_read(manifest_path)
    n=int(m['config']['n_haplotypes']);chunk=int(m['config'].get('chunk_size',n))
    with (root/'.cohort_recovery.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        while True:
            missing=[]
            for ci,chrom in enumerate(m['config']['chromosomes']):
                for start in range(0,n,chunk):
                    absent=[h for h in range(start,min(start+chunk,n)) if not
                        (root/'chromosomes'/chrom/'haps'/f'H{h:03d}'/'complete.json').exists()]
                    if absent:missing.append((ci*(n//chunk)+start//chunk,chrom,len(absent)))
            if not missing:return
            active=active_chunk_tasks(root)
            repair=[item for item in missing if item[0] not in active]
            if not repair:
                print(f'Waiting for {sum(item[2] for item in missing)} haplotypes in live cohort tasks',flush=True)
                time.sleep(30)
                continue
            from cohort_frequency_simulation import run_chunk
            history_path=root/'recovery_within_analysis.json'
            history=_read(history_path) if history_path.exists() else {'attempts':[]}
            for task,chrom,count in repair:
                if task in active_chunk_tasks(root):continue
                failed=sum(item['task']==task and not item['completed'] for item in history['attempts'])
                if failed>=3:
                    raise RuntimeError(f'Recovery task{task} exhausted its global three-attempt budget')
                for attempt in range(failed+1,4):
                    if task in active_chunk_tasks(root):break
                    print(f'Recovering task{task} {chrom}: {count} incomplete haplotypes, attempt{attempt}',flush=True)
                    entry=dict(task=task,chrom=chrom,attempt=attempt,job=os.getenv('SLURM_JOB_ID'),
                               node=os.uname().nodename,incomplete_at_start=count)
                    try:
                        run_chunk(manifest_path,task)
                        entry['completed']=True;history['attempts'].append(entry)
                        _write_json(history_path,history)
                        break
                    except Exception as exc:
                        entry.update(completed=False,error=repr(exc));history['attempts'].append(entry)
                        _write_json(history_path,history)
                        if attempt==3:raise
                        time.sleep(5)


def analyze_chrom(manifest_path, chrom, bcftools=BCFTOOLS):
    from ensure_scipy_assignment import ensure_assignment_available
    ensure_complete_cohort(manifest_path)
    ensure_assignment_available()
    manifest_path = Path(manifest_path); m = _read(manifest_path); root = manifest_path.parent
    n = int(m['config']['n_haplotypes'])
    validate_manifest(m)
    if chrom not in m['config']['chromosomes']:
        raise ValueError(f'Unexpected chromosome: {chrom}')
    sources, score_events = _source_calls(root, chrom, n)
    _validate_target_sources(chrom, sources)
    paths = _source_vcfs(root, chrom, n, sources)
    out = root / 'analysis' / chrom; out.mkdir(parents=True, exist_ok=True)
    listing = out / 'merge_inputs.txt'; listing.write_text(''.join(str(p.resolve()) + '\n' for p in paths))
    merged = out / 'merged.vcf.gz'
    command = [str(bcftools), 'merge', '-l', str(listing), '-Oz', '-o', str(merged)]
    subprocess.run(command, check=True)
    subprocess.run([str(bcftools), 'index', '-f', '-t', str(merged)], check=True)
    rows, memberships = _collect_merged(merged, chrom, n, sources)
    from cohort_frequency_rescore import rescore_thresholds
    filtered_scores = rescore_thresholds(root, chrom, m, rows, memberships, sources)
    rows, memberships = _apply_filtered_scores(rows, memberships, filtered_scores)
    events = _event_recovery(chrom, m['plans'][chrom], n, score_events, memberships, rows, filtered_scores)
    _write_tsv(out / 'merged_rows.tsv', rows, ('row_index', 'chrom'))
    _write_tsv(out / 'source_membership.tsv', memberships, ('chrom', 'row_index', 'source_id'))
    _write_tsv(out / 'event_carriers.tsv', events, ('chrom', 'event_id', 'haplotype_index'))
    result = dict(chrom=chrom, n_haplotypes=n, rows=rows, memberships=memberships, event_carriers=events,
                  command=command, merged_vcf=str(merged))
    (out / 'analysis.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result


def aggregate(manifest_path):
    manifest_path = Path(manifest_path); m = _read(manifest_path); root = manifest_path.parent
    design_validation = validate_manifest(m)
    analyses = [_read(root / 'analysis' / chrom / 'analysis.json') for chrom in m['config']['chromosomes']]
    all_rows = [r for a in analyses for r in a['rows'] if r['dedup_retained']]
    events = [r for a in analyses for r in a['event_carriers']]
    out = root / 'analysis' / 'merged'; out.mkdir(parents=True, exist_ok=True)
    grouped = []
    for population in ('all', 'af_le_0.5'):
        for chrom in ['ALL'] + list(m['config']['chromosomes']):
            for threshold in THRESHOLDS:
                selected = [dict(r, **r['rowstatus_by_maf'][str(int(threshold*100))]) for r in all_rows
                            if (chrom == 'ALL' or r['chrom'] == chrom)
                            and (population == 'all' or r['af'] <= .5) and r['maf'] >= threshold]
                for bin_name in ('<1%', '1-<5%', '5-<10%', '>=10%'):
                    rows = [r for r in selected if r['frequency_bin'] == bin_name]
                    classes = Counter(r['classification'] for r in rows)
                    n = len(rows); fp = sum(r['n_fp'] for r in rows); tp = sum(r['n_tp'] for r in rows)
                    u = sum(r['n_unresolved'] for r in rows)
                    no_tp_with_unresolved = sum(r['n_tp'] == 0 and r['n_unresolved'] > 0 for r in rows)
                    grouped.append(dict(population=population, chrom=chrom, maf_min=threshold,
                        frequency_bin=bin_name, n_rows=n, all_tp=classes['all_tp'], all_fp=classes['all_fp'],
                        mixed_tp_fp=classes['mixed_tp_fp'], any_unresolved=classes['any_unresolved'],
                        proven_mixed=sum(mixed_bounds(r)[0] for r in rows),
                        possible_mixed=sum(mixed_bounds(r)[1] for r in rows),
                        all_fp_lower=classes['all_fp']/n if n else None,
                        all_fp_upper=(classes['all_fp'] + no_tp_with_unresolved)/n if n else None,
                        any_fp_lower=sum(r['n_fp'] > 0 for r in rows)/n if n else None,
                        any_fp_upper=sum(r['n_fp'] > 0 or r['n_unresolved'] > 0 for r in rows)/n if n else None,
                        carrier_tp=tp, carrier_fp=fp, carrier_unresolved=u,
                        carrier_fp_lower=fp/(tp+fp+u) if tp+fp+u else None,
                        carrier_fp_upper=(fp+u)/(tp+fp+u) if tp+fp+u else None))
    enrichment = []
    for population in ('all', 'af_le_0.5'):
        for chrom in ['ALL'] + list(m['config']['chromosomes']):
            universe = [r for r in all_rows if (chrom == 'ALL' or r['chrom'] == chrom)
                        and (population == 'all' or r['af'] <= .5)]
            common = [r for r in universe if r['maf'] >= .05]
            for cutoff in (.01, .05):
                rare = [r for r in universe if r['maf'] < cutoff]
                for endpoint in ('all_false', 'any_false', 'mixed'):
                    def bounds(group):
                        if endpoint == 'all_false':
                            lower = sum(r['classification'] == 'all_fp' for r in group)
                            upper = sum(r['n_tp'] == 0 for r in group)
                        elif endpoint == 'any_false':
                            lower = sum(r['n_fp'] > 0 for r in group)
                            upper = sum(r['n_fp'] > 0 or r['n_unresolved'] > 0 for r in group)
                        else:
                            lower = sum(mixed_bounds(r)[0] for r in group)
                            upper = sum(mixed_bounds(r)[1] for r in group)
                        return lower, upper
                    rl, ru = bounds(rare); cl, cu = bounds(common)
                    rlo = rl/len(rare) if rare else None; rhi = ru/len(rare) if rare else None
                    clo = cl/len(common) if common else None; chi = cu/len(common) if common else None
                    enrichment.append(dict(population=population, chrom=chrom, endpoint=endpoint,
                        rare_maf_below=cutoff, rare_rows=len(rare), rare_lower_count=rl, rare_upper_count=ru,
                        common_maf_at_least=.05, common_rows=len(common), common_lower_count=cl, common_upper_count=cu,
                        rare_lower=rlo, rare_upper=rhi, common_lower=clo, common_upper=chi,
                        difference_lower=rlo-chi if rlo is not None and chi is not None else None,
                        difference_upper=rhi-clo if rhi is not None and clo is not None else None,
                        ratio_lower=rlo/chi if rlo is not None and chi else None,
                        ratio_upper=rhi/clo if rhi is not None and clo else None))
    rates = defaultdict(lambda: [0, 0, 0])
    for row in events:
        for threshold in THRESHOLDS:
            key = (row['chrom'], row['recipe'], row['true_type'], row['size'], row['true_ac'], threshold)
            vals = rates[key]; vals[0] += 1; vals[1] += row[f'recovered_maf_{int(threshold*100)}']; vals[2] += row[f'status_maf_{int(threshold*100)}'] == 'unresolved'
        key = (row['chrom'], row['recipe'], row['true_type'], row['size'], row['true_ac'], 'raw')
        vals = rates[key]; vals[0] += 1; vals[1] += row['raw_recovered']; vals[2] += row['status'] == 'unresolved'
    recovery = [dict(zip(('chrom', 'recipe', 'true_type', 'size', 'true_ac', 'filter'), key), carriers=n,
                     recovered=tp, unresolved=u, lower=tp/n, upper=(tp+u)/n)
                for key, (n, tp, u) in sorted(rates.items(), key=lambda x: str(x[0]))]
    by_event = defaultdict(list)
    for row in events:
        by_event[(row['chrom'], row['event_id'])].append(row)
    truth_events = []
    for (chrom, event_id), event_rows in sorted(by_event.items()):
        first = event_rows[0]
        entry = dict(chrom=chrom, event_id=event_id, recipe=first['recipe'], true_type=first['true_type'],
                     size=first['size'], true_ac=first['true_ac'], true_af=first['true_af'],
                     carriers=len(event_rows), raw_recovered_carriers=sum(r['raw_recovered'] for r in event_rows),
                     unresolved_carriers=sum(r['status'] == 'unresolved' for r in event_rows),
                     matched_called_afs=sorted({af for r in event_rows for af in r['matched_called_afs']}))
        entry['raw_event_discovered'] = entry['raw_recovered_carriers'] > 0
        for threshold in THRESHOLDS:
            suffix = int(threshold * 100)
            count = sum(r[f'recovered_maf_{suffix}'] for r in event_rows)
            entry[f'recovered_carriers_maf_{suffix}'] = count
            entry[f'event_discovered_maf_{suffix}'] = count > 0
        truth_events.append(entry)
    _write_tsv(out / 'frequency_metrics.tsv', grouped)
    _write_tsv(out / 'truth_carrier_recovery.tsv', recovery)
    _write_tsv(out / 'rare_vs_common.tsv', enrichment)
    _write_tsv(out / 'retained_rows.tsv', all_rows)
    _write_tsv(out / 'event_carriers.tsv', events)
    _write_tsv(out / 'truth_event_retention.tsv', truth_events)
    result = dict(chromosomes=m['config']['chromosomes'], n_haplotypes=m['config']['n_haplotypes'],
        baseline_retained_rows=len(all_rows), frequency_metrics=grouped, rare_vs_common=enrichment, truth_recovery=recovery,
        truth_event_count=len(truth_events), design_validation=design_validation,
        interpretation='Row-weighted descriptive results. Mixed and unresolved carrier assertions remain separate; no independent-record confidence intervals.')
    (out / 'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--chrom')
    parser.add_argument('--aggregate', action='store_true')
    parser.add_argument('--bcftools', default=str(BCFTOOLS))
    args = parser.parse_args()
    if bool(args.chrom) == bool(args.aggregate):
        parser.error('Choose exactly one of --chrom or --aggregate')
    result = aggregate(args.manifest) if args.aggregate else analyze_chrom(args.manifest, args.chrom, Path(args.bcftools))
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'memberships', 'event_carriers', 'frequency_metrics', 'truth_recovery')}))
