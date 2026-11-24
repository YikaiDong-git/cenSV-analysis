#!/usr/bin/env python3
"""Join existing allele/reference repeat annotations to retained panel rows."""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys

import pysam

from annotate_simulation_repeats import MEI, UPSTREAM, coverage, file_digest, parse_rm, union_bp


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def overlap(intervals, ends, start, end):
    index = bisect_right(ends, start)
    total = 0
    for left, right in intervals[index:]:
        if left >= end:
            break
        total += max(0, min(end, right) - max(start, left))
    return total


def label(length, groups, status):
    if status != 'annotated':
        return status
    te = groups.get('transposable_element', 0)
    alpha = groups.get('alpha_satellite', 0)
    sat = alpha + groups.get('other_satellite', 0)
    if te >= .8 * length:
        return 'TE_rich'
    if alpha >= .8 * length:
        return 'alpha_rich_with_TE' if te else 'alpha_rich_without_TE'
    if te:
        return 'mixed_TE_satellite' if sat else 'partial_TE'
    if sat:
        return 'partial_satellite'
    return 'other_repeat' if sum(groups.values()) else 'unmasked'


def allele_sequences(row):
    result = set()
    for alternative in row['alt'].upper().split(','):
        ref, alt, start = row['ref'].upper(), alternative, int(row['pos']) - 1
        if not ref or not alt or (set(ref) | set(alt)) - set('ACGTN'):
            result.add(('nonliteral', start, ''))
            continue
        if ref == alt:
            continue
        right_ref, right_alt, left = len(ref), len(alt), 0
        while right_ref > 1 and right_alt > 1 and ref[right_ref-1] == alt[right_alt-1]:
            right_ref -= 1
            right_alt -= 1
        while right_ref-left > 1 and right_alt-left > 1 and ref[left] == alt[left]:
            left += 1
        ref, alt, start = ref[left:right_ref], alt[left:right_alt], start + left
        if len(ref) == 1 and len(alt) > 1 and ref[0] == alt[0]:
            result.add(('INS', start + 1, alt[1:]))
        elif len(alt) == 1 and len(ref) > 1 and ref[0] == alt[0]:
            result.add(('DEL', start + 1, ref[1:]))
        else:
            result.add(('complex_replacement', start, ''))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog-manifest', type=Path, default=Path('results/panel_sv_size_distribution.json'))
    parser.add_argument('--output-prefix', type=Path, default=Path('results/panel_repeat_composition'))
    args = parser.parse_args()
    csv.field_size_limit(sys.maxsize)
    catalog = json.loads(args.catalog_manifest.read_text())
    reference = pysam.FastaFile(str(UPSTREAM / 'resources/chm13v2.0.fa.gz'))
    rows, wanted = [], {}
    for source in catalog['sources']:
        seen = set()
        with Path(source['path']).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if float(row['ac']) < 2 or float(row['af']) > .99:
                    continue
                key = tuple(row[k] for k in ['pos', 'end', 'svtype', 'svlen'])
                if key in seen:
                    continue
                seen.add(key)
                result = {k: row[k] for k in ['chrom', 'pos', 'end', 'svtype', 'svlen', 'ac', 'af']}
                result['catalog_row_id'] = ':'.join(row[k] for k in ['chrom', 'pos', 'end', 'svtype', 'svlen'])
                result['_alleles'] = []
                if row['svtype'] in {'INS', 'DEL'}:
                    sequences = allele_sequences(row)
                    for kind, start, seq in sorted(sequences):
                        allele = dict(annotation_status='nonliteral_allele' if kind=='nonliteral' else 'complex_replacement' if kind=='complex_replacement' else 'pending', annotation_source='', seq_hash='', allele_bp=0, normalized_kind=kind, normalized_start=start)
                        result['_alleles'].append(allele)
                        if not seq or set(seq) - set('ACGTN'):
                            continue
                        allele.update(seq_hash=hashlib.sha1(seq.encode()).hexdigest()[:16], allele_bp=len(seq))
                        if kind == 'INS':
                            allele['annotation_status'] = 'missing_INS_annotation'
                            sha = hashlib.sha256(seq.encode()).hexdigest()
                            if allele['seq_hash'] in wanted and wanted[allele['seq_hash']] != sha:
                                raise ValueError('Truncated sequence-hash collision')
                            wanted[allele['seq_hash']] = sha
                        else:
                            actual = reference.fetch(row['chrom'], start, start + len(seq)).upper()
                            allele['annotation_status'] = ('reference_verified' if actual == seq else 'DEL_reference_mismatch')
                rows.append(result)
    reference.close()

    shards = sorted(p for p in (MEI / 'ins_arm').glob('shard_*') if p.is_dir())
    if len(shards) != 13:
        raise ValueError('Expected 13 completed insertion annotation shards')
    available, all_input = set(), set()
    for shard in shards:
        if not (shard / '.sentinel').exists():
            raise ValueError('Incomplete insertion RepeatMasker shard')
        name, pieces = None, []
        def finish():
            if name is None:
                return
            all_input.add(name)
            if name in wanted:
                seq = ''.join(pieces).upper()
                if hashlib.sha256(seq.encode()).hexdigest() != wanted[name]:
                    raise ValueError('Actual annotation input differs from catalog allele')
                available.add(name)
        with (shard / 'query.fa').open() as handle:
            for line in handle:
                if line.startswith('>'):
                    finish()
                    name, pieces = line[1:].split()[0], []
                elif name in wanted:
                    pieces.append(line.strip())
            finish()
    ins_hits = defaultdict(list)
    for shard in shards:
        for hit in parse_rm(shard / 'query.fa.out'):
            if hit['query'] in wanted:
                ins_hits[hit['query']].append(hit)

    ref_groups, ref_windows = defaultdict(lambda: defaultdict(list)), defaultdict(list)
    ref_shards = sorted(p for p in (MEI / 'chm13_hor_rm/shards').glob('shard_*') if p.is_dir())
    if len(ref_shards) != 41:
        raise ValueError('Expected 41 completed reference annotation shards')
    for shard in ref_shards:
        if not (shard / '.sentinel').exists():
            raise ValueError('Incomplete reference RepeatMasker shard')
        name = (shard / 'q.fa').open().readline().strip()[1:]
        chrom, interval = name.split(':')
        start, end = map(int, interval.split('-'))
        ref_windows[chrom].append((start, end))
        for hit in parse_rm(shard / 'q.fa.out'):
            if hit['query'] != name:
                raise ValueError('Reference query coordinate identity mismatch')
            ref_groups[chrom][hit['group']].append((start + hit['start'], start + hit['end']))
    indexes = {}
    for chrom, groups in ref_groups.items():
        indexes[chrom] = {}
        for group, intervals in groups.items():
            intervals = merge(intervals)
            indexes[chrom][group] = (intervals, [v[1] for v in intervals])
    mei = list(csv.DictReader((MEI / 'loci/mei_by_hash.tsv').open(), delimiter='\t'))
    mei_hashes = {m['seq_hash'] for m in mei}
    totals, by_type, status_counts, allele_totals = Counter(), defaultdict(Counter), Counter(), Counter()
    allele_rows = []
    for row in rows:
        alleles = row.pop('_alleles')
        for allele in alleles:
            groups, status = {}, allele['annotation_status']
            if allele['normalized_kind'] == 'INS' and allele['seq_hash'] in available:
                status = 'annotated'
                groups = coverage(ins_hits[allele['seq_hash']], allele['allele_bp'])['group_bp']
                allele['annotation_source'] = 'exact_INS_sequence_RepeatMasker'
            elif status == 'reference_verified':
                start, end = allele['normalized_start'], allele['normalized_start'] + allele['allele_bp']
                covered = union_bp([(max(start, a), min(end, b)) for a, b in ref_windows[row['chrom']] if a < end and b > start])
                if covered == allele['allele_bp']:
                    status = 'annotated'
                    groups = {group: overlap(iv, ends, start, end) for group, (iv, ends) in indexes.get(row['chrom'], {}).items()}
                    allele['annotation_source'] = 'CHM13_array_RepeatMasker_overlap_exact_DEL_DNA_verified'
                else:
                    status = 'DEL_outside_annotation_scope'
            allele['annotation_status'] = status
            allele['composition'] = label(allele['allele_bp'], groups, status)
            allele['TE_bp'] = groups.get('transposable_element', 0)
            allele['alpha_bp'] = groups.get('alpha_satellite', 0)
            allele['TE_fraction'] = allele['TE_bp'] / allele['allele_bp'] if allele['allele_bp'] else 0
            allele['alpha_fraction'] = allele['alpha_bp'] / allele['allele_bp'] if allele['allele_bp'] else 0
            allele['MEI_predicate_positive'] = allele['normalized_kind']=='INS' and allele['seq_hash'] in mei_hashes
            allele_totals[allele['composition']] += 1
            allele_rows.append(dict(catalog_row_id=row['catalog_row_id'], svtype=row['svtype'], **allele))
        row['n_distinct_literal_alleles'] = len(alleles)
        row['n_annotated_alleles'] = sum(a['annotation_status']=='annotated' for a in alleles)
        row['n_TE_rich_alleles'] = sum(a['TE_fraction']>=.8 for a in alleles)
        row['n_alpha_rich_alleles'] = sum(a['alpha_fraction']>=.8 for a in alleles)
        row['n_TE_positive_alleles'] = sum(a['TE_bp']>0 for a in alleles)
        row['n_inserted_TE_rich_alleles'] = sum(a['normalized_kind']=='INS' and a['TE_fraction']>=.8 for a in alleles)
        row['n_deleted_TE_rich_alleles'] = sum(a['normalized_kind']=='DEL' and a['TE_fraction']>=.8 for a in alleles)
        row['MEI_predicate_positive'] = any(a['MEI_predicate_positive'] for a in alleles)
        row['all_alleles_alpha_rich'] = bool(alleles) and row['n_alpha_rich_alleles']==len(alleles)
        labels = {a['composition'] for a in alleles}
        statuses = {a['annotation_status'] for a in alleles}
        row['annotation_status'] = ('unsupported_type' if not alleles else next(iter(statuses)) if len(statuses)==1 else 'partially_annotated')
        row['composition'] = ('unsupported_type' if not labels else next(iter(labels)) if len(labels)==1 else
                              'alpha_rich_with_TE' if row['all_alleles_alpha_rich'] else
                              'heterogeneous_cluster_TE_rich_present' if row['n_TE_rich_alleles'] else 'heterogeneous_cluster')
        totals[row['composition']] += 1
        by_type[row['svtype']][row['composition']] += 1
        status_counts[row['annotation_status']] += 1
    expected = next(r['n_variants'] for r in catalog['by_chromosome'] if r['chrom']=='ALL')
    if len(rows) != expected:
        raise ValueError('Repeat tally denominator differs from original panel size summary')
    report = dict(n_catalog_rows=len(rows), counts=dict(totals), by_type=dict(by_type), annotation_status=dict(status_counts),
        n_allele_occurrences=len(allele_rows), allele_composition_counts=dict(allele_totals),
        INS_annotation_input_unique_sequences=len(all_input), existing_MEI_positive_unique_sequences=len(mei),
        existing_MEI_positive_families=dict(Counter(r['mei_class'] for r in mei)),
        retained_INS_MEI_positive_catalog_rows=sum(r['MEI_predicate_positive'] for r in rows),
        n_TE_homology_catalog_rows=sum(r['n_TE_positive_alleles'] > 0 for r in rows),
        n_TE_rich_present_catalog_rows=sum(r['n_TE_rich_alleles'] > 0 for r in rows),
        n_all_alpha_rich_catalog_rows=sum(r['all_alleles_alpha_rich'] for r in rows),
        n_any_alpha_rich_catalog_rows=sum(r['n_alpha_rich_alleles'] > 0 for r in rows),
        n_TE_rich_and_any_alpha_rich_rows=sum(r['n_TE_rich_alleles'] > 0 and r['n_alpha_rich_alleles'] > 0 for r in rows),
        scope='Same retained catalog rows and filters as panel_sv_size_distribution.json; 23 available chromosomes',
        denominator='Catalog rows, not independent events, sequence hashes, carriers, or haplotypes; comma-separated alternatives are normalized against REF and annotated individually',
        thresholds='Rich means >=80% allele bases annotated; any TE homology is a separate, permissive flag',
        polarity='INS annotates added sample DNA; DEL annotates removed CHM13 DNA. Neither establishes transposition.',
        HOR='Alpha homology does not establish whole-HOR-unit gain/loss; validated HOR count is unknown',
        missing='Symbolic types and out-of-scope/mismatched sequences remain explicitly unclassified',
        method='Each literal REF/ALT pair is trimmed right then left to remove shared padding; pure INS/DEL DNA is extracted. Complex replacements remain unclassified. INS is verified against actual existing RM query DNA; normalized DEL DNA is verified against CHM13 then intersected with same-instrument 41-array annotation. Declared row SVTYPE does not force all alternatives to have that type.',
        provenance=dict(catalog_manifest=str(args.catalog_manifest), catalog_manifest_sha256=file_digest(args.catalog_manifest),
                        source_script_sha256=file_digest(Path(__file__)), annotation_root=str(MEI)))
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix, data in [('.tsv', rows), ('_alleles.tsv', allele_rows)]:
        with Path(str(args.output_prefix) + suffix).open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter='\t')
            writer.writeheader(); writer.writerows(data)
    args.output_prefix.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
