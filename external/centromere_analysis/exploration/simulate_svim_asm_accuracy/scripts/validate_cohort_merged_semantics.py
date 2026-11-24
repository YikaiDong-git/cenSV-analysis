#!/usr/bin/env python3
"""Fail closed if merging changes a symbolic or breakend carrier assertion."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys

import pysam


ROOT = Path(__file__).resolve().parents[1]


class SemanticMismatch(ValueError):
    pass


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _sample_name(index):
    return f'SIM{index // 2:03d}_hap{index % 2 + 1}'


def _special(alt):
    return alt == '*' or alt.startswith('<') or '[' in alt or ']' in alt


def _carrier(sample):
    gt = sample.get('GT')
    return bool(gt) and all(x is not None for x in gt) and any(x > 0 for x in gt)


def _selected_alt(record, sample_name):
    gt = record.samples[sample_name].get('GT')
    if not _carrier(record.samples[sample_name]):
        raise SemanticMismatch(f'Absent or missing source carrier: {record.id} / {sample_name}')
    indices = {x for x in gt if x > 0}
    if len(indices) != 1 or not record.alts or max(indices) > len(record.alts):
        raise SemanticMismatch(f'Ambiguous selected allele: {record.id} / {sample_name}')
    index = next(iter(indices))
    return index, record.alts[index - 1]


def _allele_value(record, field, allele_index, sample_name=None):
    if sample_name is None:
        if field not in record.info:
            return None
        value = record.info[field]
        number = record.header.info[field].number
    else:
        if field not in record.format:
            return None
        value = record.samples[sample_name].get(field)
        number = record.header.formats[field].number
    if value is None:
        return None
    if number in ('A', 'R'):
        values = value if isinstance(value, (tuple, list)) else (value,)
        offset = allele_index - 1 if number == 'A' else allele_index
        if offset >= len(values):
            raise SemanticMismatch(f'Missing selected {field}: {record.id}')
        return values[offset]
    if isinstance(value, (tuple, list)):
        if len(value) != 1:
            raise SemanticMismatch(f'Ambiguous non-allele-indexed {field}: {record.id}')
        return value[0]
    return value


def _signature(record, sample_name, allele_index, alt):
    info = str(record).split('\t')[7]
    explicit = [part[4:] for part in info.split(';') if part.startswith('END=')]
    if len(explicit) > 1:
        raise SemanticMismatch(f'Duplicated END: {record.id}')
    try:
        end = int(explicit[0]) if explicit else record.stop
    except ValueError as error:
        raise SemanticMismatch(f'Invalid END: {record.id}') from error
    return dict(chrom=record.chrom, pos=record.pos, end=end,
                svtype=_allele_value(record, 'SVTYPE', allele_index),
                svlen=_allele_value(record, 'SVLEN', allele_index),
                alt=alt, cutpaste=bool(record.info.get('CUTPASTE', False)),
                cn=_allele_value(record, 'CN', allele_index, sample_name))


def validate_chromosome(root, chrom, n_haplotypes, merged_path):
    root, merged_path = Path(root), Path(merged_path)
    sources, hashes = {}, {}
    expected_samples = {_sample_name(h) for h in range(n_haplotypes)}
    for h in range(n_haplotypes):
        path = root / 'chromosomes' / chrom / 'haps' / f'H{h:03d}' / 'merge_input.vcf.gz'
        sample = _sample_name(h)
        hashes[str(path)] = _sha256(path)
        with pysam.VariantFile(str(path)) as vcf:
            if list(vcf.header.samples) != [sample]:
                raise SemanticMismatch(f'Wrong source sample header: {path}')
            for record in vcf:
                identity = record.id
                if (not identity or not re.fullmatch(rf'H{h:03d}__[^;]+', identity)
                        or identity in sources or record.chrom != chrom):
                    raise SemanticMismatch(f'Invalid source identity or chromosome: {path} / {identity}')
                index, alt = _selected_alt(record, sample)
                special = _special(alt)
                sources[identity] = dict(sample=sample, special=special,
                    signature=_signature(record, sample, index, alt) if special else None)

    observed, counts = set(), Counter()
    with pysam.VariantFile(str(merged_path)) as vcf:
        if len(vcf.header.samples) != n_haplotypes or set(vcf.header.samples) != expected_samples:
            raise SemanticMismatch(f'Wrong merged sample header: {merged_path}')
        for record in vcf:
            counts['merged_rows'] += 1
            ids = record.id.split(';') if record.id else []
            if (record.chrom != chrom or not ids or len(set(ids)) != len(ids)
                    or observed.intersection(ids) or any(i not in sources for i in ids)):
                raise SemanticMismatch(f'Untraceable merged row: {record.chrom}:{record.pos} / {record.id}')
            observed.update(ids)
            members = [sources[i] for i in ids]
            samples = [member['sample'] for member in members]
            present = {s for s in vcf.header.samples if _carrier(record.samples[s])}
            if len(set(samples)) != len(samples) or set(samples) != present:
                raise SemanticMismatch(f'Merged carrier/source mismatch: {record.id}')
            for identity, member in zip(ids, members):
                sample = member['sample']
                index, alt = _selected_alt(record, sample)
                counts['source_carrier_assertions'] += 1
                if not member['special'] and not _special(alt):
                    counts['literal_assertions_delegated_to_existing_validator'] += 1
                    continue
                if not member['special'] or not _special(alt):
                    raise SemanticMismatch(f'Literal/symbolic representation changed: {identity} / {sample}')
                actual = _signature(record, sample, index, alt)
                expected = member['signature']
                differences = {key: {'source': expected[key], 'merged': actual[key]}
                               for key in expected if expected[key] != actual[key]}
                if differences:
                    raise SemanticMismatch(f'Symbolic/BND effect changed for {identity} / {sample}: '
                                           + json.dumps(differences, sort_keys=True))
                counts['checked_symbolic_or_bnd_assertions'] += 1
                if '[' in alt or ']' in alt:
                    counts['checked_bnd_assertions'] += 1
                if expected['svtype'] == 'DUP:INT':
                    counts['trace_augmented_dup_int_assertions'] += 1
    if observed != set(sources):
        raise SemanticMismatch(f'Missing source membership in {merged_path}: {len(set(sources) - observed)}')
    for key in ('merged_rows', 'source_carrier_assertions', 'checked_symbolic_or_bnd_assertions',
                'checked_bnd_assertions', 'trace_augmented_dup_int_assertions',
                'literal_assertions_delegated_to_existing_validator'):
        counts.setdefault(key, 0)
    return dict(passed=True, chrom=chrom, n_haplotypes=n_haplotypes, counts=dict(counts),
                merged_vcf=str(merged_path), merged_vcf_sha256=_sha256(merged_path),
                source_vcf_sha256=hashes,
                scope='Symbolic/BND field preservation, not truth validation; literal DNA equivalence is delegated.',
                dup_int_note='DUP:INT destination is supplied by the original candidate trace, not the serialized VCF.')


def validate_manifest(manifest_path, chrom=None, merged_path=None):
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    root = manifest_path.parent
    chromosomes = manifest['config']['chromosomes']
    if chrom is not None and chrom not in chromosomes:
        raise SemanticMismatch(f'Chromosome absent from cohort manifest: {chrom}')
    if merged_path is not None and chrom is None:
        raise SemanticMismatch('--merged-vcf requires --chrom')
    results = []
    for target in [chrom] if chrom is not None else chromosomes:
        path = Path(merged_path) if merged_path is not None else root / 'analysis' / target / 'merged.vcf.gz'
        if not path.exists() and merged_path is None:
            path = root / 'analysis' / target / 'merged.vcf'
        results.append(validate_chromosome(root, target, int(manifest['config']['n_haplotypes']), path))
    return dict(passed=True, manifest=str(manifest_path), manifest_sha256=_sha256(manifest_path),
                validator_sha256=_sha256(__file__), chromosomes=results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=ROOT / 'results/cohort_frequency_v1/manifest.json')
    parser.add_argument('--chrom')
    parser.add_argument('--merged-vcf', type=Path)
    parser.add_argument('--output', type=Path, help='Write a JSON success or failure receipt; no inputs are modified.')
    args = parser.parse_args()
    try:
        receipt = validate_manifest(args.manifest, args.chrom, args.merged_vcf)
    except (SemanticMismatch, OSError, ValueError) as error:
        receipt = dict(passed=False, error=str(error), manifest=str(args.manifest))
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(receipt, indent=2) + '\n')
        print(json.dumps(receipt), file=sys.stderr)
        return 1
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(passed=True, chromosomes=[r['chrom'] for r in receipt['chromosomes']],
                         counts=[r['counts'] for r in receipt['chromosomes']],
                         output=str(args.output) if args.output else None)))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
