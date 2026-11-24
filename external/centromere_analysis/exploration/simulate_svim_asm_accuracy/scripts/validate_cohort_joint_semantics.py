#!/usr/bin/env python3
"""Validate every source carrier and allele in a joint cohort merge."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import sys

import pysam

from vcf_literal_equivalence import exact_padding_equivalence
from validate_cohort_merged_semantics import (SemanticMismatch, _carrier, _sample_name,
    _selected_alt, _sha256, _signature, _special)


def validate_joint(manifest_path, merged_path=None):
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    root = manifest_path.parent
    n = int(manifest['config']['n_haplotypes'])
    merged = Path(merged_path) if merged_path else root / 'joint_analysis/merged.vcf.gz'
    sources, hashes = {}, {}
    for h in range(n):
        path = root / 'joint_cohort/haps' / f'H{h:03d}' / 'merge_input.vcf.gz'
        sample = _sample_name(h)
        hashes[str(path.resolve())] = _sha256(path)
        with pysam.VariantFile(str(path)) as vcf:
            if list(vcf.header.samples) != [sample]:
                raise SemanticMismatch(f'Incorrect sample header: {path}')
            for rec in vcf:
                identity = rec.id
                if not identity or not re.fullmatch(rf'H{h:03d}__[^;]+', identity) or identity in sources:
                    raise SemanticMismatch(f'Invalid or duplicate source ID: {identity}')
                index, alt = _selected_alt(rec, sample)
                sources[identity] = dict(sample=sample, special=_special(alt),
                    signature=_signature(rec, sample, index, alt) if _special(alt) else None,
                    chrom=rec.chrom, start=rec.start, ref=rec.ref, alt=alt)
    observed, counts, chromosomes = set(), Counter(), set()
    with pysam.VariantFile(str(merged)) as vcf:
        expected_samples = {_sample_name(h) for h in range(n)}
        if len(vcf.header.samples) != n or set(vcf.header.samples) != expected_samples:
            raise SemanticMismatch('Incorrect merged sample header')
        for rec in vcf:
            chromosomes.add(rec.chrom); counts['merged_rows'] += 1
            ids = rec.id.split(';') if rec.id else []
            if not ids or len(ids) != len(set(ids)) or observed.intersection(ids) or any(i not in sources for i in ids):
                raise SemanticMismatch(f'Untraceable joint row {rec.chrom}:{rec.pos} / {rec.id}')
            observed.update(ids)
            members = [sources[i] for i in ids]
            samples = [m['sample'] for m in members]
            present = {s for s in vcf.header.samples if _carrier(rec.samples[s])}
            if len(samples) != len(set(samples)) or set(samples) != present:
                raise SemanticMismatch(f'Joint row carrier/source collision: {rec.id}')
            for identity, member in zip(ids, members):
                index, alt = _selected_alt(rec, member['sample'])
                if _special(alt) or member['special']:
                    if not _special(alt) or not member['special']:
                        raise SemanticMismatch(f'Literal/symbolic mismatch: {identity}')
                    actual = _signature(rec, member['sample'], index, alt)
                    if actual != member['signature']:
                        raise SemanticMismatch(f'Symbolic/BND effect changed: {identity}: {member["signature"]} -> {actual}')
                    counts['symbolic_or_bnd_assertions'] += 1
                    if '[' in alt or ']' in alt:
                        counts['bnd_assertions'] += 1
                else:
                    if rec.chrom != member['chrom'] or not exact_padding_equivalence(
                            rec.start, rec.ref, alt, member['start'], member['ref'], member['alt']):
                        raise SemanticMismatch(f'Literal allele changed: {identity}')
                    counts['literal_assertions'] += 1
                counts['source_carrier_assertions'] += 1
    if observed != set(sources):
        raise SemanticMismatch(f'Joint merge lost {len(set(sources)-observed)} source records')
    for key in ('merged_rows', 'source_carrier_assertions', 'literal_assertions', 'symbolic_or_bnd_assertions', 'bnd_assertions'):
        counts.setdefault(key, 0)
    return dict(passed=True, manifest_sha256=_sha256(manifest_path), validator_sha256=_sha256(__file__),
                n_haplotypes=n, truth_chromosomes=manifest['config']['chromosomes'],
                reported_chromosomes=sorted(chromosomes), merged_vcf=str(merged.resolve()),
                merged_vcf_sha256=_sha256(merged), source_vcf_sha256=hashes, counts=dict(counts),
                scope='Literal and symbolic/BND merged-allele identity; source-to-query linkage checked by joint run audit.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--merged-vcf')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output); output.parent.mkdir(parents=True, exist_ok=True)
    try:
        receipt = validate_joint(args.manifest, args.merged_vcf)
    except (OSError, ValueError) as error:
        output.write_text(json.dumps(dict(passed=False, error=str(error)), indent=2) + '\n')
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
    output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(passed=True, counts=receipt['counts'], output=str(output))))


if __name__ == '__main__':
    main()
