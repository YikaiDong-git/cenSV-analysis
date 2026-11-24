#!/usr/bin/env python3
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import time

from classify_sv_call_errors_v2 import ALGORITHM_VERSION, classify_call_errors
from complete_truth_reference import CompositeReference
from evaluate_multichrom_sv import _write_json
from multichrom_truth_scoring import _ReferenceReadCache

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/multichrom_sv_v1'
OUTPUT = ROOT / 'results/multichrom_diagnostics_v2'
SCOPES = ('raw_full', 'production_active', 'target_active')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare():
    OUTPUT.mkdir(exist_ok=True)
    m = json.loads((SOURCE / 'manifest.json').read_text())
    m['diagnostic_source_run'] = str(SOURCE)
    m['diagnostic_source_manifest_sha256'] = digest(SOURCE / 'manifest.json')
    m['diagnostic_algorithm_version'] = ALGORITHM_VERSION
    m['diagnostic_classifier_sha256'] = digest(ROOT / 'scripts/classify_sv_call_errors_v2.py')
    if (OUTPUT / 'manifest.json').exists():
        assert json.loads((OUTPUT / 'manifest.json').read_text()) == m
    else:
        _write_json(OUTPUT / 'manifest.json', m)
    return m


def run_chromosome(index):
    m = json.loads((OUTPUT / 'manifest.json').read_text())
    sha = digest(ROOT / 'scripts/classify_sv_call_errors_v2.py')
    assert sha == m['diagnostic_classifier_sha256'], 'Classifier changed after preparation'
    chrom = m['config']['chromosomes'][index]
    for bundle in m['bundles']:
        if bundle['chrom'] != chrom:
            continue
        template = SOURCE / 'bundles' / bundle['id'] / 'template.fa.gz'
        for scenario in bundle['simulations']:
            source = SOURCE / 'bundles' / bundle['id'] / scenario
            dest = OUTPUT / 'bundles' / bundle['id'] / scenario
            if (dest / 'complete.json').exists():
                completed = json.loads((dest / 'complete.json').read_text())
                assert completed['classifier_sha256'] == sha
                continue
            started = time.monotonic()
            (dest / 'scoring').mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / 'truth.json', dest / 'truth.json')
            case = json.loads((source / 'truth.json').read_text())
            strict_hashes = {}
            for scope in SCOPES:
                strict_path = source / 'scoring' / (scope + '.json')
                strict = json.loads(strict_path.read_text())
                strict_hashes[str(strict_path)] = digest(strict_path)
                strict['metadata']['classifier_sha256'] = sha
                strict['metadata']['diagnostic_classifier_script'] = 'classify_sv_call_errors_v2.py'
                strict['metadata']['frozen_strict_source_sha256'] = strict_hashes[str(strict_path)]
                _write_json(dest / 'scoring' / (scope + '.json'), strict)
                name = f'parsed_calls_{scope}.json'
                shutil.copy2(source / 'scoring' / name, dest / 'scoring' / name)
            calls = json.loads((source / 'scoring/parsed_calls_production_active.json').read_text())
            strict = json.loads((source / 'scoring/production_active.json').read_text())
            reference = CompositeReference(template, m['reference'])
            print(f"Starting {case['id']}: {len(calls)} records", flush=True)
            try:
                result = classify_call_errors(_ReferenceReadCache(reference, max_bases=32_000_000),
                                              case['truth'], calls, strict)
            finally:
                reference.close()
            result['provenance'] = dict(classifier_sha256=sha,
                classifier_script='classify_sv_call_errors_v2.py', algorithm_version=ALGORITHM_VERSION,
                strictscorer_version=strict['conventions']['exact_search_version'],
                frozen_strict_sources=strict_hashes,
                source_parsed_calls_sha256=digest(source / 'scoring/parsed_calls_production_active.json'),
                source_truth_sha256=digest(source / 'truth.json'),
                source_diagnostics_sha256=digest(source / 'scoring/call_errors.json'))
            _write_json(dest / 'scoring/call_errors.json', result)
            old = json.loads((source / 'scoring/call_errors.json').read_text())
            transition = Counter((a['category'], b['category']) for a, b in zip(old['rows'], result['rows']))
            elapsed = time.monotonic() - started
            _write_json(dest / 'complete.json', dict(completed=True, case_id=case['id'],
                classifier_sha256=sha, elapsed_seconds=elapsed, records=len(calls),
                search_bounded_before=old['summary']['search_bounded_records'],
                search_bounded_after=result['summary']['search_bounded_records'],
                transitions=[dict(before=a, after=b, count=n) for (a, b), n in transition.items()]))
            print(f"Finished {case['id']}: {elapsed:.2f}s; {dict(transition)}", flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare', 'run'])
    p.add_argument('--chromosome-index', type=int)
    a = p.parse_args()
    if a.action == 'prepare':
        prepare()
    else:
        assert a.chromosome_index is not None and 0 <= a.chromosome_index < 23
        run_chromosome(a.chromosome_index)
