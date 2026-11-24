#!/usr/bin/env python3
"""Re-score completed early cohort calls at the frozen final search capacity."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

from classify_sv_call_errors_v2 import classify_call_errors
from cohort_frequency_rescore import SNVOverlayReference
from evaluate_multichrom_sv import _write_json
from multichrom_truth_scoring import score_complete_truth


def finalize_case(root, manifest_path, receipt, receipt_path, folder, marker, truth, old, new, archive, params):
    changed = any(old[k] != new[k] for k in ('metrics', 'calls', 'events', 'components'))
    done = json.loads(marker.read_text())
    done['metrics'] = new['metrics']
    done['canonical_raw_scoring_parameters'] = params
    done['original_completion_sha256'] = hashlib.sha256((archive / 'complete.json').read_bytes()).hexdigest()
    _write_json(marker, done)
    receipt['cases'][truth['id']] = dict(chrom=truth['chrom'], haplotype=folder.name, old_metrics=old['metrics'],
        new_metrics=new['metrics'], outcomes_changed=changed, archive=str(archive.relative_to(root)))
    receipt['manifest_sha256'] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    _write_json(receipt_path, receipt)


def canonicalize(manifest_path):
    manifest_path = Path(manifest_path).resolve()
    root = manifest_path.parent
    manifest = json.loads(manifest_path.read_text())
    params = manifest['strict_scoring_parameters']
    receipt_path = root / 'raw_capacity_canonicalization.json'
    receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else {'cases': {}}
    for chrom in manifest['config']['chromosomes']:
        for marker in sorted((root / 'chromosomes' / chrom / 'haps').glob('H*/complete.json')):
            folder = marker.parent
            path = folder / 'scoring/production_active.json'
            old = json.loads(path.read_text())
            archive = folder / 'scoring_before_capacity_canonicalization'
            if all(old['conventions'].get(k) == v for k, v in params.items()):
                if old.get('metadata', {}).get('capacity_refinement') == 'canonical_raw_rescore_8_to_12':
                    truth = json.loads((folder / 'truth.json').read_text())
                    done = json.loads(marker.read_text())
                    if (truth['id'] in receipt['cases'] and done.get('canonical_raw_scoring_parameters') == params
                            and done.get('metrics') == old['metrics'] and done.get('original_completion_sha256')):
                        continue
                    original = json.loads((archive / 'production_active.json').read_text())
                    finalize_case(root, manifest_path, receipt, receipt_path, folder, marker, truth, original, old, archive, params)
                continue
            differences = {k for k, v in params.items() if old['conventions'].get(k) != v}
            if differences != {'max_component_truth'} or old['conventions']['max_component_truth'] != 8 or params['max_component_truth'] != 12:
                raise ValueError(f'Unexpected scoring discrepancy: {folder}: {differences}')
            archive.mkdir(exist_ok=True)
            for name in ('production_active.json', 'call_errors.json'):
                if not (archive / name).exists():
                    shutil.copy2(folder / 'scoring' / name, archive / name)
            if not (archive / 'complete.json').exists():
                shutil.copy2(marker, archive / 'complete.json')
            truth = json.loads((folder / 'truth.json').read_text())
            calls = json.loads((folder / 'scoring/parsed_calls_production_active.json').read_text())
            founder = root / 'chromosomes' / chrom / 'founders' / f"F{truth['founder_index']}.fa.gz"
            ref = SNVOverlayReference(founder, manifest['reference'], folder / 'private_snvs.tsv.gz', chrom)
            try:
                new = score_complete_truth(ref, truth['truth'], calls, truth['negative_windows'], **params)
                new['metadata'] = dict(old.get('metadata', {}), capacity_refinement='canonical_raw_rescore_8_to_12',
                                       original_score_sha256=hashlib.sha256((archive / 'production_active.json').read_bytes()).hexdigest())
                changed = any(old[k] != new[k] for k in ('metrics', 'calls', 'events', 'components'))
                if changed:
                    _write_json(folder / 'scoring/call_errors.json', classify_call_errors(ref, truth['truth'], calls, new))
            finally:
                ref.close()
            finalize_case(root, manifest_path, receipt, receipt_path, folder, marker, truth, old, new, archive, params)
            _write_json(path, new)
            print(f"Canonicalized {truth['id']}; outcomes_changed={changed}", flush=True)
    print(f"Canonicalized total: {len(receipt['cases'])}", flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', default='results/cohort_frequency_v1/manifest.json')
    canonicalize(parser.parse_args().manifest)
