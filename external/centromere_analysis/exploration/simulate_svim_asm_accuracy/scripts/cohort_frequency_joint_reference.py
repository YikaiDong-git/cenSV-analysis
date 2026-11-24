"""Known four-chromosome haplotypes for complete-reference cohort evaluation."""
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from cohort_frequency_rescore import SNVOverlayReference
from evaluate_multichrom_sv import _write_json


def verify_joint_code_identity(root):
    path = Path(root) / 'joint_execution_code_identity.json'
    identity = json.loads(path.read_text())
    for source, expected in identity.items():
        if hashlib.sha256(Path(source).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Joint execution code changed after freezing: {source}')
    return hashlib.sha256(path.read_bytes()).hexdigest()


class JointSNVOverlayReference:
    def __init__(self, root, manifest, haplotype_index):
        self.root = Path(root)
        self.h = haplotype_index
        self.references = {}
        self.truth_records = {}
        for chrom in manifest['config']['chromosomes']:
            folder = self.root / 'chromosomes' / chrom / 'haps' / f'H{self.h:03d}'
            truth = json.loads((folder / 'truth.json').read_text())
            self.truth_records[chrom] = truth
            founder = self.root / 'chromosomes' / chrom / 'founders' / f"F{truth['founder_index']}.fa.gz"
            self.references[chrom] = SNVOverlayReference(founder, manifest['reference'],
                                                        folder / 'private_snvs.tsv.gz', chrom)
        self.fallback = next(iter(self.references.values()))

    def fetch(self, chrom, start=None, end=None):
        return self.references.get(chrom, self.fallback).fetch(chrom, start, end)

    def get_reference_length(self, chrom):
        return self.references.get(chrom, self.fallback).get_reference_length(chrom)

    def verify_templates(self):
        result = {}
        for chrom, expected in self.truth_records.items():
            sequence = self.fetch(chrom)
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            if digest != expected['template_sha256']:
                raise ValueError(f'Joint background identity mismatch: {chrom}/H{self.h:03d}')
            result[chrom] = digest
        return result

    def close(self):
        for reference in self.references.values():
            reference.close()


def joint_truth(root, manifest, h):
    from cohort_frequency_metrics import sample_name
    root = Path(root)
    events, windows, source_hashes, templates, queries = [], [], {}, {}, {}
    for chrom in manifest['config']['chromosomes']:
        folder = root / 'chromosomes' / chrom / 'haps' / f'H{h:03d}'
        path = folder / 'truth.json'
        record = json.loads(path.read_text())
        if not json.loads((folder / 'complete.json').read_text()).get('completed'):
            raise ValueError(f'Incomplete source case: {folder}')
        if not json.loads((folder / 'oracle.json').read_text()).get('passed'):
            raise ValueError(f'Source oracle failed: {folder}')
        if record['chrom'] != chrom or record['haplotype_index'] != h or record['sample_name'] != sample_name(h):
            raise ValueError(f'Source haplotype identity differs: {folder}')
        for e in record['truth']:
            events.append(dict(e, id=f"{chrom}__{e['id']}", original_id=e['id']))
        for window in record['negative_windows']:
            windows.append(dict(window, id=f"{chrom}__{window['id']}"))
        source_hashes[chrom] = hashlib.sha256(path.read_bytes()).hexdigest()
        templates[chrom] = record['template_sha256']
        queries[chrom] = record['query_sha256']
    return dict(id=f'joint_H{h:03d}', haplotype_index=h, sample_name=sample_name(h), truth=events,
                negative_windows=windows, source_truth_sha256=source_hashes,
                template_sha256_by_chromosome=templates, query_sha256_by_chromosome=queries,
                query_scope='Four chromosome query contigs; full CHM13 reference evaluation')


def _rescore_one_joint_haplotype(args):
    from score_cohort_joint_truth import score_joint_truth
    root, joint_root, manifest, h, selections = args
    root, joint_root = Path(root), Path(joint_root)
    identity_sha256 = verify_joint_code_identity(root)
    folder = joint_root / 'haps' / f'H{h:03d}'
    truth = json.loads((folder / 'truth.json').read_text())
    parsed = json.loads((folder / 'scoring/parsed_calls_production_active.json').read_text())
    reference = JointSNVOverlayReference(root, manifest, h)
    try:
        identities = reference.verify_templates()
        cache, result = {}, {}
        for percent in (0, 1, 5, 10):
            kept = selections[percent]
            calls = [c for c in parsed if f'H{h:03d}__' + c['id'] in kept]
            key = tuple(c['id'] for c in calls)
            if key not in cache:
                cache[key] = score_joint_truth(reference, truth['truth'], calls,
                                              truth['negative_windows'], **manifest['strict_scoring_parameters'])
            score = dict(cache[key], metadata=dict(cache[key].get('metadata', {}),
                joint_execution_code_identity_sha256=identity_sha256, haplotype_index=h,
                maf_filter_percent=percent, template_sha256_by_chromosome=identities))
            _write_json(folder / f'scoring/filtered_maf_{percent}.json', score)
            result[percent] = score
    finally:
        reference.close()
    return h, identities, result


def rescore_joint_thresholds(root, joint_root, manifest, rows, memberships, sources):
    root, joint_root = Path(root), Path(joint_root)
    n = manifest['config']['n_haplotypes']
    by_hap = {h: {p: set() for p in (0, 1, 5, 10)} for h in range(n)}
    for member in memberships:
        row = rows[member['row_index']]
        if not row['dedup_retained']:
            continue
        for percent in (0, 1, 5, 10):
            if row['maf'] >= percent / 100:
                by_hap[member['haplotype_index']][percent].add(member['source_id'])
    result = {p: {} for p in (0, 1, 5, 10)}
    seen = {chrom: set() for chrom in manifest['config']['chromosomes']}
    tasks = [(str(root), str(joint_root), manifest, h, by_hap[h]) for h in range(n)]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for h, identities, scores in pool.map(_rescore_one_joint_haplotype, tasks):
            for chrom, digest in identities.items():
                if digest in seen[chrom]:
                    raise ValueError(f'Duplicate joint chromosome background: {chrom}/H{h:03d}')
                seen[chrom].add(digest)
            for percent in (0, 1, 5, 10):
                result[percent][h] = scores[percent]
            if h % 16 == 0:
                print(f'Joint filtered rescore H{h:03d}', flush=True)
    return result
