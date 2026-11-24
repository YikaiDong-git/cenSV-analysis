"""Re-evaluate exact structural truth after cohort frequency filtering."""
from bisect import bisect_left
import csv
import gzip
import hashlib
import json
from pathlib import Path

from complete_truth_reference import CompositeReference
from evaluate_multichrom_sv import _write_json
from multichrom_truth_scoring import score_complete_truth


class SNVOverlayReference:
    def __init__(self, founder, reference, ledger, chrom):
        self.base=CompositeReference(founder,reference);self.chrom=chrom
        self.positions=[];self.old=[];self.new=[]
        with gzip.open(ledger,'rt') as f:
            for row in csv.DictReader(f,delimiter='\t'):
                pos=int(row['position_0based'])
                if row['chrom']!=chrom or (self.positions and pos<=self.positions[-1]):
                    raise ValueError('Private SNP ledger is unsorted, repeated or wrong chromosome')
                self.positions.append(pos);self.old.append(ord(row['ref']));self.new.append(ord(row['alt']))

    def get_reference_length(self,chrom):
        return self.base.get_reference_length(chrom)

    def fetch(self,chrom,start=None,end=None):
        start=0 if start is None else start
        end=self.get_reference_length(chrom) if end is None else end
        text=self.base.fetch(chrom,start,end)
        if chrom!=self.chrom:return text
        left=bisect_left(self.positions,start);right=bisect_left(self.positions,end)
        if left==right:return text
        sequence=bytearray(text.encode())
        for i in range(left,right):
            relative=self.positions[i]-start
            if sequence[relative]!=self.old[i]:
                raise ValueError('Private SNP reference base disagrees with founder')
            sequence[relative]=self.new[i]
        return sequence.decode()

    def close(self):
        self.base.close()


def rescore_thresholds(root,chrom,manifest,rows,memberships,sources):
    root=Path(root);n=manifest['config']['n_haplotypes']
    thresholds={0:0.0,1:.01,5:.05,10:.10}
    by_haplotype={h:{p:set() for p in thresholds} for h in range(n)}
    for member in memberships:
        row=rows[member['row_index']];h=member['haplotype_index']
        if not row['dedup_retained']:continue
        for percent,threshold in thresholds.items():
            if row['maf']>=threshold:by_haplotype[h][percent].add(member['source_id'])
    result={p:{} for p in thresholds}
    seen_templates=set()
    for h in range(n):
        folder=root/'chromosomes'/chrom/'haps'/f'H{h:03d}'
        truth=json.loads((folder/'truth.json').read_text())
        if truth['template_sha256'] in seen_templates:
            raise ValueError(f'Duplicated chromosome background in cohort: {chrom} H{h:03d}')
        seen_templates.add(truth['template_sha256'])
        founder=root/'chromosomes'/chrom/'founders'/f"F{truth['founder_index']}.fa.gz"
        ledger=folder/'private_snvs.tsv.gz'
        reference=SNVOverlayReference(founder,manifest['reference'],ledger,chrom)
        try:
            rebuilt=reference.fetch(chrom)
            if hashlib.sha256(rebuilt.encode()).hexdigest()!=truth['template_sha256']:
                raise ValueError(f'Reconstructed background identity differs: {chrom} H{h:03d}')
            del rebuilt
            parsed=json.loads((folder/'scoring/parsed_calls_production_active.json').read_text())
            cache={}
            for percent in thresholds:
                retained=by_haplotype[h][percent]
                calls=[c for c in parsed if f'H{h:03d}__'+c['id'] in retained]
                key=tuple(c['id'] for c in calls)
                if key not in cache:
                    cache[key]=score_complete_truth(reference,truth['truth'],calls,truth['negative_windows'],
                                                    **manifest['strict_scoring_parameters'])
                score=cache[key]
                _write_json(folder/f'scoring/filtered_maf_{percent}.json',score)
                result[percent][h]=score
        finally:reference.close()
        if h%32==0:print(f'Rescored {chrom} H{h:03d}',flush=True)
    return result
