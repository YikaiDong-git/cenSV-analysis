import os
import csv

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'svim_asm'))
TEMPLATE = 'mapq{mapq}_minsvlen50_maxpart10_partdist2000_maxedit1000_maxsv100000_qgap50_rgap50'
CHROM = 'chr6'
MAPQS = [3,5,10,20]


def load(folder, chrom):
    path = os.path.join(BASE, folder, chrom, f"{chrom}_variants.tsv")
    data = {}
    if not os.path.exists(path):
        return data
    with open(path,'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for r in reader:
            k = (r['chrom'], r['pos'], r['end'], r['svtype'])
            data[k] = r
    return data

ref = load(TEMPLATE.format(mapq=5), CHROM)
for m in [3,10,20]:
    tgt = load(TEMPLATE.format(mapq=m), CHROM)
    keys = set(ref.keys()) | set(tgt.keys())
    added = len(set(tgt.keys()) - set(ref.keys()))
    missing = len(set(ref.keys()) - set(tgt.keys()))
    modified = 0
    for k in set(ref.keys()) & set(tgt.keys()):
        r = ref[k]; t = tgt[k]
        diffs = [c for c in ['ac','svlen'] if r.get(c) != t.get(c)]
        if diffs:
            modified += 1
    total = len(keys)
    print(f"mapq{m} vs mapq5: total={total}, added={added} ({added/total:.2f}), missing={missing} ({missing/total:.2f}), modified={modified} ({modified/total:.2f})")
