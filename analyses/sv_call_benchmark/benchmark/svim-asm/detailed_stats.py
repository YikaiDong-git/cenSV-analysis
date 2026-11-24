import csv
import collections

REF_FILE = "../svim_asm/mapq5_minsvlen50_maxpart10_partdist2000_maxedit1000_maxsv100000_qgap50_rgap50/chr17/chr17_variants.tsv"
TGT_FILE = "../svim_asm/mapq5_minsvlen10_maxpart30_partdist500_maxedit200_maxsv100000_qgap50_rgap50/chr17/chr17_variants.tsv"

def load_data(filepath):
    data = {}
    header = None
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            header = reader.fieldnames
            for row in reader:
                key = (row['chrom'], row['pos'], row['end'], row['svtype'])
                data[key] = row
    except FileNotFoundError:
        pass
    return data, header


def analyze():
    ref_data, header = load_data(REF_FILE)
    tgt_data, _ = load_data(TGT_FILE)

    ref_keys = set(ref_data.keys())
    tgt_keys = set(tgt_data.keys())
    all_keys = ref_keys | tgt_keys
    total_union = len(all_keys)

    missing_in_tgt = ref_keys - tgt_keys
    added_in_tgt = tgt_keys - ref_keys
    intersection = ref_keys & tgt_keys

    diff_col_counts = collections.defaultdict(int)
    modified_keys = []
    identical_keys = []

    compare_cols = [c for c in header if c != 'id'] if header else []

    for k in intersection:
        r = ref_data[k]
        t = tgt_data[k]
        is_diff = False
        diffs_here = []
        for col in compare_cols:
            if r[col] != t[col]:
                is_diff = True
                diffs_here.append(col)
        if is_diff:
            modified_keys.append(k)
            for d in diffs_here:
                diff_col_counts[d] += 1
        else:
            identical_keys.append(k)

    out = {
        'total_union': total_union,
        'identical': len(identical_keys),
        'added': len(added_in_tgt),
        'missing': len(missing_in_tgt),
        'modified': len(modified_keys),
        'diff_by_column': dict(diff_col_counts)
    }

    for k, v in out.items():
        print(f"{k}: {v}")

if __name__ == '__main__':
    analyze()
