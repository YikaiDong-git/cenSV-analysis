import os
import csv
import re

import argparse

PARAM_REGEX = re.compile(
    r"mapq(\d+)_minsvlen(\d+)_maxpart(\d+)_partdist(\d+)_maxedit(\d+)_maxsv(\d+)_qgap(\d+)_rgap(\d+)"
)

def parse_params(folder_name):
    m = PARAM_REGEX.match(folder_name)
    if m:
        return {
            'mapq': int(m.group(1)),
            'minsvlen': int(m.group(2)),
            'raw': folder_name
        }
    return None

def load(path):
    data = {}
    if not os.path.exists(path):
        return data
    with open(path,'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for r in reader:
            k = (r['chrom'], r['pos'], r['end'], r['svtype'])
            data[k] = r
    return data

def main():
    parser = argparse.ArgumentParser(description='Global analysis of SVIM-ASM results')
    parser.add_argument('--input-dir', required=True, help='Directory containing results')
    args = parser.parse_args()
    
    BASE_DIR = args.input_dir
    
    if not os.path.exists(BASE_DIR):
        print("Base directory not found")
        return

    folders = [d for d in os.listdir(BASE_DIR) if os.path.isdir(os.path.join(BASE_DIR, d))]
    parsed = []
    for f in folders:
        p = parse_params(f)
        if p:
            parsed.append(p)

    if not parsed:
        print("No valid parameter folders found.")
        return

    ref_run = next((p for p in parsed if p['mapq']==5 and p['minsvlen']==50), None)
    if not ref_run:
        ref_run = parsed[0]
        print(f"Default reference (Q5 L50) not found, using {ref_run['raw']}")
    else:
        print(f"Using Reference: {ref_run['raw']}")

    tgt_run = next((p for p in parsed if p['minsvlen']==10 and p['mapq']==ref_run['mapq']), None)
    if not tgt_run:
        tgt_run = next((p for p in parsed if p['raw'] != ref_run['raw']), None)
    
    if not tgt_run:
        print("No differing runs found to compare.")
        return

    print(f"Comparing Reference vs Target: {tgt_run['raw']}")
    print("-" * 80)
    
    ref_folder = os.path.join(BASE_DIR, ref_run['raw'])
    chroms = [d for d in os.listdir(ref_folder) if d.startswith('chr') and os.path.isdir(os.path.join(ref_folder, d))]
    
    for chrom in sorted(chroms):
        refp = os.path.join(ref_folder, chrom, f"{chrom}_variants.tsv")
        tgtp = os.path.join(BASE_DIR, tgt_run['raw'], chrom, f"{chrom}_variants.tsv")
        
        refd = load(refp)
        tgtd = load(tgtp)
        
        keys = set(refd.keys()) | set(tgtd.keys())
        added = len(set(tgtd.keys()) - set(refd.keys()))
        missing = len(set(refd.keys()) - set(tgtd.keys()))
        modified = 0
        
        common = set(refd.keys()) & set(tgtd.keys())
        for k in common:
            r = refd[k]; t = tgtd[k]
            for col in ['ac','svlen']:
                if r.get(col) != t.get(col):
                    modified += 1
                    break
                    
        total = len(keys)
        pct = (added+missing+modified)/total*100 if total else 0
        print(f"{chrom:<6}: Total={total:<6} Diff={pct:>6.2f}% (Add={added:<4} Miss={missing:<4} Mod={modified:<4})")

if __name__ == '__main__':
    main()
