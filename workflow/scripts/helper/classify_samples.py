import re
import csv
import os
import subprocess
from collections import defaultdict
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

INPUT_FILE = CENSV_ROOT + '/workflow/scripts/25.12.8_explore_assemblies.txt'
OUTPUT_FILE = CENSV_ROOT + '/workflow/scripts/sample_classification.csv'

def classify_sample(sample, files):
    groups = []
    notes = []
    
    if any(f.endswith('hprc_r2_v1.0.1.fa.gz') for f in files): groups.append("Group 1")
    if any(f.endswith('hprc_r2_v1.fa.gz') for f in files):     groups.append("Group 2")
    has_genbank = any('for_genbank.fa.gz' in f for f in files)
    has_polished = any('.polished.fasta' in f for f in files)
    if has_genbank or has_polished:
        groups.append("Group 3")
        if has_polished: notes.append("Polished")
        if has_genbank: notes.append("GenBank")

    if sample == 'HG06807':
        groups.append("Group 4")
        notes.append("RepeatMasker masked")

    std_files = [f for f in files if 'hprc_r2_v1' in f]
    hap_type = "mat_pat" if any('_mat_' in f or '_pat_' in f for f in std_files) else "hap1_hap2"

    return ", ".join(groups), hap_type, "; ".join(notes)

def fetch_s3_dataset():
    if os.path.exists(INPUT_FILE):
        print(f"{INPUT_FILE} already exists. Skipping fetch.")
        return True

    cmd = [
        "aws", "s3", "ls",
        "s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/",
        "--no-sign-request", "--recursive"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        with open(INPUT_FILE, "w") as f:
            f.write(result.stdout)
        return True
    return False

def main():
    fetch_s3_dataset()
    
    print(f"Reading {INPUT_FILE}...")
    sample_files = defaultdict(set)
    pattern = re.compile(r'submissions/[^/]+/([^/]+)/assemblies/')

    with open(INPUT_FILE, 'r') as f:
        for line in f:
            if match := pattern.search(line):
                sample_files[match.group(1)].add(line.strip().split('/')[-1])

    results = sorted([
        [sample, *classify_sample(sample, files)] 
        for sample, files in sample_files.items()
    ], key=lambda x: x[0])

    print(f"Writing to {OUTPUT_FILE}...")
    with open(OUTPUT_FILE, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Sample', 'Available_Groups', 'Haplotype_Type_Group1_2', 'Notes'])
        writer.writerows(results)
    print("Done.")

if __name__ == "__main__":
    main()
