import pandas as pd
import argparse

def main(my_samples_file, g1k_samples_file, output_file, target_groups=None):
    my_df = pd.read_csv(my_samples_file)
    if target_groups:
        print(f"Filtering for groups: {target_groups}")
        def has_overlap(group_str):
            if pd.isna(group_str): return False
            available = [g.strip() for g in str(group_str).split(',')]
            return any(t in available for t in target_groups)
        my_df = my_df[my_df['Available_Groups'].apply(has_overlap)]
    my_samples = set(my_df['Sample'].astype(str).str.strip().tolist())

    g1k_df = pd.read_csv(g1k_samples_file, sep=r'\s+')
    g1k_samples = set(g1k_df['SampleID'].astype(str).str.strip().tolist())
    overlap = my_samples.intersection(g1k_samples)
    print(f"Found {len(overlap)} overlapping samples.")
    with open(output_file, 'w') as f:
        for sample in sorted(overlap):
            f.write(f"{sample}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--my_samples", required=True, help="Path to local sample classification CSV")
    parser.add_argument("--g1k_samples", required=True, help="Path to 1000G samples PED file")
    parser.add_argument("--output", required=True, help="Path to output text file")
    parser.add_argument("--groups", help="Comma-separated list of target groups to filter by")
    args = parser.parse_args()
    groups_list = [g.strip() for g in args.groups.split(',')] if args.groups else None
    
    main(args.my_samples, args.g1k_samples, args.output, groups_list)
