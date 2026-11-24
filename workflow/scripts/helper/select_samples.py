import csv
import os

def get_samples_by_group(csv_file, target_groups, max_samples=None):
    samples = {}
    
    group_versions = {
        "Group 1": "hprc_r2_v1.0.1",
        "Group 2": "hprc_r2_v1",
        "Group 3": "for_genbank", 
        "Group 4": "v1.0.genbank" 
    }

    if not os.path.exists(csv_file):
        print(f"Error: Classification file not found at {csv_file}")
        return {}

    count = 0
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if max_samples is not None and count >= max_samples:
                break
                
            sample = row['Sample']
            available_groups_str = row['Available_Groups']
            available_groups = [t.strip() for t in available_groups_str.split(',')]
            hap_type_group1_2 = row.get('Haplotype_Type_Group1_2', 'hap1_hap2')
            
            selected_group = None
            for group in target_groups:
                if group in available_groups:
                    selected_group = group
                    break
            
            if selected_group:
                samples[sample] = {
                    'group': selected_group,
                    'version': group_versions.get(selected_group, "unknown"),
                    'hap_type_group1_2': hap_type_group1_2
                }
                count += 1
                    
    return samples
