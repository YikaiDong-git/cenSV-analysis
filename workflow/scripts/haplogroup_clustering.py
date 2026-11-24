#!/usr/bin/env python3
"""Centromere haplogroup clustering with phylogenetic SV decomposition."""
import argparse
import logging
import sys
import re
from pathlib import Path

import cyvcf2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2
from sklearn.metrics import silhouette_score

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

METHODS = ['jaccard_nj', 'jaccard_ward', 'cosine_nj', 'cosine_ward', 'correlation_nj', 'correlation_ward']

def _compute_distance_and_linkage(matrix, samples, method):
    nj_tree = None

    if method == 'jaccard_nj':
        dist_vector = pdist(matrix.astype(bool), metric='jaccard')
        dist_vector = np.nan_to_num(dist_vector, nan=1.0)
        dist_matrix = squareform(dist_vector)
        linkage_matrix = linkage(dist_vector, method='average')

        from skbio import DistanceMatrix as SkbioDistanceMatrix
        from skbio.tree import nj
        dm = SkbioDistanceMatrix(dist_matrix, ids=samples)
        nj_tree = nj(dm)
        logger.info(f"    NJ tree built with {len(samples)} tips")

    elif method == 'jaccard_ward':
        dist_vector = pdist(matrix.astype(bool), metric='jaccard')
        dist_vector = np.nan_to_num(dist_vector, nan=1.0)
        dist_matrix = squareform(dist_vector)
        linkage_matrix = linkage(dist_vector, method='ward')
        logger.info(f"    Jaccard+Ward linkage built with {len(samples)} samples")

    elif method == 'cosine_nj':
        dist_vector = np.nan_to_num(pdist(matrix, metric='cosine'), nan=1.0)
        dist_matrix = squareform(dist_vector)
        linkage_matrix = linkage(dist_vector, method='average')

        from skbio import DistanceMatrix as SkbioDistanceMatrix
        from skbio.tree import nj
        dm = SkbioDistanceMatrix(dist_matrix, ids=samples)
        nj_tree = nj(dm)
        logger.info(f"    Cosine NJ tree built with {len(samples)} tips")

    elif method == 'cosine_ward':
        dist_vector = np.nan_to_num(pdist(matrix, metric='cosine'), nan=1.0)
        dist_matrix = squareform(dist_vector)
        linkage_matrix = linkage(dist_vector, method='ward')

    elif method == 'correlation_nj':
        dist_vector = np.nan_to_num(pdist(matrix, metric='correlation'), nan=1.0)
        dist_matrix = squareform(dist_vector)
        linkage_matrix = linkage(dist_vector, method='average')

        from skbio import DistanceMatrix as SkbioDistanceMatrix
        from skbio.tree import nj
        dm = SkbioDistanceMatrix(dist_matrix, ids=samples)
        nj_tree = nj(dm)
        logger.info(f"    Correlation NJ tree built with {len(samples)} tips")

    elif method == 'correlation_ward':
        dist_vector = np.nan_to_num(pdist(matrix, metric='correlation'), nan=1.0)
        dist_matrix = squareform(dist_vector)
        linkage_matrix = linkage(dist_vector, method='ward')
        logger.info(f"    Correlation+Ward linkage built with {len(samples)} samples")

    else:
        raise ValueError(f"Unknown method '{method}'. Choose from: {METHODS}")

    return dist_vector, dist_matrix, linkage_matrix, nj_tree


def extract_branch_vectors(nj_tree, samples, min_tips=5):
    rooted = nj_tree.root_at_midpoint()
    sample_set = set(samples)
    
    branches = []
    for i, node in enumerate(rooted.traverse(include_self=False)):
        if node.is_tip(): continue
        tips_under = node.subset() & sample_set
        n_tips = len(tips_under)
        if min_tips <= n_tips <= len(samples) - min_tips:
            branches.append({'branch_id': f"B{i}", 'length': node.length or 0.0, 'tips': tips_under})

    if not branches:
        logger.warning("    No internal branches found")
        return []

    result = []
    for b in branches:
        ct = np.array([1 if s in b['tips'] else 0 for s in samples], dtype=np.int8)
        result.append({'branch_id': b['branch_id'], 'length': b['length'], 'n_tips': len(b['tips']), 'vector': ct})

    logger.info(f"    Selected {len(result)} candidate branches for decomposition (tips>={min_tips})")
    return result


def extract_branch_vectors_from_linkage(linkage_matrix, samples, min_tips=5):
    n = len(samples)

    def _get_leaves(node_id):
        if node_id < n:
            return {int(node_id)}
        row = int(node_id - n)
        left = int(linkage_matrix[row, 0])
        right = int(linkage_matrix[row, 1])
        return _get_leaves(left) | _get_leaves(right)

    branches = []
    for row_idx in range(len(linkage_matrix)):
        left = int(linkage_matrix[row_idx, 0])
        right = int(linkage_matrix[row_idx, 1])
        merge_dist = linkage_matrix[row_idx, 2]

        for side_id, child in [('L', left), ('R', right)]:
            leaves = _get_leaves(child)
            n_tips = len(leaves)
            if min_tips <= n_tips <= n - min_tips:
                ct = np.zeros(n, dtype=np.int8)
                for leaf_idx in leaves:
                    ct[leaf_idx] = 1
                branches.append({
                    'branch_id': f"H{row_idx}{side_id}",
                    'length': merge_dist,
                    'n_tips': n_tips,
                    'vector': ct,
                })

    if not branches:
        logger.warning("    No internal branches found from linkage")
        return []

    seen = set()
    unique_branches = []
    for b in branches:
        key = tuple(b['vector'])
        if key not in seen:
            seen.add(key)
            unique_branches.append(b)

    logger.info(f"    Selected {len(unique_branches)} candidate branches from linkage "
                f"(tips>={min_tips})")
    return unique_branches


def _r2_binary(x, y):
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1] ** 2)


def decompose_svs(sv_matrix, branch_vectors, variants_info, samples,
                  r2_threshold=0.3, min_carriers=5, min_haplogroup_size=10):
    if not branch_vectors:
        logger.warning("    No branch vectors; skipping decomposition")
        return [], None, None, {}

    n_samples, n_svs = sv_matrix.shape
    Ct_mat = np.column_stack([b['vector'] for b in branch_vectors])
    branch_ids = [b['branch_id'] for b in branch_vectors]
    n_branches = len(branch_ids)

    records = []
    decomposed_vectors = []
    decomposed_variants = []

    stats = {
        'original_svs': n_svs,
        'decomposed_svs': 0,
        'monophyletic_svs': 0,
        'non_decomposed_svs': 0,
        'residual_svs': 0,
        'original_decomposed_svs': 0
    }

    for j in range(n_svs):
        sv_orig_id = variants_info[j].get('id', f"SV_{j}")
        sv_info = variants_info[j].copy()
        sv_original = sv_matrix[:, j].copy()
        sv = sv_original.copy()
        iteration = 0
        is_decomposed = False

        while True:
            ac = int(sv.sum())
            if ac < min_carriers:
                if ac > 0:
                    records.append({
                        'sv_id': sv_orig_id, 'iteration': iteration,
                        'type': 'residual_small', 'branch_id': None,
                        'r2': None, 'n_carriers': ac,
                        'carrier_indices': np.where(sv > 0)[0].tolist()
                    })

                    if is_decomposed:
                        new_id = f"{sv_orig_id}_residual"
                        new_info = sv_info.copy()
                        new_info.update({
                            'id': new_id,
                            'original_sv_id': sv_orig_id,
                            'mapped_haplogroup': 'residual',
                            'r2_score': None,
                            'sv_category': 'residual',
                            'ac': ac
                        })
                        decomposed_vectors.append(sv.copy())
                        decomposed_variants.append(new_info)
                        stats['residual_svs'] += 1
                    else:
                        new_info = sv_info.copy()
                        new_info.update({
                            'original_sv_id': sv_orig_id,
                            'mapped_haplogroup': None,
                            'r2_score': None,
                            'sv_category': 'non_decomposed',
                            'ac': int(sv_original.sum())
                        })
                        decomposed_vectors.append(sv_original.copy())
                        decomposed_variants.append(new_info)
                        stats['non_decomposed_svs'] += 1
                elif not is_decomposed:
                    new_info = sv_info.copy()
                    new_info.update({
                        'original_sv_id': sv_orig_id,
                        'mapped_haplogroup': None,
                        'r2_score': None,
                        'sv_category': 'non_decomposed',
                        'ac': 0
                    })
                    decomposed_vectors.append(sv_original.copy())
                    decomposed_variants.append(new_info)
                    stats['non_decomposed_svs'] += 1
                break

            r2_vals = np.array([_r2_binary(sv, Ct_mat[:, t]) for t in range(n_branches)])
            best_t = int(np.argmax(r2_vals))
            best_r2 = r2_vals[best_t]

            if best_r2 < r2_threshold:
                records.append({
                    'sv_id': sv_orig_id, 'iteration': iteration,
                    'type': 'paraphyletic', 'branch_id': None,
                    'r2': float(best_r2), 'n_carriers': ac,
                    'carrier_indices': np.where(sv > 0)[0].tolist()
                })

                if is_decomposed:
                    new_id = f"{sv_orig_id}_residual"
                    stats['residual_svs'] += 1

                    new_info = sv_info.copy()
                    new_info.update({
                        'id': new_id,
                        'original_sv_id': sv_orig_id,
                        'mapped_haplogroup': 'residual',
                        'r2_score': float(best_r2),
                        'sv_category': 'residual',
                        'ac': ac
                    })
                    decomposed_vectors.append(sv.copy())
                    decomposed_variants.append(new_info)
                else:
                    new_info = sv_info.copy()
                    new_info.update({
                        'original_sv_id': sv_orig_id,
                        'mapped_haplogroup': None,
                        'r2_score': float(best_r2),
                        'sv_category': 'non_decomposed',
                        'ac': int(sv_original.sum())
                    })
                    decomposed_vectors.append(sv_original.copy())
                    decomposed_variants.append(new_info)
                    stats['non_decomposed_svs'] += 1

                break

            ct_best = Ct_mat[:, best_t]
            bid = branch_ids[best_t]
            mono = sv & ct_best
            n_mono = int(mono.sum())

            if iteration == 0 and n_mono < min_haplogroup_size:
                new_info = sv_info.copy()
                new_info.update({
                    'original_sv_id': sv_orig_id,
                    'mapped_haplogroup': None,
                    'r2_score': float(best_r2),
                    'sv_category': 'non_decomposed',
                    'ac': int(sv_original.sum())
                })
                decomposed_vectors.append(sv_original.copy())
                decomposed_variants.append(new_info)
                stats['non_decomposed_svs'] += 1
                break

            if n_mono < min_carriers:
                records.append({
                    'sv_id': sv_orig_id, 'iteration': iteration,
                    'type': 'paraphyletic', 'branch_id': bid,
                    'r2': float(best_r2), 'n_carriers': ac,
                    'carrier_indices': np.where(sv > 0)[0].tolist()
                })

                if is_decomposed:
                    new_id = f"{sv_orig_id}_residual"
                    stats['residual_svs'] += 1

                    new_info = sv_info.copy()
                    new_info.update({
                        'id': new_id,
                        'original_sv_id': sv_orig_id,
                        'mapped_haplogroup': 'residual',
                        'r2_score': float(best_r2),
                        'sv_category': 'residual',
                        'ac': ac
                    })
                    decomposed_vectors.append(sv.copy())
                    decomposed_variants.append(new_info)
                else:
                    new_info = sv_info.copy()
                    new_info.update({
                        'original_sv_id': sv_orig_id,
                        'mapped_haplogroup': None,
                        'r2_score': float(best_r2),
                        'sv_category': 'non_decomposed',
                        'ac': int(sv_original.sum())
                    })
                    decomposed_vectors.append(sv_original.copy())
                    decomposed_variants.append(new_info)
                    stats['non_decomposed_svs'] += 1

                break

            records.append({
                'sv_id': sv_orig_id, 'iteration': iteration,
                'type': 'monophyletic', 'branch_id': bid,
                'r2': float(best_r2), 'n_carriers': n_mono,
                'carrier_indices': np.where(mono > 0)[0].tolist()
            })

            new_id = f"{sv_orig_id}_{bid}"
            new_info = sv_info.copy()
            new_info.update({
                'id': new_id,
                'original_sv_id': sv_orig_id,
                'mapped_haplogroup': bid,
                'r2_score': float(best_r2),
                'sv_category': 'monophyletic',
                'ac': n_mono
            })

            decomposed_vectors.append(mono)
            decomposed_variants.append(new_info)
            stats['monophyletic_svs'] += 1

            if not is_decomposed:
                is_decomposed = True
                stats['decomposed_svs'] += 1

                orig_info = sv_info.copy()
                orig_ac = int(sv_original.sum())
                orig_info.update({
                    'id': f"{sv_orig_id}_original",
                    'original_sv_id': sv_orig_id,
                    'mapped_haplogroup': None,
                    'r2_score': None,
                    'sv_category': 'original_decomposed',
                    'ac': orig_ac
                })
                decomposed_vectors.append(sv_original.copy())
                decomposed_variants.append(orig_info)
                stats['original_decomposed_svs'] += 1

            sv = sv & (~ct_best.astype(bool)).astype(np.int8)
            iteration += 1

    if decomposed_vectors:
        decomposed_matrix = np.column_stack(decomposed_vectors)
    else:
        decomposed_matrix = np.empty((n_samples, 0))

    return records, decomposed_matrix, decomposed_variants, stats


def parse_vcf_region(vcf_path, region, min_samples=2, max_freq=1.0, overlap_samples=None, dedup=True):
    reader = cyvcf2.VCF(vcf_path, gts012=True)
    samples = list(reader.samples)
    
    if overlap_samples:
        valid_samples = []
        valid_indices = []
        for i, s in enumerate(samples):
            info = extract_sample_info(s)
            if info and info[0] in overlap_samples:
                valid_samples.append(s)
                valid_indices.append(i)
        
        if not valid_samples:
            logger.warning(f"No overlap samples found in VCF for {region}")
            return None, None, []
            
        samples = valid_samples
        sample_indices = np.array(valid_indices)
    else:
        sample_indices = np.arange(len(samples))

    n_samples = len(samples)
    max_count = int(n_samples * max_freq)
    presence_list, variants_info = [], []

    for v in reader(region):
        gt_types = v.gt_types[sample_indices]
        is_present = (gt_types == 1) | (gt_types == 2)
        ac = np.sum(is_present)
        if ac < min_samples or ac > max_count:
            continue
        presence_list.append(is_present.astype(np.int8))
        variants_info.append({
            'chrom': v.CHROM, 'pos': v.POS,
            'end': v.INFO.get('END', v.POS),
            'id': v.ID or f"{v.CHROM}_{v.POS}",
            'svtype': v.INFO.get('SVTYPE', 'UNK'),
            'svlen': abs(v.INFO.get('SVLEN', 0) or 0),
            'ref': v.REF,
            'alt': ",".join(v.ALT) if v.ALT else '.',
            'ac': ac,
            'an': n_samples,
            'af': round(ac / n_samples, 6) if n_samples > 0 else 0,
            'maf': round(min(ac / n_samples, 1.0 - (ac / n_samples)), 6) if n_samples > 0 else 0,
            'n_miss': 0
        })
    reader.close()

    if not presence_list:
        logger.warning(f"No variants found in {region}")
        return None, None, samples

    if dedup:
        presence_list, variants_info, _ = deduplicate_svs(
            presence_list, variants_info, logger=logger
        )
        if not presence_list:
            logger.warning(f"No variants remain after deduplication in {region}")
            return None, None, samples

    return np.vstack(presence_list).T, variants_info, samples


def deduplicate_svs(genotype_list, variants_info, logger=None):
    from collections import defaultdict

    groups = defaultdict(list)
    for i, v in enumerate(variants_info):
        key = (v['pos'], v.get('end', v['pos']), v.get('svtype', 'UNK'), v.get('svlen', 0))
        groups[key].append(i)

    keep_indices = []
    n_groups_with_dups = 0
    for key, indices in groups.items():
        if len(indices) == 1:
            keep_indices.append(indices[0])
        else:
            n_groups_with_dups += 1
            best_idx = max(indices, key=lambda i: variants_info[i].get('ac', 0))
            keep_indices.append(best_idx)

    keep_indices.sort()

    dedup_gt_list = [genotype_list[i] for i in keep_indices]
    dedup_variants = [variants_info[i] for i in keep_indices]

    stats = {
        'n_before': len(variants_info),
        'n_after': len(dedup_variants),
        'n_removed': len(variants_info) - len(dedup_variants),
        'n_groups_with_dups': n_groups_with_dups,
    }

    if logger and stats['n_removed'] > 0:
        logger.info(f"    Dedup: {stats['n_before']} -> {stats['n_after']} SVs "
                    f"({stats['n_removed']} duplicates removed from "
                    f"{stats['n_groups_with_dups']} duplicate groups)")

    return dedup_gt_list, dedup_variants, stats


def perform_clustering(matrix, samples, k_range, min_cluster_size=5, method='jaccard_nj'):
    valid_indices = np.sum(matrix, axis=1) > 0
    if not np.any(valid_indices):
        logger.warning("No samples have any variant presence. Skipping.")
        return None

    matrix = matrix[valid_indices]
    samples = [s for i, s in enumerate(samples) if valid_indices[i]]

    dist_vector, dist_matrix, linkage_matrix, nj_tree = _compute_distance_and_linkage(
        matrix, samples, method
    )

    k_min, k_max = k_range
    safe_max = min(k_max, len(samples) // min_cluster_size)
    best_k, best_score = 2, -1.0
    k_scores = {}

    for k in range(k_min, safe_max + 1):
        labels = fcluster(linkage_matrix, k, criterion='maxclust')
        unique, counts = np.unique(labels, return_counts=True)
        if np.min(counts) < min_cluster_size:
            logger.info(f"    K={k} Rejected: min cluster size {np.min(counts)}")
            break
        if len(unique) < 2:
            continue
        score = silhouette_score(dist_matrix, labels, metric='precomputed')
        k_scores[k] = score
        logger.info(f"    K={k} Silhouette={score:.4f}")
        if score > best_score:
            best_score = score
            best_k = k

    if not k_scores:
        logger.warning(f"  No valid K found. Defaulting to K={best_k}.")

    return {
        'samples': samples,
        'labels': fcluster(linkage_matrix, best_k, criterion='maxclust'),
        'linkage': linkage_matrix,
        'k_scores': k_scores,
        'best_k': best_k,
        'matrix': matrix,
        'nj_tree': nj_tree,
        'method': method
    }


def extract_sample_info(sample_name):
    match = re.match(r'^([A-Za-z0-9]+)_(hap[12]|mat|pat)', sample_name)
    if match:
        s_id = match.group(1)
        hap_str = match.group(2)
        hap = 'hap1' if hap_str in ['hap1', 'mat'] else 'hap2'
        return s_id, hap


def convert_to_diploid(haploid_matrix, haploid_samples):
    indiv_to_idx = {}
    for i, s in enumerate(haploid_samples):
        info = extract_sample_info(s)
        if info:
            s_id = info[0]
            indiv_to_idx.setdefault(s_id, []).append(i)
    
    diploid_samples = sorted(list(indiv_to_idx.keys()))
    n_indiv = len(diploid_samples)
    n_vars = haploid_matrix.shape[1]
    
    diploid_matrix = np.zeros((n_indiv, n_vars), dtype=np.int8)
    
    for i, s_id in enumerate(diploid_samples):
        indices = indiv_to_idx[s_id]
        diploid_matrix[i, :] = np.sum(haploid_matrix[indices, :], axis=0)
        
    return diploid_matrix, diploid_samples


def plot_sv_pie_chart(variants_info, chrom, output_path):
    df = pd.DataFrame(variants_info)
    sv_counts = df['svtype'].value_counts()
    if len(sv_counts) == 0:
        return
    percentages = 100 * sv_counts / sv_counts.sum()
    labels = [f"{svtype} ({pct:.1f}%)" for svtype, pct in zip(sv_counts.index, percentages)]
    plt.figure(figsize=(8, 8))
    plt.pie(sv_counts, labels=labels, autopct='', startangle=90,
            colors=plt.cm.terrain(np.linspace(0, 1, len(sv_counts))))
    plt.title(f"Distribution of Structural Variant Types\n{chrom}")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


YLABEL_MAP = {
    'jaccard_nj': 'Jaccard Distance (UPGMA Linkage)',
    'jaccard_ward': 'Jaccard Distance (Ward Linkage)',
    'cosine_nj': 'Cosine Distance (UPGMA Linkage)',
    'cosine_ward': 'Cosine Distance (Ward Linkage)',
    'correlation_nj': 'Correlation Distance (UPGMA Linkage)',
    'correlation_ward': 'Correlation Distance (Ward Linkage)',
}


def save_outputs(results, variants_info, matrix, chrom, output_dir,
                 decomp_records=None, branch_vectors=None,
                 decomp_matrix=None, decomp_variants=None):
    out_path = Path(output_dir) / chrom
    out_path.mkdir(parents=True, exist_ok=True)
    method = results['method']
    n_vars = len(variants_info)

    pd.DataFrame(variants_info).to_csv(out_path / f"{chrom}_variants.tsv", sep='\t', index=False)
    
    if decomp_matrix is not None and decomp_variants is not None:
        ns, nsv = decomp_matrix.shape
        if nsv > 0:
            diploid_matrix, diploid_samples = convert_to_diploid(decomp_matrix, results['samples'])
            
            n_indiv = len(diploid_samples)
            an = 2 * n_indiv
            for i, v in enumerate(decomp_variants):
                ac = int(np.sum(diploid_matrix[:, i]))
                af = ac / an if an > 0 else 0
                maf = min(af, 1.0 - af)
                v['ac'] = ac
                v['an'] = an
                v['af'] = round(af, 6)
                v['maf'] = round(maf, 6)
            
            np.savez_compressed(out_path / f"{chrom}_sv_genotypes.npz", genotypes=diploid_matrix.T)
            
            with open(out_path / f"{chrom}_sv_samples.txt", 'w') as fh:
                fh.write('\n'.join(diploid_samples) + '\n')
            
            pd.DataFrame(decomp_variants).to_csv(out_path / f"{chrom}_sv_variants.tsv", sep='\t', index=False)
            logger.info(f"    Saved {nsv} decomposed pure SVs (diploid) to {chrom}_sv_genotypes.npz")
        else:
            logger.warning("    Decomposition was run but yielded no monophyletic SVs.")

    nj_tree = results.get('nj_tree')

    if decomp_records:
        df_decomp = pd.DataFrame(decomp_records)
        mono_counts = df_decomp[df_decomp['type'] == 'monophyletic']['branch_id'].value_counts()
        
        validated_hgs = []
        if branch_vectors:
             for b in branch_vectors:
                bid = b['branch_id']
                count = mono_counts.get(bid, 0)
                if count > 10:
                    vector = b['vector']
                    carriers = [s for i, s in enumerate(results['samples']) if vector[i] == 1]
                    
                    validated_hgs.append({
                        'haplogroup_id': bid,
                        'monophyletic_svs': count, 
                        'n_carriers': len(carriers),
                        'carriers': ",".join(carriers)
                    })
        
        if validated_hgs:
             pd.DataFrame(validated_hgs).sort_values('monophyletic_svs', ascending=False).to_csv(out_path / f"{chrom}_validated_haplogroups.tsv", sep='\t', index=False)
             logger.info(f"    Saved {len(validated_hgs)} validated haplogroups to {chrom}_validated_haplogroups.tsv")
        else:
             pd.DataFrame(columns=['haplogroup_id', 'monophyletic_svs', 'n_carriers', 'carriers']).to_csv(
                 out_path / f"{chrom}_validated_haplogroups.tsv", sep='\t', index=False)
             logger.info(f"    No validated haplogroups (>10 mono SVs) found. Created empty output.")

        mono_counts_dict = mono_counts.to_dict()
        if nj_tree:
             visualize_branch_support(nj_tree, mono_counts_dict, out_path / f"{chrom}_branch_support.png")

    data = []
    for sample, label in zip(results['samples'], results['labels']):
        s_id, hap = extract_sample_info(sample)
        data.append({
            'sample': sample, 'sample_id': s_id, 'haplotype': hap,
            'chromosome': chrom,
            'haplogroup': f"C{chrom.replace('chr','')}-aHOR-HG{label}",
            'haplogroup_id': f"{chrom}_{label}"
        })
    df_hg = pd.DataFrame(data)
    df_hg['cluster_size'] = df_hg['haplogroup_id'].map(df_hg['haplogroup_id'].value_counts())
    df_hg.to_csv(out_path / f"{chrom}_haplogroups.tsv", sep='\t', index=False)

    plot_sv_pie_chart(variants_info, chrom, out_path / f"{chrom}_sv_distribution.pdf")

    best_k = results['best_k']
    Z = results['linkage']
    color_threshold = 0
    if best_k > 1 and Z.shape[0] >= best_k:
        color_threshold = 0.5 * (Z[-best_k, 2] + Z[-(best_k-1), 2])
    elif Z.shape[0] > 0:
        color_threshold = 0.5 * Z[0, 2]
    dendrogram(Z, leaf_rotation=90, color_threshold=color_threshold, no_labels=True)
    plt.title(f"{chrom} Haplogroups (K={best_k}, {method})\nVariants: {n_vars}")
    plt.xlabel("Samples")
    plt.ylabel(YLABEL_MAP.get(method, method))
    plt.tight_layout()
    plt.savefig(out_path / f"{chrom}_dendrogram.png")
    plt.close()

    nj_tree_from_results = results.get('nj_tree')
    if nj_tree_from_results is not None:
        nj_tree_from_results.write(str(out_path / f"{chrom}_nj.nwk"), format='newick')
        
        if branch_vectors and decomp_records:
            from skbio.tree import TreeNode
            
            mono_counts = df_decomp[df_decomp['type'] == 'monophyletic']['branch_id'].value_counts().to_dict()
            
            rooted = nj_tree_from_results.root_at_midpoint()
            
            for i, node in enumerate(rooted.traverse(include_self=False)):
                bid = f"B{i}"
                support = mono_counts.get(bid, 0)
                if hasattr(node, 'support'):
                     node.support = support
                else:
                     node.name = f"{node.name}:{support}" if node.name else str(support)

            rooted.write(str(out_path / f"{chrom}_nj_annotated.nwk"), format='newick')
            logger.info(f"    Saved annotated tree to {chrom}_nj_annotated.nwk")

    if branch_vectors:
        bv_data = [{'branch_id': b['branch_id'], 'length': b['length'], 'n_tips': b['n_tips']}
                   for b in branch_vectors]
        pd.DataFrame(bv_data).to_csv(out_path / f"{chrom}_branches.tsv", sep='\t', index=False)

    if decomp_records:
        df_decomp = pd.DataFrame(decomp_records)
        df_decomp['carrier_indices'] = df_decomp['carrier_indices'].apply(
            lambda x: ','.join(map(str, x)) if x else '')
        df_decomp.to_csv(out_path / f"{chrom}_decomposition.tsv", sep='\t', index=False)
        for t, c in df_decomp['type'].value_counts().items():
            logger.info(f"  Decomposition: {t} = {c}")

    if n_vars > 1000:
        idx = np.argsort(np.var(matrix, axis=0))[-1000:]
        matrix_plot = matrix[:, idx]
        title_suffix = f"(top 1000/{n_vars})"
    else:
        matrix_plot = matrix
        title_suffix = f"({n_vars} variants)"
    sort_idx = np.argsort(results['labels'])
    sns.heatmap(matrix_plot[sort_idx], cmap="YlGnBu", xticklabels=False,
                yticklabels=False, cbar_kws={'label': 'Presence'})
    plt.title(f"{chrom} Heatmap {title_suffix}\n(Sorted by Haplogroup)")
    plt.savefig(out_path / f"{chrom}_heatmap.png")
    plt.close()

    if results['k_scores']:
        ks = sorted(results['k_scores'].keys())
        plt.plot(ks, [results['k_scores'][k] for k in ks], 'b-o')
        plt.axvline(best_k, color='red', linestyle='--')
        plt.xlabel("K")
        plt.ylabel("Silhouette Score")
        plt.title(f"{chrom} K={best_k}")
        plt.savefig(out_path / f"{chrom}_silhouette.png")
        plt.close()

    return df_hg, pd.DataFrame(variants_info)


def visualize_branch_support(nj_tree, support_dict, output_path):
    sorted_branches = sorted(support_dict.items(), key=lambda x: x[1], reverse=True)
    top_20 = sorted_branches[:20]
    
    if not top_20:
        return

    b_ids, counts = zip(*top_20)
    
    plt.figure(figsize=(10, 6))
    plt.bar(b_ids, counts, color='teal')
    plt.xlabel('Branch ID')
    plt.ylabel('Number of Monophyletic SVs')
    plt.title('Top Validated Haplogroups (Branch Support)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Centromere Haplogroup Clustering")
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--chrom', required=True, help="Chromosome to process (e.g. chr17)")
    parser.add_argument('--method', default='jaccard_nj', choices=METHODS)
    parser.add_argument('--min-samples', type=int, default=5)
    parser.add_argument('--max-freq', type=float, default=0.99)
    parser.add_argument('--k-range', default='2,20')
    parser.add_argument('--min-cluster-size', type=int, default=5)
    parser.add_argument('--decompose', action='store_true',
                        help="Run iterative SV decomposition (jaccard_nj only)")
    parser.add_argument('--r2-threshold', type=float, default=0.3)
    parser.add_argument('--min-carriers', type=int, default=5)
    parser.add_argument('--min-tips', type=int, default=5,
                        help="Minimum number of tips for decomposition branches (default: 5)")
    parser.add_argument('--overlap', action='store_true', help="Filter samples by overlap list")
    parser.add_argument('--overlap-file', type=str, default="resources/overlap_samples.txt", help="Path to overlap samples file")

    args = parser.parse_args()
    k_min, k_max = map(int, args.k_range.split(','))
    chrom = args.chrom

    overlap_samples = None
    if args.overlap:
        try:
            with open(args.overlap_file) as f:
                overlap_samples = set(line.strip() for line in f if line.strip() and not line.startswith('#'))
            logger.info(f"Loaded {len(overlap_samples)} overlap samples from {args.overlap_file}")
        except Exception as e:
            logger.error(f"Failed to load overlap samples from {args.overlap_file}: {e}")
            sys.exit(1)

    logger.info(f"Processing {chrom} with method={args.method}")
    matrix, variants, samples = parse_vcf_region(
        args.vcf, chrom, args.min_samples, args.max_freq, overlap_samples
    )
    
    if matrix is None:
        logger.error(f"No usable variants found for {chrom}")
        sys.exit(0)

    results = perform_clustering(
        matrix, samples, (k_min, k_max), args.min_cluster_size, args.method
    )
    
    if results is None:
        logger.error(f"Clustering failed for {chrom}")
        sys.exit(0)

    branch_vectors, decomp_records = [], []
    mk_decomp, vk_decomp = None, None

    if args.decompose and results.get('nj_tree') is not None:
        logger.info(f"  Running SV decomposition for {chrom}")
        branch_vectors = extract_branch_vectors(
            results['nj_tree'], results['samples'],
            min_tips=args.min_tips
        )
        decomp_records, mk_decomp, vk_decomp, stats = decompose_svs(
            results['matrix'], branch_vectors, variants, results['samples'],
            r2_threshold=args.r2_threshold, min_carriers=args.min_carriers
        )
        
        if stats:
            logger.info(f"  Decomposition stats for {chrom}:")
            logger.info(f"    Original SVs: {stats['original_svs']}")
            logger.info(f"    Decomposed SVs: {stats['decomposed_svs']}")
            logger.info(f"    Monophyletic SVs created: {stats['monophyletic_svs']}")
            logger.info(f"    Non-decomposed SVs kept: {stats['non_decomposed_svs']}")
            logger.info(f"    Residual SVs kept: {stats['residual_svs']}")
            logger.info(f"    Total SVs after decomposition: {len(vk_decomp)}")

    df_hg, df_var = save_outputs(
        results, variants, results['matrix'], chrom, args.output,
        decomp_records=decomp_records, branch_vectors=branch_vectors,
        decomp_matrix=mk_decomp, decomp_variants=vk_decomp
    )
    
    out = Path(args.output)
    df_hg.to_csv(out / f"{chrom}_all_haplogroups.tsv", sep='\t', index=False)
    df_var.to_csv(out / f"{chrom}_all_variants.tsv", sep='\t', index=False)
    
    if decomp_records:
        df_d = pd.DataFrame(decomp_records)
        df_d['chrom'] = chrom
        df_d['carrier_indices'] = df_d['carrier_indices'].apply(
            lambda x: ','.join(map(str, x)) if isinstance(x, list) else x)
        df_d.to_csv(out / f"{chrom}_all_decomposition.tsv", sep='\t', index=False)

    if vk_decomp:
         pd.DataFrame(vk_decomp).to_csv(out / f"{chrom}_all_sv_variants.tsv", sep='\t', index=False)
        
    logger.info(f"Finished {chrom}: K={results['best_k']}, Variants={len(variants)}")


if __name__ == "__main__":
    main()
