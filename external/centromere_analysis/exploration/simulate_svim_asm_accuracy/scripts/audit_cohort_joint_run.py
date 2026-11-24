#!/usr/bin/env python3
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import re
from pathlib import Path
import tempfile

import pysam

from audit_cohort_frequency_run import (AuditFailure, audit_aggregate, audit_resources, audit_score,
    classify, fields, identity_catalog, mixed, read_json, require, same, sha, table_matches)
from cohort_frequency_joint_reference import joint_truth

THRESHOLDS = (0, 1, 5, 10)


JOINT_IDENTITY_REQUIRED = {
    'prepare_joint_cohort.py', 'cohort_frequency_joint_reference.py',
    'score_cohort_joint_truth.py', 'cohort_bnd_junction_validation.py',
    'multichrom_truth_scoring.py', 'cohort_frequency_joint_metrics.py',
    'validate_cohort_joint_semantics.py', 'audit_cohort_joint_run.py',
}


def validate_joint_identity(root, evidence):
    root = Path(root)
    path = root / 'joint_execution_code_identity.json'
    identity = read_json(path, evidence)
    require(isinstance(identity, dict) and identity and
            all(isinstance(v, str) and re.fullmatch(r'[0-9a-f]{64}', v) for v in identity.values()),
            'Invalid joint execution identity map')
    require(JOINT_IDENTITY_REQUIRED <= {Path(name).name for name in identity},
            'Joint execution identity omits required source files')
    archived_hashes = set()
    for folder in (root/'code_snapshots', root/'joint_code_archive', root/'joint_pilot'):
        if folder.exists():
            for suffix in ('*.py', '*.sh', '*.json'):
                archived_hashes.update(sha(p) for p in folder.rglob(suffix) if p.is_file())
    for source, digest in identity.items():
        source = Path(source)
        require((source.is_file() and sha(source) == digest) or digest in archived_hashes,
                f'Joint execution source missing or changed without archive: {source}')
    return sha(path)


def vcf_fields_by_id(path):
    result = {}
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            fields = str(record).rstrip('\n').split('\t')
            if not record.id or record.id in result:
                raise ValueError(f'Invalid VCF ID in {path}: {record.id}')
            result[record.id] = tuple(fields[:2] + fields[3:])
    return result



def validate_false_junction_receipt(check, call_id, truth, query_sha256_by_chromosome=None):
    require(check.get('call_id') == call_id and check.get('status') == 'definite_false_junction' and
            check.get('definite_fp') is True and check.get('within_domain_contradiction') is True and
            check.get('complete_endpoint_chromosome_absence') is True,
            'BND receipt lacks definite local and whole-chromosome contradiction')
    endpoints = check.get('endpoints', [])
    require(len(endpoints) == 2 and all(len(e) == 3 for e in endpoints), 'Invalid BND endpoints')
    windows = check.get('windows', [])
    require(len(windows) == 2 and [w['chrom'] for w in windows] == [e[0] for e in endpoints],
            'BND local search does not cover both endpoints')
    require(all(w.get('forward_match_offset') == -1 and w.get('reverse_match_offset') == -1 and
                re.fullmatch(r'[0-9a-f]{64}', w.get('edited_sha256', '')) for w in windows),
            'BND local absence evidence is incomplete')
    junction = check.get('junction_sequence', '')
    require(bool(junction) and hashlib.sha256(junction.encode()).hexdigest() == check.get('junction_sha256'),
            'BND junction sequence digest differs')
    complete = check.get('complete_chromosomes', [])
    endpoint_chroms = {e[0] for e in endpoints}
    require(len(complete) == len(endpoint_chroms) and {c['chrom'] for c in complete} == endpoint_chroms,
            'BND whole-chromosome scan does not cover all endpoint chromosomes')
    for chrom in complete:
        events = [e for e in truth if e['chrom'] == chrom['chrom']]
        expected_length = chrom['reference_length'] + sum(len(e['alt']) - (e['end']-e['start']) for e in events)
        require(chrom.get('forward_match_offset') == -1 and chrom.get('reverse_match_offset') == -1 and
                chrom.get('ambiguous_bases_present') is False and chrom.get('edited_length') == expected_length and
                set(chrom.get('truth_event_ids', [])) == {e['id'] for e in events} and
                re.fullmatch(r'[0-9a-f]{64}', chrom.get('edited_sha256', '')),
                'BND whole-chromosome absence evidence is incomplete')
        if query_sha256_by_chromosome and chrom['chrom'] in query_sha256_by_chromosome:
            same(chrom['edited_sha256'], query_sha256_by_chromosome[chrom['chrom']],
                 'BND endpoint full edited chromosome hash')


def validate_joint_resources(joint_run, manifest):
    n = manifest['config']['n_haplotypes']
    chroms = manifest['config']['chromosomes']
    fields(joint_run, dict(completed=True, haplotypes=n, source_queries=n*len(chroms)), 'Joint resource receipt')
    job, node, cpus = joint_run.get('job'), joint_run.get('node'), joint_run.get('requested_cpus')
    require(isinstance(job, str) and re.fullmatch(r'[0-9]+', job) and isinstance(node, str) and
            node.split('.')[0] in set(manifest['config']['node_allowlist']) and
            isinstance(cpus, str) and cpus.isdigit() and 0 < int(cpus) <= 4,
            'Invalid joint job/node/CPU receipt')
    return job, node, cpus

def run_audit(root, output, semantic_receipt=None):
    root, output = Path(root).resolve(), Path(output).resolve()
    receipt = dict(passed=False, status='in_progress', source_sha256={},
                   started_utc=datetime.now(timezone.utc).isoformat())
    output.parent.mkdir(parents=True, exist_ok=True)
    def save():
        output.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    save()
    try:
        evidence = receipt['source_sha256']
        manifest = read_json(root/'manifest.json', evidence)
        chroms = manifest['config']['chromosomes']; n = manifest['config']['n_haplotypes']
        require(n == 256 and set(chroms) == {'chr1','chr7','chr13','chr19'}, 'Unexpected cohort design')
        require(read_json(root/'batch_equivalence.json', evidence).get('passed') is True, 'Mapping gate failed')
        epochs = identity_catalog(root, manifest, evidence)
        joint_identity = validate_joint_identity(root, evidence)
        joint_run = read_json(root/'joint_cohort/complete.json', evidence)
        job, node, cpus = validate_joint_resources(joint_run, manifest)
        plans = {c: {e['id']: e for e in manifest['plans'][c]['events']} for c in chroms}
        require(all(len(p)==32 for p in plans.values()), 'Incorrect site count')
        recipe_counts = defaultdict(Counter)
        for plan in plans.values():
            for event in plan.values():
                require(len(event['carriers']) == len(set(event['carriers'])) == event['carrier_count'] and
                        all(isinstance(h,int) and 0<=h<n for h in event['carriers']), 'Invalid truth carriers')
                recipe_counts[event['recipe']][event['carrier_count']] += 1
        require(len(recipe_counts)==16 and all(v==Counter([1,2,4,8,16,32,64,128])
                                                for v in recipe_counts.values()), 'Unbalanced event-frequency design')
        for chrom in chroms:
            founder_dir = root/'chromosomes'/chrom/'founders'
            fields(read_json(founder_dir/'complete.json', evidence), dict(completed=True,founders=8), 'Founder completion')
            for index in range(8):
                founder = read_json(founder_dir/f'F{index}.json', evidence)
                require(founder.get('ledger_replay_verified') is True, 'Founder SNP replay proof missing')
                require(re.fullmatch(r'[0-9a-f]{64}', founder['sequence_sha256']), 'Invalid founder sequence digest')
        global_data = read_json(root/'joint_analysis/analysis.json', evidence)
        rows, members, event_carriers = (global_data[k] for k in ('rows','memberships','event_carriers'))
        require([r['row_index'] for r in rows] == list(range(len(rows))), 'Noncontiguous global rows')
        grouped, ids, retained = defaultdict(list), set(), defaultdict(set)
        for member in members:
            identity, h, index = member['source_id'], member['haplotype_index'], member['row_index']
            require(identity not in ids and identity.startswith(f'H{h:03d}__') and 0 <= h < n and 0<=index<len(rows),
                    'Invalid or duplicate joint source membership')
            ids.add(identity); grouped[index].append(member)
            row = rows[index]
            require(member['chrom'] == row['chrom'], 'Membership reported chromosome differs')
            for p in THRESHOLDS:
                if row['dedup_retained'] and row['maf'] >= p/100:
                    retained[h,p].add(identity.split('__',1)[1])
        nodes, epoch_counts, all_source_ids = {node.split('.')[0]}, Counter(), set()
        event_scores, call_scores, template_hashes = {}, {}, {c:set() for c in chroms}
        for h in range(n):
            joint = root/'joint_cohort/haps'/f'H{h:03d}'
            done, truth = read_json(joint/'complete.json', evidence), read_json(joint/'truth.json', evidence)
            fields(done, dict(completed=True,haplotype_index=h,manifest_sha256=sha(root/'manifest.json'),
                              joint_execution_code_identity_sha256=joint_identity,
                              job=job, node=node, requested_cpus=cpus), 'Joint completion')
            same(truth, joint_truth(root, manifest, h), f'Joint truth union H{h:03d}')
            union_calls = []
            union_vcf_fields = {}
            for chrom in chroms:
                source = root/'chromosomes'/chrom/'haps'/f'H{h:03d}'
                source_done = read_json(source/'complete.json', evidence)
                source_truth = read_json(source/'truth.json', evidence)
                oracle = read_json(source/'oracle.json', evidence)
                require(source_done['completed'] is True and oracle.get('passed') is True, f'Incomplete source {source}')
                oracle_checks = ('all_truth_recovered', 'no_false_or_unresolved_calls', 'no_unresolved_events',
                                 'whole_sequence_reconstructed', 'no_false_negative_windows')
                require(all(oracle.get('checks', {}).get(key) is True for key in oracle_checks),
                        f'Source oracle checks incomplete: {source}')
                require(source_truth.get('ledger_replay_verified') is True, f'Missing SNV replay {source}')
                fields(source_truth, dict(id=f'{chrom}_H{h:03d}',chrom=chrom,haplotype_index=h,
                    sample_name=f'SIM{h//2:03d}_hap{h%2+1}'), 'Source truth identity')
                require(all(re.fullmatch(r'[0-9a-f]{64}', source_truth[k])
                            for k in ('template_sha256', 'query_sha256')), 'Invalid source sequence digest')
                fields(source_done, dict(case_id=source_truth['id'],query_sha256=source_truth['query_sha256'],
                                         template_sha256=source_truth['template_sha256']), 'Source sequence identity')
                digest = source_truth['template_sha256']
                require(digest not in template_hashes[chrom], f'Duplicate source template {chrom}/H{h:03d}')
                template_hashes[chrom].add(digest)
                same(done['template_sha256_by_chromosome'][chrom], digest, 'Joint verified template identity')
                wanted = {name for name,e in plans[chrom].items() if h in e['carriers']}
                require(len(source_truth['truth']) == len(wanted) and
                        {e['id'] for e in source_truth['truth']} == wanted, 'Source carrier assignment differs')
                for event in source_truth['truth']:
                    fields(event, {k: v for k, v in plans[chrom][event['id']].items() if k != 'carriers'},
                           'Frozen structural event design')
                    same(hashlib.sha256(event['alt'].encode()).hexdigest(), event['inserted_sha256'],
                         'Inserted payload SHA-256')
                fields(oracle['metrics'], dict(n_truth=len(wanted), tp_events=len(wanted), fp_calls=0,
                                               unresolved_calls=0, full_sequence_reconstruction=True),
                       'Source perfect-call oracle')
                original_calls = read_json(source/'scoring/parsed_calls_production_active.json', evidence)
                original_vcf = source/'active.vcf.gz'
                fields(done['source_hashes'][chrom], dict(vcf_sha256=sha(original_vcf),
                    parsed_sha256=sha(source/'scoring/parsed_calls_production_active.json')), 'Source union hashes')
                union_calls.extend(dict(c,id=f"{chrom}__{c['id']}",source_query_chrom=chrom,
                                        source_record_id=c['id']) for c in original_calls)
                for original_id, value in vcf_fields_by_id(original_vcf).items():
                    identity=f'{chrom}__{original_id}'
                    require(identity not in union_vcf_fields, 'Duplicate namespaced VCF source')
                    union_vcf_fields[identity]=value
                nodes.add(source_done['node'].split('.')[0])
                epoch = source_done.get('execution_code_identity_sha256') or source_done.get('pilot_generation_code_identity_sha256')
                require(epoch in epochs, f'Unarchived execution identity: {source}')
                epoch_counts[epoch] += 1
            parsed = read_json(joint/'scoring/parsed_calls_production_active.json', evidence)
            same(parsed, union_calls, f'Joint parsed union H{h:03d}')
            same(vcf_fields_by_id(joint/'active.vcf.gz'), union_vcf_fields,
                 f'Joint VCF preserves every source field H{h:03d}')
            expected_prefixed = {f'H{h:03d}__{key}': value for key, value in union_vcf_fields.items()}
            same(vcf_fields_by_id(joint/'merge_input.vcf.gz'), expected_prefixed,
                 f'Merge input preserves joint source fields H{h:03d}')
            calls = [c['id'] for c in parsed]
            require(len(calls)==len(set(calls))==done['records'], 'Invalid joint source IDs')
            all_source_ids.update(f'H{h:03d}__{i}' for i in calls)
            truth_ids = {e['id'] for e in truth['truth']}
            for p in ('raw',*THRESHOLDS):
                name = 'production_active.json' if p=='raw' else f'filtered_maf_{p}.json'
                score = read_json(joint/'scoring'/name, evidence)
                metadata = score.get('metadata', {})
                fields(metadata, dict(joint_execution_code_identity_sha256=joint_identity,
                                      template_sha256_by_chromosome=truth['template_sha256_by_chromosome']),
                       'Joint score execution identity')
                if p != 'raw':
                    fields(metadata, dict(haplotype_index=h, maf_filter_percent=p), 'Joint filtered score identity')
                call_ids = calls if p=='raw' else retained[h,p]
                call_scores[h,p], event_scores[h,p] = audit_score(score,truth_ids,call_ids,manifest['strict_scoring_parameters'])
                require(score['conventions'].get('junction_contradiction_version') == 'exact_bnd_junction_contradiction_v1',
                        'Joint BND junction validation version missing')
                junction = score.get('junction_validation', {})
                excluded = junction.get('excluded_from_recovery_search', [])
                require(len(excluded) == len(set(excluded)) == junction.get('n_proven_false_records') and
                        set(excluded) <= set(call_ids), 'Invalid proven-false junction count')
                for call_id in excluded:
                    check = junction.get('receipts', {}).get(call_id, {})
                    require(call_scores[h,p][call_id] == 'fp', 'Proven BND was not scored false')
                    validate_false_junction_receipt(check, call_id, truth['truth'], truth['query_sha256_by_chromosome'])
                if p=='raw':
                    same(done['metrics'],score['metrics'],'Joint raw metrics')
        require(ids == all_source_ids, 'Merged/source membership coverage differs')
        winners={}
        for row in rows:
            group=grouped[row['row_index']]
            require(set(row['source_ids']) == {m['source_id'] for m in group}, 'Row source IDs disagree')
            require(len(group)==len({m['haplotype_index'] for m in group})==row['ac'], 'Joint carrier AC differs')
            fields(row, dict(an=n,af=row['ac']/n,maf=min(row['ac']/n,1-row['ac']/n)), 'Joint row frequency')
            raw=[call_scores[m['haplotype_index'],'raw'][m['source_id'].split('__',1)[1]] for m in group]
            fields(row,{'raw_'+k:v for k,v in classify(raw).items()},'Raw joint row status')
            passes=row['ac']>=2 and row['ac']<=int(n*.99)
            same(row['passes_ac_af'],passes,'Production AC/AF')
            if passes:
                key=(row['chrom'],row['pos'],row['end'],row['svtype'],row['svlen'])
                if key not in winners or row['ac']>rows[winners[key]]['ac']:winners[key]=row['row_index']
            for p in THRESHOLDS:
                survives=row['dedup_retained'] and row['maf']>=p/100
                require((str(p) in row['rowstatus_by_maf'])==survives,'Wrong threshold row status presence')
                if survives:
                    statuses=[call_scores[m['haplotype_index'],p][m['source_id'].split('__',1)[1]] for m in group]
                    fields(row['rowstatus_by_maf'][str(p)],classify(statuses),'Filtered row status')
                    if p==0: fields(row,classify(statuses),'Baseline row status')
                for member in group:
                    value=call_scores[member['haplotype_index'],p][member['source_id'].split('__',1)[1]] if survives else None
                    same(member[f'status_maf_{p}'],value,'Filtered membership status')
        require({r['row_index'] for r in rows if r['dedup_retained']} == set(winners.values()), 'Joint dedup differs')
        wanted={(chrom,f'{chrom}__{e["id"]}',h) for chrom in chroms for e in plans[chrom].values() for h in e['carriers']}
        actual={(r['chrom'],r['event_id'],r['haplotype_index']) for r in event_carriers}
        require(len(actual)==len(event_carriers)==len(wanted) and actual==wanted, 'Truth event carrier coverage differs')
        for e in event_carriers:
            h,event=e['haplotype_index'],e['event_id']
            fields(e,dict(status=event_scores[h,'raw'][event]['status'],raw_recovered=event_scores[h,'raw'][event]['recovered']), 'Raw truth recovery')
            for p in THRESHOLDS:
                fields(e,{f'status_maf_{p}':event_scores[h,p][event]['status'],
                          f'recovered_maf_{p}':event_scores[h,p][event]['recovered']}, 'Filtered truth recovery')
        for name, data in [('merged_rows.tsv',rows),('source_membership.tsv',members),('event_carriers.tsv',event_carriers)]:
            table_matches(root/'joint_analysis'/name,data,evidence)
        for chrom in sorted(set(chroms)|{r['chrom'] for r in rows}):
            part=read_json(root/'joint_analysis'/chrom/'analysis.json',evidence)
            same(part['rows'],[r for r in rows if r['chrom']==chrom], 'Reported-chromosome export')
        analyses=[]
        for chrom in chroms:
            selected=[r for r in rows if r['chrom']==chrom]
            if chrom==chroms[0]:selected += [r for r in rows if r['chrom'] not in chroms]
            analyses.append(dict(rows=selected,event_carriers=[r for r in event_carriers if r['chrom']==chrom]))
        with tempfile.TemporaryDirectory(prefix='audit_joint_') as tmp:
            t=Path(tmp); (t/'analysis').mkdir(); (t/'analysis/merged').symlink_to(root/'joint_analysis/merged',target_is_directory=True)
            receipt['reconciliation']=audit_aggregate(t,manifest,analyses,evidence)
        semantics=read_json(Path(semantic_receipt) if semantic_receipt else root/'joint_analysis/merged_semantics_validation.json',evidence)
        require(semantics.get('passed') is True, 'Joint merged semantics failed')
        fields(semantics,dict(manifest_sha256=sha(root/'manifest.json'),n_haplotypes=n,
            merged_vcf_sha256=sha(root/'joint_analysis/merged.vcf.gz'),
            validator_sha256=sha(Path(__file__).with_name('validate_cohort_joint_semantics.py'))), 'Joint semantic identity')
        require(set(semantics['reported_chromosomes'])=={r['chrom'] for r in rows}, 'Semantic reported chromosomes differ')
        fields(semantics['counts'],dict(merged_rows=len(rows),source_carrier_assertions=len(members)), 'Semantic source counts')
        expected={str((root/'joint_cohort/haps'/f'H{h:03d}'/'merge_input.vcf.gz').resolve()) for h in range(n)}
        require(set(semantics['source_vcf_sha256'])==expected, 'Incomplete semantic source VCFs')
        for path,digest in semantics['source_vcf_sha256'].items():same(sha(path),digest,'Stale semantic source VCF')
        receipt['resources']=audit_resources(read_json(root/'submissions.json',evidence),manifest['config'],nodes)
        receipt.update(passed=True,status='complete',source_cases=n*len(chroms),joint_cases=n,
                       reported_chromosomes=sorted({r['chrom'] for r in rows}),execution_epochs=dict(epoch_counts),joint_execution_identity_sha256=joint_identity,
                       auditor_sha256=sha(__file__),finished_utc=datetime.now(timezone.utc).isoformat())
        save(); return receipt
    except Exception as error:
        receipt.update(passed=False,status='failed',error=f'{type(error).__name__}: {error}',finished_utc=datetime.now(timezone.utc).isoformat())
        save(); raise AuditFailure(receipt['error']) from error


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',default='results/cohort_frequency_v1')
    parser.add_argument('--output')
    parser.add_argument('--semantic-receipt')
    args=parser.parse_args()
    result=run_audit(args.root,args.output or str(Path(args.root)/'joint_analysis/final_run_audit.json'),args.semantic_receipt)
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}))


if __name__=='__main__':main()
