"""Production records interpreted on a fully known background-SNV template."""
import json
from pathlib import Path

import pysam

from native_truth import reverse_complement
from score_sv_recovery import parse_calls


def parse_complete_calls(vcf_path,template,trace_path=None):
    calls=parse_calls(vcf_path,template)
    traces={x['vcf_id']:x for x in json.loads(Path(trace_path).read_text())} if trace_path else {}
    with pysam.VariantFile(str(vcf_path)) as vcf:
        raw={r.id:dict(reported_chrom=r.chrom,reported_pos=r.pos,reported_end=r.stop,
                      reported_svlen=r.info.get('SVLEN'),literal_ref=r.ref,literal_alt=r.alts)
             for r in vcf}
    for call in calls:
        call.update(raw[call['id']]);trace=traces.get(call['id'],{})
        call['literal_edit']=call['edit'];call['candidate_trace']=trace
        call['allele_interpretation']='literal_allele'
        if call['kind']=='INV':
            a=trace.get('source_start',call['reported_pos']-1);b=trace.get('source_end',call['reported_end'])
            call['edit']=(a,b,reverse_complement(template.fetch(call['chrom'],a,b)))
            call['allele_interpretation']='background_SNV_restored_inversion'
        elif call['kind']=='DUP:TANDEM':
            call['allele_interpretation']='background_SNV_restored_tandem_source'
        elif call['kind']=='DUP:INT':
            if trace and not trace.get('cutpaste') and trace.get('trace_key_multiplicity',1)==1:
                source=template.fetch(trace['source_contig'],trace['source_start'],trace['source_end'])
                call['chrom']=trace['dest_contig'];p=trace['dest_start'];call['edit']=(p,p,source)
                call['allele_interpretation']='background_SNV_restored_traced_dispersed_source'
            else:
                call['edit']=None;call['allele_interpretation']='unresolved_missing_destination_or_cutpaste'
        elif call['kind']=='BND':call['allele_interpretation']='junction_record_not_reconstructed'
        if call['edit']:
            a,b,alt=call['edit'];call['length']=len(alt)-(b-a)
    return calls
