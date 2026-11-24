#!/usr/bin/env python3
"""Instrument output only: preserve the production VCF and retain candidate origins."""
from collections import defaultdict,deque
import json
import os
from pathlib import Path
import runpy
import sys
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]


def entry_key(entry):
    fields=entry.rstrip('\n').split('\t')
    fields[2]='*'
    return tuple(fields)


def capture_write(original,*args,**kwargs):
    (intdup,inv,tandem,deletions,insertions,breakends,version,names,lengths,types,reference,options)=args
    if options.tandem_duplications_as_insertions or options.interspersed_duplications_as_insertions:
        raise ValueError('Trace adapter expects unchanged production output flags')
    sequence=not options.symbolic_alleles
    lookup=defaultdict(deque)
    for kind,candidates in [('DEL',deletions),('INV',inv),('INS',insertions),
                            ('DUP:TANDEM',tandem),('DUP:INT',intdup),('BND',breakends)]:
        if kind not in types:continue
        for candidate_index,candidate in enumerate(candidates):
            data={k:v for k,v in vars(candidate).items() if isinstance(v,(str,int,float,bool,list,tuple,type(None)))}
            data['candidate_kind']=kind
            data['candidate_group']=f'{kind}:{candidate_index}'
            if kind in ['DEL','INV','INS']:
                entries=[candidate.get_vcf_entry(sequence,reference,options.query_names)]
            elif kind in ['DUP:TANDEM','DUP:INT']:
                entries=[candidate.get_vcf_entry_as_dup(options.query_names)]
            else:
                entries=[candidate.get_vcf_entry(options.query_names),candidate.get_vcf_entry_reverse(options.query_names)]
            for entry in entries:lookup[entry_key(entry)].append(dict(data))
    multiplicities={k:len(v) for k,v in lookup.items()}
    original(*args,**kwargs)
    records=[]
    with (Path(options.working_dir)/'variants.vcf').open() as handle:
        for line in handle:
            if line.startswith('#'):continue
            key=entry_key(line)
            if not lookup[key]:raise ValueError('Trace entry does not match production VCF')
            trace=lookup[key].popleft()
            trace.update(vcf_id=line.split('\t')[2],trace_key_multiplicity=multiplicities[key])
            records.append(trace)
    if any(lookup.values()):raise ValueError('Unwritten trace candidates')
    (Path(options.working_dir)/'candidate_trace.json').write_text(json.dumps(records,indent=2)+'\n')


def main():
    source=Path(os.environ.get('SVIM_ASM_SOURCE',CENSV_ROOT + '/software/svim-asm/src'))
    sys.path.insert(0,str(source))
    import svim_asm.SVIM_COMBINE as combine
    original=combine.write_final_vcf
    combine.write_final_vcf=lambda *a,**k:capture_write(original,*a,**k)
    program=source/'svim_asm/svim-asm'
    sys.argv[0]=str(program)
    runpy.run_path(str(program),run_name='__main__')


if __name__=='__main__':main()
