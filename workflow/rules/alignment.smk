import os as _censv_os  # cenSV-analysis: cluster paths read from the environment
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
SCRATCH_DIR = config.get("scratch_dir", CENSV_ROOT + "/scratch/minimap2")

rule align_to_chm13:
    input:
        fasta=os.path.join("input", "HPRC", "{sample}_{hap}.fa.gz"),
        ref=os.path.join("resources", "chm13v2.0.fa.gz")
    output:
        "alignments/{mode}/{sample}_{hap}_to_chm13.bam"
    benchmark:
        "benchmarks/align_to_chm13/align_to_chm13_{mode}_{sample}_{hap}.tsv"
    resources:
        mem_mb=48000,
        runtime=120,
        cpus=4
    threads: 3
    shell:
        """
        scratch_file={SCRATCH_DIR}/{wildcards.sample}_{wildcards.hap}.fa.gz
        cp {input.fasta} $scratch_file # scratch has better IO

        mode="{wildcards.mode}"
        if [ "$mode" == "sniffles" ]; then
             minimap2 -x asm5 -a -t {threads} {input.ref} $scratch_file | samtools sort -@ {threads} -o {output} --write-index
        elif [ "$mode" == "svim_asm" ]; then
             minimap2 -a -x asm5 --cs -r2k -t {threads} {input.ref} $scratch_file | samtools sort -@ {threads} -o {output} --write-index
        elif [ "$mode" == "dipcall_sam" ]; then
             minimap2 -a -x asm5 --cs -t {threads} {input.ref} $scratch_file | gzip > {output}
        elif [ "$mode" == "dipcall_paf" ]; then
             minimap2 -c --paf-no-hit -x asm5 --cs -t {threads} {input.ref} $scratch_file | gzip > {output}
        fi
        """
