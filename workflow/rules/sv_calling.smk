rule filter_to_active:
    input:
        bam="alignments/sniffles/{sample}_{hap}_to_chm13.bam",
        bed=os.path.join("resources", "chm13v2.0_censat_v2.1.active.filtered.bed")
    output:
        "alignments/sniffles/{sample}_{hap}_to_chm13.active.bam"
    benchmark:
        "benchmarks/filter_to_active/filter_to_active_{sample}_{hap}.tsv"
    resources:
        mem_mb=16000,
        runtime=10,
        cpus=2
    shell:
        """
        samtools view --region-file {input.bed} -b {input.bam} | samtools sort -@ {threads} -o {output} --write-index
        """

def get_tr_bed_input(wildcards):
    if wildcards.tr_mode == "censat_filtered":
        return "resources/chm13v2.0_censat_v2.1.active.filtered.bed"
    elif wildcards.tr_mode == "censat":
        return "resources/chm13v2.0_censat_v2.1.bed"
    return []

rule sniffles_candidate:
    input:
        bam = "alignments/sniffles/{sample}_{hap}_to_chm13.active.bam",
        ref = rules.download_reference.output.fasta,
        tr_bed = get_tr_bed_input
    output:
        "variant_calling/sniffles/{sample}_{hap}_to_chm13.active_trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}.snf"
    benchmark:
        "benchmarks/sniffles_candidate/sniffles_candidate_{sample}_{hap}_trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}.tsv"
    resources:
        mem_mb=4000,
        runtime=60,
        cpus=2
    threads:
        2
    params:
        tr_arg = lambda wildcards: f"--tandem-repeats {get_tr_bed_input(wildcards)}" if get_tr_bed_input(wildcards) else ""
    shell:
        """
        sniffles --input {input.bam} --reference {input.ref} --threads {threads} --snf {output} --mapq {wildcards.mapq} --minsvlen {wildcards.minsvlen} --minsupport {wildcards.minsupport} {params.tr_arg}
        """

def get_snf_inputs(wildcards):
    return expand("variant_calling/sniffles/{sample}_{hap}_to_chm13.active_trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}.snf",
                  sample=SAMPLES, hap=HAPLOTYPES,
                  mapq=[wildcards.mapq], minsvlen=[wildcards.minsvlen],
                  minsupport=[wildcards.minsupport], tr_mode=[wildcards.tr_mode])

rule sniffles_call:
    input:
        snf = get_snf_inputs,
        ref = rules.download_reference.output.fasta
    output:
        "variant_calling/sniffles/cohort_sniffles_trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}.vcf.gz"
    benchmark:
        "benchmarks/sniffles_call/sniffles_call_trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}.tsv"
    resources:
        mem_mb = 1000,
        runtime = 10,
        cpus = 2
    threads: 2
    shell:
        """
        sniffles --input {input.snf} --threads {threads} --reference {input.ref} \
            --vcf {output} --combine-max-inmemory-results 600
        """

rule svim_asm_call_filter_zip:
    wildcard_constraints:
        hap="hap1|hap2"
    input:
        bam = "alignments/svim_asm/{sample}_{hap}_to_chm13.bam",
        ref = rules.download_reference.output.fasta,
        bed = "resources/chm13v2.0_censat_v2.1.active.bed"
    output:
        vcf = "variant_calling/svim_asm/{sample}_{hap}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.vcf.gz",
        tbi = "variant_calling/svim_asm/{sample}_{hap}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.vcf.gz.tbi"
    params:
        sample_name = "{sample}_{hap}",
        outdir = "variant_calling/svim_asm/{sample}_{hap}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}_tmp"
    benchmark:
        "benchmarks/svim_asm_call_filter_zip/svim_asm_call_{sample}_{hap}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.tsv"
    resources:
        mem_mb=16000,
        runtime=60,
        cpus=2
    shell:
        """
        export PYTHONPATH=software/svim-asm/src${{PYTHONPATH:+:$PYTHONPATH}}
        python3 software/svim-asm/src/svim_asm/svim-asm haploid {params.outdir} {input.bam} {input.ref} \
            --sample {params.sample_name} \
            --min_mapq {wildcards.mapq} \
            --min_sv_size {wildcards.minsvlen} \
            --max_partition_size {wildcards.maxpart} \
            --partition_max_distance {wildcards.partdist} \
            --max_edit_distance {wildcards.maxedit} \
            --max_sv_size {wildcards.maxsv} \
            --query_gap_tolerance {wildcards.qgap} \
            --reference_gap_tolerance {wildcards.rgap}

        bcftools view -T {input.bed} -O z -o {output.vcf} \
            {params.outdir}/variants.vcf

        tabix -p vcf {output.vcf}

        rm -rf {params.outdir}
        """

rule svim_asm_merge:
    input:
        vcfs = lambda wildcards: expand("variant_calling/svim_asm/{sample}_{hap}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.vcf.gz",
                      sample=SAMPLES, hap=HAPLOTYPES, mapq=[wildcards.mapq], minsvlen=[wildcards.minsvlen], maxpart=[wildcards.maxpart],
                      partdist=[wildcards.partdist], maxedit=[wildcards.maxedit], maxsv=[wildcards.maxsv], qgap=[wildcards.qgap], rgap=[wildcards.rgap]),
        tbis = lambda wildcards: expand("variant_calling/svim_asm/{sample}_{hap}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.vcf.gz.tbi",
                      sample=SAMPLES, hap=HAPLOTYPES, mapq=[wildcards.mapq], minsvlen=[wildcards.minsvlen], maxpart=[wildcards.maxpart],
                      partdist=[wildcards.partdist], maxedit=[wildcards.maxedit], maxsv=[wildcards.maxsv], qgap=[wildcards.qgap], rgap=[wildcards.rgap])
    output:
        vcf = "variant_calling/svim_asm/cohort_svim_asm_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.vcf.gz",
        tbi = "variant_calling/svim_asm/cohort_svim_asm_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.vcf.gz.tbi"
    benchmark:
        "benchmarks/svim_asm_merge/svim_asm_merge_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.tsv"
    resources:
        mem_mb=16000,
        runtime=60
    shell:
        """
        bcftools merge {input.vcfs} -O z -o {output.vcf}
        tabix -p vcf {output.vcf}
        """
