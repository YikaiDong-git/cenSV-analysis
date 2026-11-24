rule all_short_reads:
    """Aggregate target for short-reads related outputs (peri-centromeric SNPs)."""
    input:
        expand("resources/1000G/peri_cen_snps.chr{chr}.vcf.gz", chr=[str(c) for c in range(1, 23)] + ["X"])

def get_1000g_vcf_url(wildcards):
    base_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
    if wildcards.chr == "X":
        return f"{base_url}/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    else:
        return f"{base_url}/1kGP_high_coverage_Illumina.chr{wildcards.chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

rule download_1000g_phased:
    output:
        vcf = "resources/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
        tbi = "resources/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"
    params:
        vcf_url = get_1000g_vcf_url,
        tbi_url = lambda w: get_1000g_vcf_url(w) + ".tbi"
    benchmark:
        "benchmarks/download_1000g_phased/download_1000g_phased_{chr}.tsv"
    resources:
        mem_mb=2000,
        runtime=120,
        cpus=1
    shell:
        """
        wget {params.vcf_url} -O {output.vcf}
        wget {params.tbi_url} -O {output.tbi}
        """

rule generate_pericentromeric_bed:
    input:
        centromeres = "resources/grch38_centromeres.bed",
        het_gaps = "resources/grch38_heterochromatin_gaps.bed"
    benchmark:
        "benchmarks/generate_pericentromeric_bed/generate_pericentromeric_bed_{fksize}.tsv"
    output:
        "resources/peri_centromeric_regions_grch38_{fksize}.bed"
    params:
        script = "workflow/scripts/generate_pericen_bed.py"
    shell:
        """
        python3 {params.script} \
            --centromere-bed {input.centromeres} \
            --het-bed {input.het_gaps} \
            --fksize {wildcards.fksize} \
            --output {output}
        """

rule extract_pericentromeric_snps:
    input:
        vcf = "resources/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
        tbi = "resources/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi",
        regions = "resources/peri_centromeric_regions_grch38_{fksize}.bed",
        samples = "resources/overlap_samples.txt"
    benchmark:
        "benchmarks/extract_pericentromeric_snps/extract_pericentromeric_snps_{chr}_{fksize}.tsv"
    output:
        vcf = "resources/1000G/fksize{fksize}_peri_cen_snps.chr{chr}.vcf.gz",
        tbi = "resources/1000G/fksize{fksize}_peri_cen_snps.chr{chr}.vcf.gz.tbi"
    shell:
        """
                bcftools view -R {input.regions} -S {input.samples} -Ou {input.vcf} | \
                bcftools view -v snps -m2 -M2 --min-af 0.05:minor -Oz -o {output.vcf} - && \
                bcftools index -t {output.vcf}
        """

rule extract_snps_python:
    input:
        vcf = "resources/1000G/fksize{fksize}_peri_cen_snps.chr{chr}.vcf.gz",
        tbi = "resources/1000G/fksize{fksize}_peri_cen_snps.chr{chr}.vcf.gz.tbi",
        samples = "resources/overlap_samples.txt",
        script = "workflow/scripts/snp_extraction.py"
    output:
        npz = "snp_extraction/fksize{fksize}/chr{chr}_snp_genotypes.npz",
        tsv = "snp_extraction/fksize{fksize}/chr{chr}_snp_variants.tsv",
        txt = "snp_extraction/fksize{fksize}/chr{chr}_snp_samples.txt"
    benchmark:
        "benchmarks/extract_snps_python/extract_snps_python_{chr}_{fksize}.tsv"
    shell:
        """
        python3 {input.script} \
            --vcf-dir resources/1000G \
            --vcf-pattern "fksize{wildcards.fksize}_peri_cen_snps.chr{{chrom}}.vcf.gz" \
            --samples {input.samples} \
            --output snp_extraction/fksize{wildcards.fksize} \
            --chromosomes {wildcards.chr}
        """
