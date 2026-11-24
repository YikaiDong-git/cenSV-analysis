rule download_reference:
    output:
        fasta = os.path.join("resources", "chm13v2.0.fa.gz"),
        gzi = os.path.join("resources", "chm13v2.0.fa.gz.gzi")
    benchmark:
        "benchmarks/download_reference/download_reference.tsv"
    shell:
        """
        aws s3 cp --no-sign-request s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz {output.fasta}
        aws s3 cp --no-sign-request s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz.gzi {output.gzi}
        """

rule download_censat:
    output:
        all_bed =      os.path.join("resources", "chm13v2.0_censat_v2.1.bed"),
        active_bed = os.path.join("resources", "chm13v2.0_censat_v2.1.active.bed")
    params:
        pattern = "H1L"
    benchmark:
        "benchmarks/download_censat/download_censat.tsv"
    shell:
        """
        aws s3 cp --no-sign-request s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed {output.all_bed}
        grep {params.pattern} {output.all_bed} > {output.active_bed}
        """

rule filter_active_bed:
    input:
        os.path.join("resources", "chm13v2.0_censat_v2.1.active.bed")
    output:
        os.path.join("resources", "chm13v2.0_censat_v2.1.active.filtered.bed")
    benchmark:
        "benchmarks/filter_active_bed/filter_active_bed.tsv"
    shell:
        r"""
        awk -F'\t' '{{ if(match($4,/\(([^)]*)\)/,m) && m[1] ~ /\// && m[1] ~ /S[0-9]+C/){{ next }} print }}' {input} > {output}
        """

def get_assembly_url(wildcards):
    if wildcards.sample not in SAMPLES_INFO:
        return f"s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/{wildcards.sample}/assemblies/freeze_2/{wildcards.sample}_{wildcards.hap}_hprc_r2_v1.fa.gz" # v1.0

    info = SAMPLES_INFO[wildcards.sample]
    group = info.get('group')
    version = info.get('version')
    hap_type_group1_2 = info.get('hap_type_group1_2', '')

    base_url = "s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION"

    if group == "Group 3":
        remote_hap = wildcards.hap
        return f"{base_url}/{wildcards.sample}/assemblies/freeze_2/assembly_pipeline/ncbi_upload/{wildcards.sample}.{remote_hap}_{version}.fa.gz"
    else:
        if hap_type_group1_2 == 'mat_pat':
            remote_hap = 'mat' if wildcards.hap == 'hap1' else 'pat'
        else:
            remote_hap = wildcards.hap

        return f"{base_url}/{wildcards.sample}/assemblies/freeze_2/{wildcards.sample}_{remote_hap}_{version}.fa.gz"

def get_rerun_assembly_url(wildcards):
    url = get_assembly_url(wildcards)
    if 'assembly_pipeline/ncbi_upload' in url:
        return url.replace('assembly_pipeline/ncbi_upload', 'assembly_pipeline/rerun/ncbi_upload')
    return ""

rule prepare_assembly_url_soup_and_classify_groups:
    input:
        script="workflow/scripts/helper/classify_samples.py"
    output:
        csv=classification_csv
    benchmark:
        "benchmarks/prepare_assembly_url_soup_and_classify_groups/prepare_assembly_url_soup_and_classify_groups.tsv"
    shell:
        """
        python3 {input.script}
        """

rule download_assemblies:
    output:
        os.path.join("input", "HPRC", "{sample}_{hap}.fa.gz")
    params:
        url = get_assembly_url,
        url_rerun = get_rerun_assembly_url
    benchmark:
        "benchmarks/download_assemblies/download_assemblies_{sample}_{hap}.tsv"
    shell:
        """
        aws s3 cp --no-sign-request {params.url} {output} || ( [ -n "{params.url_rerun}" ] && aws s3 cp --no-sign-request {params.url_rerun} {output} ) || exit 1
        """

rule download_1000g_metadata:
    output:
        "resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
    benchmark:
        "benchmarks/download_1000g_metadata/download_1000g_metadata.tsv"
    shell:
        """
        wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt -O {output}
        """

rule identify_overlap_samples:
    input:
        my_samples = "workflow/scripts/helper/sample_classification.csv",
        g1k_samples = "resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
    output:
        "resources/overlap_samples.txt"
    params:
        groups = ",".join(TARGET_GROUPS)
    benchmark:
        "benchmarks/identify_overlap_samples/identify_overlap_samples.tsv"
    shell:
        """
        python3 workflow/scripts/helper/identify_overlap_samples.py \
            --my_samples {input.my_samples} \
            --g1k_samples {input.g1k_samples} \
            --output {output} \
            --groups "{params.groups}"
        """

rule create_centromeres_bed:
    """UCSC hg38 cytoBand `acen` bands merged per chromosome, with the four recorded
    adjustments applied and chrY dropped.

    The inline awk this replaces emitted the plain merge, so the rule did not reproduce its
    own committed output: a clean re-run moved the centromere boundary on chr1, chr2, chr6
    and chrX, and with it every Strategy 1 window. workflow/scripts/make_centromere_bed.py carries
    the four adjustments explicitly and asserts the acen merge each was derived from, so a
    change upstream at UCSC fails the build rather than being applied silently.

    Verify the committed resource still regenerates:
        python3 workflow/scripts/make_centromere_bed.py --check resources/grch38_centromeres.bed
    """
    output:
        "resources/grch38_centromeres.bed"
    params:
        script = "workflow/scripts/make_centromere_bed.py"
    benchmark:
        "benchmarks/create_centromeres_bed/create_centromeres_bed.tsv"
    shell:
        """
        python3 {params.script} --out {output}
        """

rule create_heterochromatin_gaps_bed:
    """Download heterochromatin and short_arm gap annotations from UCSC hg38 gap table."""
    output:
        "resources/grch38_heterochromatin_gaps.bed"
    benchmark:
        "benchmarks/create_heterochromatin_gaps_bed/create_heterochromatin_gaps_bed.tsv"
    shell:
        """
        curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz" | \
          gunzip -c | \
          awk -F'\\t' '($8=="heterochromatin" || $8=="short_arm") {{print $2"\\t"$3"\\t"$4"\\t"$8}}' | \
          sort -k1,1V -k2,2n > {output}
        """
