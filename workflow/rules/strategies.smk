FKSIZE_STRATEGIES = {"fksize100k": 100000, "fksize200k": 200000, "fksize500k": 500000}
EDGE_STRATEGIES   = {"edge500": 500, "edge1000": 1000, "edge2000": 2000}
MIDPOINT_STRATEGIES = {"midpoint500": 500, "midpoint1000": 1000}

ALL_STRATEGIES = list(FKSIZE_STRATEGIES.keys()) + list(EDGE_STRATEGIES.keys()) + list(MIDPOINT_STRATEGIES.keys())

STRAT_SVIM_SUFFIX = "mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}"

STRAT_EN_BASE  = "results/elastic_net/decompose/cosine_nj"
STRAT_IMP_BASE = "results/imputation/decompose/cosine_nj"

rule extract_snps_edge:
    input:
        vcf_dir = "resources/1000G",
        samples = "resources/overlap_samples.txt",
        centromere_bed = "resources/grch38_centromeres.bed",
        het_bed = "resources/grch38_heterochromatin_gaps.bed",
        script = "workflow/scripts/snp_extraction_edge.py"
    output:
        npz = "snp_extraction/edge{n_snps}/chr{chr}_snp_genotypes.npz",
        tsv = "snp_extraction/edge{n_snps}/chr{chr}_snp_variants.tsv",
        txt = "snp_extraction/edge{n_snps}/chr{chr}_snp_samples.txt"
    params:
        outdir = "snp_extraction/edge{n_snps}",
        vcf_pattern = lambda wildcards: "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    wildcard_constraints:
        n_snps = r"\d+"
    benchmark:
        "benchmarks/extract_snps_edge/extract_snps_edge_{n_snps}_chr{chr}.tsv"
    resources:
        mem_mb = 8000,
        runtime = 60,
        cpus = 1
    shell:
        """
        python3 {input.script} \
            --vcf-dir {input.vcf_dir} \
            --vcf-pattern "{params.vcf_pattern}" \
            --samples {input.samples} \
            --centromere-bed {input.centromere_bed} \
            --het-bed {input.het_bed} \
            --output {params.outdir} \
            --n-snps-per-side {wildcards.n_snps} \
            --chromosomes {wildcards.chr} \
            --min-maf 0.05
        """

rule extract_snps_midpoint:
    input:
        vcf_dir = "resources/1000G",
        samples = "resources/overlap_samples.txt",
        centromere_bed = "resources/grch38_centromeres.bed",
        script = "workflow/scripts/snp_extraction_midpoint.py"
    output:
        npz = "snp_extraction/midpoint{n_snps}/chr{chr}_snp_genotypes.npz",
        tsv = "snp_extraction/midpoint{n_snps}/chr{chr}_snp_variants.tsv",
        txt = "snp_extraction/midpoint{n_snps}/chr{chr}_snp_samples.txt"
    params:
        outdir = "snp_extraction/midpoint{n_snps}",
        vcf_pattern = lambda wildcards: "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    wildcard_constraints:
        n_snps = r"\d+"
    benchmark:
        "benchmarks/extract_snps_midpoint/extract_snps_midpoint_{n_snps}_chr{chr}.tsv"
    resources:
        mem_mb = 8000,
        runtime = 60,
        cpus = 1
    shell:
        """
        python3 {input.script} \
            --vcf-dir {input.vcf_dir} \
            --vcf-pattern "{params.vcf_pattern}" \
            --samples {input.samples} \
            --centromere-bed {input.centromere_bed} \
            --output {params.outdir} \
            --n-snps-per-side {wildcards.n_snps} \
            --chromosomes {wildcards.chr} \
            --min-maf 0.05
        """

FKSIZE_NAME_TO_VALUE = {"fksize100k": "100000", "fksize200k": "200000", "fksize500k": "500000"}

def get_strategy_snp_dir(wildcards):
    """Return the SNP extraction directory for a given strategy."""
    strat = wildcards.strategy
    if strat in FKSIZE_NAME_TO_VALUE:
        return f"snp_extraction/fksize{FKSIZE_NAME_TO_VALUE[strat]}"
    elif strat.startswith("edge"):
        n = strat.replace("edge", "")
        return f"snp_extraction/edge{n}"
    elif strat.startswith("midpoint"):
        n = strat.replace("midpoint", "")
        return f"snp_extraction/midpoint{n}"
    else:
        raise ValueError(f"Unknown strategy: {strat}")

def get_strategy_snp_npz(wildcards):
    """Return the SNP npz file path for a given strategy + chromosome."""
    snp_dir = get_strategy_snp_dir(wildcards)
    return f"{snp_dir}/chr{wildcards.chr}_snp_genotypes.npz"

def get_strategy_sv_dir(wildcards):
    """SV directory for decompose mode (per-chromosome subdirs)."""
    return (f"sv_extraction/decompose/mapq{wildcards.mapq}_minsvlen{wildcards.minsvlen}"
            f"_maxpart{wildcards.maxpart}_partdist{wildcards.partdist}"
            f"_maxedit{wildcards.maxedit}_maxsv{wildcards.maxsv}"
            f"_qgap{wildcards.qgap}_rgap{wildcards.rgap}/chr{wildcards.chr}")

rule elastic_net_train_strategy:
    input:
        snp_matrix = get_strategy_snp_npz,
        sv_summary = f"sv_extraction/decompose/{STRAT_SVIM_SUFFIX}/allchr_sv_variants.tsv",
        script = "workflow/scripts/elastic_net_train.py"
    output:
        summary = STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/sv_metrics.tsv",
        trained_models = STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/trained_models.pkl",
        model_coefficients = STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/model_coefficients.tsv"
    params:
        snp_dir = get_strategy_snp_dir,
        sv_dir  = get_strategy_sv_dir,
        outdir  = STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX,
        min_snp_freq = 0.05,
        min_sv_freq = 0.05,
        max_snp_freq = 0.99,
        max_sv_freq = 0.99,
        sv_filter = "all"
    wildcard_constraints:
        strategy = "|".join(ALL_STRATEGIES)
    benchmark:
        "benchmarks/elastic_net_strategy/en_{strategy}_" + STRAT_SVIM_SUFFIX + "_chr{chr}.tsv"
    resources:
        mem_mb = 16000,
        runtime = 1440,
        cpus = 12
    threads: 12
    shell:
        """
        python3 {input.script} \
            --snp-dir {params.snp_dir} \
            --sv-dir {params.sv_dir} \
            --output {params.outdir} \
            --chromosome chr{wildcards.chr} \
            --n-folds 4 \
            --n-jobs {threads} \
            --min-snp-freq {params.min_snp_freq} \
            --min-sv-freq {params.min_sv_freq} \
            --max-snp-freq {params.max_snp_freq} \
            --max-sv-freq {params.max_sv_freq} \
            --sv-filter {params.sv_filter}
        """

rule elastic_net_train_chunk:
    input:
        snp_matrix = get_strategy_snp_npz,
        sv_summary = f"sv_extraction/decompose/{STRAT_SVIM_SUFFIX}/allchr_sv_variants.tsv",
        script = "workflow/scripts/elastic_net_train_chunk.py"
    output:
        summary = STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chunk_{chunk_id}/sv_metrics.tsv",
        trained_models = STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chunk_{chunk_id}/trained_models.pkl"
    params:
        snp_dir = get_strategy_snp_dir,
        sv_dir  = get_strategy_sv_dir,
        outdir  = STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chunk_{chunk_id}"
    wildcard_constraints:
        strategy = "|".join(ALL_STRATEGIES),
        chunk_id = r"\d+"
    benchmark:
        "benchmarks/elastic_net_chunk/en_{strategy}_" + STRAT_SVIM_SUFFIX + "_chr{chr}_chunk{chunk_id}.tsv"
    resources:
        mem_mb = 16000,
        runtime = 1440,
        cpus = 12
    threads: 12
    shell:
        """
        python3 {input.script} \
            --snp-dir {params.snp_dir} \
            --sv-dir {params.sv_dir} \
            --output {params.outdir} \
            --chromosome chr{wildcards.chr} \
            --sv-start $(( {wildcards.chunk_id} * 1000 )) \
            --sv-end $(( ({wildcards.chunk_id} + 1) * 1000 )) \
            --n-folds 4 \
            --n-jobs {threads} \
            --min-snp-freq 0.05 \
            --min-sv-freq 0.05 \
            --sv-filter all
        """

rule impute_1kg_strategy:
    input:
        en_metrics = STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/sv_metrics.tsv",
        vcf = "resources/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
        samples = IMPUTATION_SAMPLES,
        script = "workflow/scripts/impute_1kg.py"
    output:
        dosages = STRAT_IMP_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chr{chr}_imputed_dosages.tsv",
        summary = STRAT_IMP_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chr{chr}_imputation_summary.tsv",
        snp_matching = STRAT_IMP_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chr{chr}_snp_matching.tsv",
        sample_list = STRAT_IMP_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chr{chr}_imputed_samples.txt"
    params:
        model_dir = STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}"
    wildcard_constraints:
        strategy = "|".join(ALL_STRATEGIES)
    benchmark:
        "benchmarks/impute_strategy/impute_{strategy}_" + STRAT_SVIM_SUFFIX + "_chr{chr}.tsv"
    resources:
        mem_mb = 16000,
        runtime = 60,
        cpus = 1
    shell:
        """
        python3 {input.script} \
            --model-dir {params.model_dir} \
            --vcf {input.vcf} \
            --samples {input.samples} \
            --output $(dirname {output.dosages}) \
            --chromosome chr{wildcards.chr} \
            --r2-threshold {IMPUTATION_R2_THRESHOLD}
        """

rule select_best_strategy:
    input:
        en_metrics = expand(
            STRAT_EN_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/sv_metrics.tsv",
            strategy=ALL_STRATEGIES, chr=SV_CHRS,
            mapq=MIN_MAPQ_SVIM, minsvlen=MINSVLEN, maxpart=MAX_PARTITION_SIZE,
            partdist=PARTITION_MAX_DISTANCE, maxedit=MAX_EDIT_DISTANCE,
            maxsv=MAX_SV_SIZE, qgap=QUERY_GAP_TOLERANCE, rgap=REFERENCE_GAP_TOLERANCE),
        imp_dosages = expand(
            STRAT_IMP_BASE + "/{strategy}_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chr{chr}_imputed_dosages.tsv",
            strategy=ALL_STRATEGIES, chr=SV_CHRS,
            mapq=MIN_MAPQ_SVIM, minsvlen=MINSVLEN, maxpart=MAX_PARTITION_SIZE,
            partdist=PARTITION_MAX_DISTANCE, maxedit=MAX_EDIT_DISTANCE,
            maxsv=MAX_SV_SIZE, qgap=QUERY_GAP_TOLERANCE, rgap=REFERENCE_GAP_TOLERANCE),
        script = "workflow/scripts/select_best_strategy.py"
    output:
        best8_flag = expand(
            STRAT_IMP_BASE + "/best8_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chr{chr}_imputed_dosages.tsv",
            chr=SV_CHRS,
            mapq=MIN_MAPQ_SVIM, minsvlen=MINSVLEN, maxpart=MAX_PARTITION_SIZE,
            partdist=PARTITION_MAX_DISTANCE, maxedit=MAX_EDIT_DISTANCE,
            maxsv=MAX_SV_SIZE, qgap=QUERY_GAP_TOLERANCE, rgap=REFERENCE_GAP_TOLERANCE),
        best6_flag = expand(
            STRAT_IMP_BASE + "/best6_" + STRAT_SVIM_SUFFIX + "/chr{chr}/chr{chr}_imputed_dosages.tsv",
            chr=SV_CHRS,
            mapq=MIN_MAPQ_SVIM, minsvlen=MINSVLEN, maxpart=MAX_PARTITION_SIZE,
            partdist=PARTITION_MAX_DISTANCE, maxedit=MAX_EDIT_DISTANCE,
            maxsv=MAX_SV_SIZE, qgap=QUERY_GAP_TOLERANCE, rgap=REFERENCE_GAP_TOLERANCE)
    params:
        en_base = STRAT_EN_BASE,
        imp_base = STRAT_IMP_BASE,
        suffix = lambda wc: f"_mapq{MIN_MAPQ_SVIM[0]}_minsvlen{MINSVLEN[0]}_maxpart{MAX_PARTITION_SIZE[0]}_partdist{PARTITION_MAX_DISTANCE[0]}_maxedit{MAX_EDIT_DISTANCE[0]}_maxsv{MAX_SV_SIZE[0]}_qgap{QUERY_GAP_TOLERANCE[0]}_rgap{REFERENCE_GAP_TOLERANCE[0]}",
        chroms = ",".join(["chr" + c for c in SV_CHRS])
    benchmark:
        "benchmarks/select_best_strategy/select_best_strategy.tsv"
    resources:
        mem_mb = 16000,
        runtime = 60,
        cpus = 1
    shell:
        """
        python3 {input.script} \
            --en-base-dir {params.en_base} \
            --imp-base-dir {params.imp_base} \
            --output {params.imp_base} \
            --suffix "{params.suffix}" \
            --chromosomes {params.chroms} \
            --mode both
        """
