rule extract_svs_python:
    input:
        vcf = "variant_calling/svim_asm/cohort_svim_asm_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.vcf.gz",
        samples = "resources/overlap_samples.txt",
        script = "workflow/scripts/sv_extraction.py"
    benchmark:
        "benchmarks/extract_svs_python/extract_svs_python_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.tsv"
    output:
        summary = f"sv_extraction/standard/mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/extraction_summary.tsv",
        npzs = expand(f"sv_extraction/standard/mapq{{{{mapq}}}}_minsvlen{{{{minsvlen}}}}_maxpart{{{{maxpart}}}}_partdist{{{{partdist}}}}_maxedit{{{{maxedit}}}}_maxsv{{{{maxsv}}}}_qgap{{{{qgap}}}}_rgap{{{{rgap}}}}/chr{{chr}}_sv_genotypes.npz", chr=SV_CHRS),
        variants = expand(f"sv_extraction/standard/mapq{{{{mapq}}}}_minsvlen{{{{minsvlen}}}}_maxpart{{{{maxpart}}}}_partdist{{{{partdist}}}}_maxedit{{{{maxedit}}}}_maxsv{{{{maxsv}}}}_qgap{{{{qgap}}}}_rgap{{{{rgap}}}}/chr{{chr}}_sv_variants.tsv", chr=SV_CHRS),
        samples = expand(f"sv_extraction/standard/mapq{{{{mapq}}}}_minsvlen{{{{minsvlen}}}}_maxpart{{{{maxpart}}}}_partdist{{{{partdist}}}}_maxedit{{{{maxedit}}}}_maxsv{{{{maxsv}}}}_qgap{{{{qgap}}}}_rgap{{{{rgap}}}}/chr{{chr}}_sv_samples.txt", chr=SV_CHRS)
    params:
        outdir = f"sv_extraction/standard/mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}",
        chroms = ",".join(["chr" + c for c in SV_CHRS])
    shell:
        """
        python3 {input.script} \
            --vcf {input.vcf} \
            --samples {input.samples} \
            --output {params.outdir} \
            --chromosomes {params.chroms}
        """

rule extract_svs_python_all:
    input:
        expand(
            "sv_extraction/standard/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/extraction_summary.tsv",
            mapq=MIN_MAPQ_SVIM, minsvlen=MINSVLEN, maxpart=MAX_PARTITION_SIZE,
            partdist=PARTITION_MAX_DISTANCE, maxedit=MAX_EDIT_DISTANCE, maxsv=MAX_SV_SIZE,
            qgap=QUERY_GAP_TOLERANCE, rgap=REFERENCE_GAP_TOLERANCE
        )

rule elastic_net_train:
    input:
        snp_matrix = "snp_extraction/fksize{fksize}/chr{chr}_snp_genotypes.npz",
        sv_summary = lambda wildcards: f"sv_extraction/{SV_MODE}/mapq{wildcards.mapq}_minsvlen{wildcards.minsvlen}_maxpart{wildcards.maxpart}_partdist{wildcards.partdist}_maxedit{wildcards.maxedit}_maxsv{wildcards.maxsv}_qgap{wildcards.qgap}_rgap{wildcards.rgap}/" + ("extraction_summary.tsv" if SV_MODE == "standard" else "allchr_sv_variants.tsv"),
        script  = "workflow/scripts/elastic_net_train.py"
    output:
        summary = f"results/elastic_net/{SV_MODE}/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/chr{{chr}}/sv_metrics.tsv",
        trained_models = f"results/elastic_net/{SV_MODE}/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/chr{{chr}}/trained_models.pkl",
        model_coefficients = f"results/elastic_net/{SV_MODE}/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/chr{{chr}}/model_coefficients.tsv"
    benchmark:
        "benchmarks/elastic_net_train/elastic_net_train_fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}_chr{chr}.tsv"
    params:
        snp_dir = "snp_extraction/fksize{fksize}",
        sv_dir  = lambda wildcards: f"sv_extraction/{SV_MODE}/mapq{wildcards.mapq}_minsvlen{wildcards.minsvlen}_maxpart{wildcards.maxpart}_partdist{wildcards.partdist}_maxedit{wildcards.maxedit}_maxsv{wildcards.maxsv}_qgap{wildcards.qgap}_rgap{wildcards.rgap}" + (f"/chr{wildcards.chr}" if SV_MODE == "decompose" else ""),
        outdir  = f"results/elastic_net/{SV_MODE}/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}",
        min_snp_freq = 0.05,
        min_sv_freq = 0.05,
        max_snp_freq = 0.99,
        max_sv_freq = 0.99,
        sv_filter = "all"
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
