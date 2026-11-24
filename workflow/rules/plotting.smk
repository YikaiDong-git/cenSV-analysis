rule plot_variant_distribution:
    input:
        snp_variants = expand("snp_extraction/fksize{{fksize}}/chr{chr}_snp_variants.tsv", chr=[str(c) for c in range(1, 23)] + ["X"]),
        sv_summary = lambda wildcards: f"sv_extraction/{SV_MODE}/mapq{wildcards.mapq}_minsvlen{wildcards.minsvlen}_maxpart{wildcards.maxpart}_partdist{wildcards.partdist}_maxedit{wildcards.maxedit}_maxsv{wildcards.maxsv}_qgap{wildcards.qgap}_rgap{wildcards.rgap}/" + ("extraction_summary.tsv" if SV_MODE == "standard" else "allchr_sv_variants.tsv"),
        script = "workflow/scripts/plot/plot_variant_distribution.py"
    output:
        "plots/snp_distribution/fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/snp_counts.png",
        "plots/snp_distribution/fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/snp_density.png",
        "plots/snp_distribution/fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/snp_rarity_bar.png",
        "plots/snp_distribution/fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/snp_ac_cdf.png",
        "plots/snp_distribution/fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/sv_counts.png",
        "plots/snp_distribution/fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/sv_density.png",
        "plots/snp_distribution/fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/sv_rarity_bar.png",
        "plots/snp_distribution/fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/sv_ac_cdf.png",
        expand("plots/snp_distribution/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/snp_chr{chr}_correlation.png", chr=[str(c) for c in range(1, 23)] + ["X"]),
        expand("plots/snp_distribution/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/snp_chr{chr}_clustering.png", chr=[str(c) for c in range(1, 23)] + ["X"]),
        expand("plots/snp_distribution/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/sv_chr{chr}_correlation.png", chr=SV_CHRS),
        expand("plots/snp_distribution/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/sv_chr{chr}_clustering.png", chr=SV_CHRS)
    benchmark:
        "benchmarks/plot_variant_distribution/plot_variant_distribution_fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.tsv"
    params:
        snp_dir = "snp_extraction/fksize{fksize}",
        sv_dir = lambda wildcards: f"sv_extraction/{SV_MODE}/mapq{wildcards.mapq}_minsvlen{wildcards.minsvlen}_maxpart{wildcards.maxpart}_partdist{wildcards.partdist}_maxedit{wildcards.maxedit}_maxsv{wildcards.maxsv}_qgap{wildcards.qgap}_rgap{wildcards.rgap}",
        output_dir = "plots/snp_distribution/fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}"
    resources:
        mem_mb = 4000,
        runtime = 10
    shell:
        """
        python3 {input.script} \
            --snp-dir {params.snp_dir} \
            --sv-dir {params.sv_dir} \
            --output-dir {params.output_dir}
        """

rule EN_result_analysis:
    input:
        metrics = f"results/elastic_net/{SV_MODE}/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/chr{{chr}}/sv_metrics.tsv",
        coefs = f"results/elastic_net/{SV_MODE}/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/chr{{chr}}/model_coefficients.tsv",
        snp_vars = "snp_extraction/fksize{fksize}/chr{chr}_snp_variants.tsv",
        sv_vars = lambda wildcards: f"sv_extraction/{SV_MODE}/mapq{wildcards.mapq}_minsvlen{wildcards.minsvlen}_maxpart{wildcards.maxpart}_partdist{wildcards.partdist}_maxedit{wildcards.maxedit}_maxsv{wildcards.maxsv}_qgap{wildcards.qgap}_rgap{wildcards.rgap}/" + (f"chr{wildcards.chr}/chr{wildcards.chr}_sv_variants.tsv" if SV_MODE == "decompose" else f"chr{wildcards.chr}_sv_variants.tsv"),
        script = "workflow/scripts/plot/visualize_elastic_net.py"
    output:
        f"results/elastic_net/{SV_MODE}/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/chr{{chr}}/plots/{{chr}}_performance_distribution.png",
        f"results/elastic_net/{SV_MODE}/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/chr{{chr}}/plots/{{chr}}_distance_histogram.png"
    benchmark:
        "benchmarks/EN_result_analysis/EN_result_analysis_fksize{fksize}_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}_{chr}.tsv"
    params:
        output_dir = f"results/elastic_net/{SV_MODE}/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/chr{{chr}}/plots",
        prefix = "{chr}"
    resources:
        mem_mb = 4000,
        runtime = 10
    shell:
        """
        python3 {input.script} \
            --metrics {input.metrics} \
            --coefs {input.coefs} \
            --snp-variants {input.snp_vars} \
            --sv-variants {input.sv_vars} \
            --output-dir {params.output_dir} \
            --prefix {params.prefix}
        """

rule EN_global_analysis:
    input:
        metrics = expand(f"results/elastic_net/{SV_MODE}/fksize{{{{fksize}}}}_mapq{{{{mapq}}}}_minsvlen{{{{minsvlen}}}}_maxpart{{{{maxpart}}}}_partdist{{{{partdist}}}}_maxedit{{{{maxedit}}}}_maxsv{{{{maxsv}}}}_qgap{{{{qgap}}}}_rgap{{{{rgap}}}}/chr{{chr}}/sv_metrics.tsv", chr=SV_CHRS),
        coefs = expand(f"results/elastic_net/{SV_MODE}/fksize{{{{fksize}}}}_mapq{{{{mapq}}}}_minsvlen{{{{minsvlen}}}}_maxpart{{{{maxpart}}}}_partdist{{{{partdist}}}}_maxedit{{{{maxedit}}}}_maxsv{{{{maxsv}}}}_qgap{{{{qgap}}}}_rgap{{{{rgap}}}}/chr{{chr}}/model_coefficients.tsv", chr=SV_CHRS),
        snp_vars = expand("snp_extraction/fksize{{fksize}}/chr{chr}_snp_variants.tsv", chr=SV_CHRS),
        sv_vars = expand(f"sv_extraction/{SV_MODE}/mapq{{{{mapq}}}}_minsvlen{{{{minsvlen}}}}_maxpart{{{{maxpart}}}}_partdist{{{{partdist}}}}_maxedit{{{{maxedit}}}}_maxsv{{{{maxsv}}}}_qgap{{{{qgap}}}}_rgap{{{{rgap}}}}/" + ("chr{chr}/chr{chr}_sv_variants.tsv" if SV_MODE == "decompose" else "chr{chr}_sv_variants.tsv"), chr=SV_CHRS),
        script = "workflow/scripts/plot/visualize_elastic_net.py"
    output:
        f"results/elastic_net/{SV_MODE}/global_plots/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/global_distance_histogram.png",
        f"results/elastic_net/{SV_MODE}/global_plots/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/global_per_chr_distances.png"
    params:
        output_dir = f"results/elastic_net/{SV_MODE}/global_plots/fksize{{fksize}}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}",
        prefix = "global"
    resources:
        mem_mb = 8000,
        runtime = 10
    shell:
        """
        python3 {input.script} \
            --metrics {input.metrics} \
            --coefs {input.coefs} \
            --snp-variants {input.snp_vars} \
            --sv-variants {input.sv_vars} \
            --output-dir {params.output_dir} \
            --prefix {params.prefix}
        """
