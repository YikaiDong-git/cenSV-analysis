rule haplogroup_clustering_sniffles:
    input:
        vcf="variant_calling/sniffles/cohort_sniffles_trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}.vcf.gz"
    output:
        hg="results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}/chr{chrom}/chr{chrom}_haplogroups.tsv",
        var="results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}/chr{chrom}/chr{chrom}_variants.tsv",
        decomp="results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}/chr{chrom}/chr{chrom}_decomposition.tsv",
        val="results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}/chr{chrom}/chr{chrom}_validated_haplogroups.tsv"
    params:
        outdir="results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}",
        script="workflow/scripts/haplogroup_clustering.py",
        k_range="2,30",
        min_samples=2
    log:
        "results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}/chr{chrom}/run.log"
    benchmark:
        "benchmarks/haplogroup_clustering_sniffles/haplogroup_clustering_trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}_chr{chrom}.tsv"
    resources:
        mem_mb=8000,
        runtime=10,
        cpus=1
    shell:
        """
        python3 {params.script} \
            --vcf {input.vcf} \
            --output {params.outdir} \
            --chrom chr{wildcards.chrom} \
            --min-samples {params.min_samples} \
            --k-range {params.k_range} \
            --decompose \
            2> {log} 1> {log}
        """

rule merge_sniffles_outputs:
    input:
        hgs=expand("results/haplogroup_clustering/sniffles/trf_{{tr_mode}}_minsupp{{minsupport}}_mapq{{mapq}}_minsvlen{{minsvlen}}/chr{chrom}/chr{chrom}_haplogroups.tsv", chrom=SV_CHRS),
        vals=expand("results/haplogroup_clustering/sniffles/trf_{{tr_mode}}_minsupp{{minsupport}}_mapq{{mapq}}_minsvlen{{minsvlen}}/chr{chrom}/chr{chrom}_validated_haplogroups.tsv", chrom=SV_CHRS),
        vars=expand("results/haplogroup_clustering/sniffles/trf_{{tr_mode}}_minsupp{{minsupport}}_mapq{{mapq}}_minsvlen{{minsvlen}}/chr{chrom}/chr{chrom}_variants.tsv", chrom=SV_CHRS),
        decomps=expand("results/haplogroup_clustering/sniffles/trf_{{tr_mode}}_minsupp{{minsupport}}_mapq{{mapq}}_minsvlen{{minsvlen}}/chr{chrom}/chr{chrom}_decomposition.tsv", chrom=SV_CHRS)
    output:
        all_hg="results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}/allchr_haplogroups.tsv",
        all_val="results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}/allchr_validated_haplogroups.tsv",
        all_var="results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}/allchr_variants.tsv",
        all_decomp="results/haplogroup_clustering/sniffles/trf_{tr_mode}_minsupp{minsupport}_mapq{mapq}_minsvlen{minsvlen}/allchr_decomposition.tsv"
    shell:
        """
        head -n 1 {input.hgs[0]} > {output.all_hg}
        for f in {input.hgs}; do tail -n +2 $f >> {output.all_hg}; done
        head -n 1 {input.vals[0]} > {output.all_val}
        for f in {input.vals}; do tail -n +2 $f >> {output.all_val}; done
        head -n 1 {input.vars[0]} > {output.all_var}
        for f in {input.vars}; do tail -n +2 $f >> {output.all_var}; done
        head -n 1 {input.decomps[0]} > {output.all_decomp}
        for f in {input.decomps}; do tail -n +2 $f >> {output.all_decomp}; done
        """

rule haplogroup_clustering_svim:
    input:
        vcf = "variant_calling/svim_asm/cohort_svim_asm_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}.vcf.gz",
        script = "workflow/scripts/haplogroup_clustering.py"
    output:
        hg=SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/chr{chrom}/chr{chrom}_haplogroups.tsv",
        var=SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/chr{chrom}/chr{chrom}_sv_variants.tsv",
        decomp=SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/chr{chrom}/chr{chrom}_decomposition.tsv",
        val=SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/chr{chrom}/chr{chrom}_validated_haplogroups.tsv"
    params:
        outdir=SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}",
        k_range="2,30",
        min_samples=2
    log:
        SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/chr{chrom}/run.log"
    benchmark:
        "benchmarks/haplogroup_clustering_svim/haplogroup_clustering_svim_asm_mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}_chr{chrom}.tsv"
    resources:
        mem_mb=8000,
        runtime=10,
        cpus=1
    shell:
        """
        python3 {input.script} \
            --vcf {input.vcf} \
            --output {params.outdir} \
            --chrom chr{wildcards.chrom} \
            --min-samples {params.min_samples} \
            --k-range {params.k_range} \
            --decompose \
            2> {log} 1> {log}
        """

rule merge_svim_outputs:
    input:
        hgs=lambda w: expand(SV_DIR + "/mapq" + str(w.mapq) + "_minsvlen" + str(w.minsvlen) + "_maxpart" + str(w.maxpart) + "_partdist" + str(w.partdist) + "_maxedit" + str(w.maxedit) + "_maxsv" + str(w.maxsv) + "_qgap" + str(w.qgap) + "_rgap" + str(w.rgap) + "/chr{chrom}/chr{chrom}_haplogroups.tsv", chrom=SV_CHRS),
        vals=lambda w: expand(SV_DIR + "/mapq" + str(w.mapq) + "_minsvlen" + str(w.minsvlen) + "_maxpart" + str(w.maxpart) + "_partdist" + str(w.partdist) + "_maxedit" + str(w.maxedit) + "_maxsv" + str(w.maxsv) + "_qgap" + str(w.qgap) + "_rgap" + str(w.rgap) + "/chr{chrom}/chr{chrom}_validated_haplogroups.tsv", chrom=SV_CHRS),
        vars=lambda w: expand(SV_DIR + "/mapq" + str(w.mapq) + "_minsvlen" + str(w.minsvlen) + "_maxpart" + str(w.maxpart) + "_partdist" + str(w.partdist) + "_maxedit" + str(w.maxedit) + "_maxsv" + str(w.maxsv) + "_qgap" + str(w.qgap) + "_rgap" + str(w.rgap) + "/chr{chrom}/chr{chrom}_sv_variants.tsv", chrom=SV_CHRS),
        decomps=lambda w: expand(SV_DIR + "/mapq" + str(w.mapq) + "_minsvlen" + str(w.minsvlen) + "_maxpart" + str(w.maxpart) + "_partdist" + str(w.partdist) + "_maxedit" + str(w.maxedit) + "_maxsv" + str(w.maxsv) + "_qgap" + str(w.qgap) + "_rgap" + str(w.rgap) + "/chr{chrom}/chr{chrom}_decomposition.tsv", chrom=SV_CHRS)
    output:
        all_hg=SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/allchr_haplogroups.tsv",
        all_val=SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/allchr_validated_haplogroups.tsv",
        all_var=SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/allchr_sv_variants.tsv",
        all_decomp=SV_DIR + "/mapq{mapq}_minsvlen{minsvlen}_maxpart{maxpart}_partdist{partdist}_maxedit{maxedit}_maxsv{maxsv}_qgap{qgap}_rgap{rgap}/allchr_decomposition.tsv"
    shell:
        """
        head -n 1 {input.hgs[0]} > {output.all_hg}
        for f in {input.hgs}; do tail -n +2 $f >> {output.all_hg}; done
        head -n 1 {input.vals[0]} > {output.all_val}
        for f in {input.vals}; do tail -n +2 $f >> {output.all_val}; done
        head -n 1 {input.vars[0]} > {output.all_var}
        for f in {input.vars}; do tail -n +2 $f >> {output.all_var}; done
        head -n 1 {input.decomps[0]} > {output.all_decomp}
        for f in {input.decomps}; do tail -n +2 $f >> {output.all_decomp}; done
        """

rule analyze_clustering_params:
    input:
        tsv = expand(f"sv_extraction/{SV_MODE}/mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}/allchr_haplogroups.tsv",
                      mapq=MIN_MAPQ_SVIM, minsvlen=MINSVLEN, maxpart=MAX_PARTITION_SIZE, partdist=PARTITION_MAX_DISTANCE, maxedit=MAX_EDIT_DISTANCE, maxsv=MAX_SV_SIZE, qgap=QUERY_GAP_TOLERANCE, rgap=REFERENCE_GAP_TOLERANCE)
    output:
       "results/haplogroup_clustering/para_cluster_influence/clusters_per_chrom.png",
       "results/haplogroup_clustering/para_cluster_influence/num_of_cluster_correlation_per_chrom.png",
       "results/haplogroup_clustering/para_cluster_influence/ac_per_chrom.png",
       "results/haplogroup_clustering/para_cluster_influence/ac_correlation_per_chrom.png",
       "results/haplogroup_clustering/para_cluster_influence/svlen_per_chrom.png",
       "results/haplogroup_clustering/para_cluster_influence/svlen_correlation_per_chrom.png",
       "results/haplogroup_clustering/para_cluster_influence/Interaction_per_chrom.png"
    params:
        root_dir = f"sv_extraction/{SV_MODE}",
        output_dir = "results/haplogroup_clustering/para_cluster_influence",
        pattern = r"mapq(?P<min_mapq_svim>\d+)_minsvlen(?P<minsvlen>\d+)_maxpart(?P<maxpart>\d+)_partdist(?P<partdist>\d+)_maxedit(?P<maxedit>\d+)_maxsv(?P<max_sv_size>\d+)_qgap(?P<query_gap_tolerance>\d+)_rgap(?P<reference_gap_tolerance>\d+)",
        params_to_analyze = "min_mapq_svim minsvlen max_sv_size query_gap_tolerance reference_gap_tolerance",
        script = "workflow/scripts/understand_data/analyze_clustering_params.py"
    resources:
        mem_mb = 4000
    shell:
        """
        python3 {params.script} \
            --root_dir {params.root_dir} \
            --output_dir {params.output_dir} \
            --pattern '{params.pattern}' \
            --params_to_analyze {params.params_to_analyze}
        """
