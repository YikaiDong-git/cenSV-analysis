CIS_WINDOWS = config.get("cis_windows", [10000000, 0])   # 10 Mb + whole-chrom
EQTL_PEER = config.get("eqtl_peer", 60)
SQTL_PEER = config.get("sqtl_peer", 15)
GENO_PC = config.get("geno_pc", 5)
INCLUDE_SEX = config.get("include_sex", True)

COV_STRING = f"sex{INCLUDE_SEX}_gpc{GENO_PC}_peer"

MAGE_EXPRESSION_BED = config.get("mage_expression_bed", "resources/MAGE/expression_bed.bed.gz")
MAGE_SPLICING_BED = config.get("mage_splicing_bed", "resources/MAGE/splicing_bed.bed.gz")
MAGE_PHENOTYPE_GROUPS = config.get("mage_phenotype_groups", "resources/MAGE/splicing_phenotype_groups.txt")
MAGE_EQTL_COV = config.get("mage_eqtl_covariates", "resources/MAGE/eqtl_covariates.tab.gz")
MAGE_SQTL_COV = config.get("mage_sqtl_covariates", "resources/MAGE/sqtl_covariates.tab.gz")

EN_STRATEGY = config.get("en_strategy", "best8")

IMP_BASE = f"results/imputation/decompose/cosine_nj/{EN_STRATEGY}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}"
EQTL_BASE = f"results/eqtl_overlap/decompose/cosine_nj/{EN_STRATEGY}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}"
SQTL_BASE = f"results/sqtl_overlap/decompose/cosine_nj/{EN_STRATEGY}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}"

rule cis_eqtl:
    input:
        dosage = IMP_BASE + "/chr{chr}/chr{chr}_imputed_dosages.tsv",
        centromere_bed = "resources/grch38_centromeres.bed",
        expression = MAGE_EXPRESSION_BED,
        covariates = MAGE_EQTL_COV,
        script = "workflow/scripts/qtl_overlap.py"
    output:
        EQTL_BASE + f"/chr{{chr}}/{COV_STRING}{EQTL_PEER}/window_{{window}}/cis_qtl_results.tsv"
    params:
        imputation_dir = IMP_BASE + "/chr{chr}",
        outdir = EQTL_BASE + f"/chr{{chr}}/{COV_STRING}{EQTL_PEER}/window_{{window}}"
    wildcard_constraints:
        window = r"\d+"
    benchmark:
        f"benchmarks/cis_eqtl/cis_eqtl_{EN_STRATEGY}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}_chr{{chr}}_w{{window}}.tsv"
    resources:
        mem_mb = 16000,
        runtime = 120,
        cpus = 1
    shell:
        """
        python3 {input.script} \
            --imputation-dir {params.imputation_dir} \
            --centromere-bed {input.centromere_bed} \
            --output {params.outdir} \
            --chromosome chr{wildcards.chr} \
            --qtl-type eqtl \
            --expression {input.expression} \
            --eqtl-covariates {input.covariates} \
            --regmode within-chrom \
            --regmode-para {wildcards.window}
        """

rule cis_sqtl:
    input:
        dosage = IMP_BASE + "/chr{chr}/chr{chr}_imputed_dosages.tsv",
        centromere_bed = "resources/grch38_centromeres.bed",
        splicing = MAGE_SPLICING_BED,
        phenotype_groups = MAGE_PHENOTYPE_GROUPS,
        covariates = MAGE_SQTL_COV,
        script = "workflow/scripts/qtl_overlap.py"
    output:
        SQTL_BASE + f"/chr{{chr}}/{COV_STRING}{SQTL_PEER}/window_{{window}}/cis_qtl_results.tsv"
    params:
        imputation_dir = IMP_BASE + "/chr{chr}",
        outdir = SQTL_BASE + f"/chr{{chr}}/{COV_STRING}{SQTL_PEER}/window_{{window}}"
    wildcard_constraints:
        window = r"\d+"
    benchmark:
        f"benchmarks/cis_sqtl/cis_sqtl_{EN_STRATEGY}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}_chr{{chr}}_w{{window}}.tsv"
    resources:
        mem_mb = 16000,
        runtime = 120,
        cpus = 1
    shell:
        """
        python3 {input.script} \
            --imputation-dir {params.imputation_dir} \
            --centromere-bed {input.centromere_bed} \
            --output {params.outdir} \
            --chromosome chr{wildcards.chr} \
            --qtl-type sqtl \
            --splicing-bed {input.splicing} \
            --phenotype-groups {input.phenotype_groups} \
            --sqtl-covariates {input.covariates} \
            --regmode within-chrom \
            --regmode-para {wildcards.window}
        """

rule trans_eqtl:
    input:
        dosage = IMP_BASE + "/chr{chr}/chr{chr}_imputed_dosages.tsv",
        centromere_bed = "resources/grch38_centromeres.bed",
        expression = MAGE_EXPRESSION_BED,
        covariates = MAGE_EQTL_COV,
        script = "workflow/scripts/qtl_overlap.py"
    output:
        EQTL_BASE + f"/chr{{chr}}/{COV_STRING}{EQTL_PEER}/trans/trans_qtl_results.tsv"
    params:
        imputation_dir = IMP_BASE + "/chr{chr}",
        outdir = EQTL_BASE + f"/chr{{chr}}/{COV_STRING}{EQTL_PEER}/trans",
        other_chroms = lambda wildcards: ",".join([c for c in SV_CHRS if c != wildcards.chr])
    benchmark:
        f"benchmarks/trans_eqtl/trans_eqtl_{EN_STRATEGY}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}_chr{{chr}}.tsv"
    resources:
        mem_mb = 32000,
        runtime = 360,
        cpus = 1
    shell:
        """
        python3 {input.script} \
            --imputation-dir {params.imputation_dir} \
            --centromere-bed {input.centromere_bed} \
            --output {params.outdir} \
            --chromosome chr{wildcards.chr} \
            --qtl-type eqtl \
            --expression {input.expression} \
            --eqtl-covariates {input.covariates} \
            --regmode cross-chrom \
            --regmode-para {params.other_chroms}
        """

rule trans_sqtl:
    input:
        dosage = IMP_BASE + "/chr{chr}/chr{chr}_imputed_dosages.tsv",
        centromere_bed = "resources/grch38_centromeres.bed",
        splicing = MAGE_SPLICING_BED,
        phenotype_groups = MAGE_PHENOTYPE_GROUPS,
        covariates = MAGE_SQTL_COV,
        script = "workflow/scripts/qtl_overlap.py"
    output:
        SQTL_BASE + f"/chr{{chr}}/{COV_STRING}{SQTL_PEER}/trans/trans_qtl_results.tsv"
    params:
        imputation_dir = IMP_BASE + "/chr{chr}",
        outdir = SQTL_BASE + f"/chr{{chr}}/{COV_STRING}{SQTL_PEER}/trans",
        other_chroms = lambda wildcards: ",".join([c for c in SV_CHRS if c != wildcards.chr])
    benchmark:
        f"benchmarks/trans_sqtl/trans_sqtl_{EN_STRATEGY}_mapq{{mapq}}_minsvlen{{minsvlen}}_maxpart{{maxpart}}_partdist{{partdist}}_maxedit{{maxedit}}_maxsv{{maxsv}}_qgap{{qgap}}_rgap{{rgap}}_chr{{chr}}.tsv"
    resources:
        mem_mb = 32000,
        runtime = 360,
        cpus = 1
    shell:
        """
        python3 {input.script} \
            --imputation-dir {params.imputation_dir} \
            --centromere-bed {input.centromere_bed} \
            --output {params.outdir} \
            --chromosome chr{wildcards.chr} \
            --qtl-type sqtl \
            --splicing-bed {input.splicing} \
            --phenotype-groups {input.phenotype_groups} \
            --sqtl-covariates {input.covariates} \
            --regmode cross-chrom \
            --regmode-para {params.other_chroms}
        """
