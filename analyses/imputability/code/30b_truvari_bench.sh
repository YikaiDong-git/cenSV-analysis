#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail

PROJECT=${CENSV_ROOT}
WORK=${PROJECT}/analyses/imputability/output/hg002_q100_truth
LOGS=${PROJECT}/analyses/imputability/logs

SUFFIX=mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200
PROD_H1=${PROJECT}/variant_calling/svim_asm/HG002_hap1_${SUFFIX}.vcf.gz
PROD_H2=${PROJECT}/variant_calling/svim_asm/HG002_hap2_${SUFFIX}.vcf.gz
TRUTH_P=${WORK}/q100_PATERNAL_${SUFFIX}.vcf.gz
TRUTH_M=${WORK}/q100_MATERNAL_${SUFFIX}.vcf.gz

if [[ ! -s "${TRUTH_P}" || ! -s "${TRUTH_M}" ]]; then
    echo "[ERR] Q100 truth VCFs missing; Phase A3 must finish first."
    exit 1
fi

mkdir -p "${LOGS}"

SBATCH=${WORK}/run_truvari.sh
cat > "${SBATCH}" << BATCH
#SBATCH --job-name=truvari_bench
#SBATCH --output=${LOGS}/truvari_bench_%j.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
cd ${PROJECT}
export PATH=${PROJECT}/.pixi/envs/default/bin:\$PATH

echo "[\$(date -Iseconds)] renaming + merging haps..."
echo HG002 > ${WORK}/_sample_name.txt
for f in ${PROD_H1} ${PROD_H2} ${TRUTH_P} ${TRUTH_M}; do
    bn=\$(basename \$f .vcf.gz)
    bcftools reheader -s ${WORK}/_sample_name.txt \$f -o ${WORK}/_\${bn}_renamed.vcf.gz
    tabix -p vcf -f ${WORK}/_\${bn}_renamed.vcf.gz
done
bcftools concat -a -D ${WORK}/_HG002_hap1_*.vcf.gz ${WORK}/_HG002_hap2_*.vcf.gz \\
    -O z -o ${WORK}/hg002_prod_union.vcf.gz
tabix -p vcf -f ${WORK}/hg002_prod_union.vcf.gz
bcftools concat -a -D ${WORK}/_q100_PATERNAL_*.vcf.gz ${WORK}/_q100_MATERNAL_*.vcf.gz \\
    -O z -o ${WORK}/hg002_q100_truth_union.vcf.gz
tabix -p vcf -f ${WORK}/hg002_q100_truth_union.vcf.gz

TRUVARI_OUT=${WORK}/truvari_bench_out
rm -rf \${TRUVARI_OUT}
echo "[\$(date -Iseconds)] truvari bench..."
truvari bench \\
    -b ${WORK}/hg002_q100_truth_union.vcf.gz \\
    -c ${WORK}/hg002_prod_union.vcf.gz \\
    -o \${TRUVARI_OUT} \\
    --refdist 500 \\
    --pctseq 0.7 \\
    --pctsize 0.7 \\
    --pick single \\
    --reference ${PROJECT}/resources/chm13v2.0.fa.gz \\
    --includebed ${PROJECT}/resources/chm13v2.0_censat_v2.1.active.bed

echo "[\$(date -Iseconds)] truvari bench summary:"
cat \${TRUVARI_OUT}/summary.json 2>&1 | head -50
BATCH

JOBID=$(sbatch --parsable "${SBATCH}")
echo "[OK] truvari_bench JOBID=${JOBID}"
echo "  log: ${LOGS}/truvari_bench_${JOBID}.log"
