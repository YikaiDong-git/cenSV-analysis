#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail

PROJECT=${CENSV_ROOT}
CODE=${PROJECT}/analyses/imputability/code
LOGS=${PROJECT}/analyses/imputability/logs
OUT=${PROJECT}/analyses/imputability/output/mapq_per_chrom
mkdir -p "${LOGS}" "${OUT}"

JOBIDS=()
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do
    CHROM="chr${c}"
    JID=$(sbatch --parsable \
        --job-name=mapq_${CHROM} \
        --output=${LOGS}/mapq_${CHROM}_%j.log \
        --cpus-per-task=4 \
        --mem=8G \
        --time=04:00:00 \
        --wrap="cd ${PROJECT} && ~/.pixi/bin/pixi run python3 ${CODE}/32_mapq_feature_extract.py --chrom ${CHROM} --out-dir ${OUT}")
    JOBIDS+=("${JID}")
    echo "submitted ${CHROM}: ${JID}"
done

echo
echo "All 23 jobs submitted. JOBIDS: ${JOBIDS[*]}"
echo "Monitor: ssh hgcc 'squeue -u ydon268 -n mapq_chr1,mapq_chr2,...'"
