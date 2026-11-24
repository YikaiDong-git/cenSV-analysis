#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail

PROJECT=${CENSV_ROOT}
Q100_DIR=${PROJECT}/resources/hg002_q100
WORK=${PROJECT}/analyses/imputability/output/hg002_q100_truth
LOGS=${PROJECT}/analyses/imputability/logs
mkdir -p "${WORK}" "${LOGS}"

SBATCH_SCRIPT=${WORK}/run.sh
cat > "${SBATCH_SCRIPT}" << 'BATCH'
#SBATCH --job-name=q100_align_call
#SBATCH --output=DISCOVERY_LOGS/q100_align_call_%j.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=05:00:00
set -euo pipefail
PROJECT=${CENSV_ROOT}
Q100_DIR=${PROJECT}/resources/hg002_q100
WORK=${PROJECT}/analyses/imputability/output/hg002_q100_truth
cd ${PROJECT}

export PATH=${PROJECT}/.pixi/envs/default/bin:$PATH
export PYTHONPATH=${PROJECT}/software/svim-asm/src${PYTHONPATH:+:$PYTHONPATH}

SVIM=${PROJECT}/software/svim-asm/src/svim_asm/svim-asm
REF=${PROJECT}/resources/chm13v2.0.fa.gz
Q100=${Q100_DIR}/hg002v1.1.fasta.gz

echo "[$(date -Iseconds)] PROVENANCE"
echo "  repo commit: $(git rev-parse HEAD 2>/dev/null || echo NA)"
echo "  minimap2:    $(minimap2 --version 2>&1 | head -1)"
echo "  samtools:    $(samtools --version 2>&1 | head -1)"
echo "  bcftools:    $(bcftools --version 2>&1 | head -1)"
echo "  python3:     $(python3 --version 2>&1)"
echo "  svim-asm:    $(python3 ${SVIM} --version 2>&1 | head -1 || true)"
echo "  Q100 md5:    $(md5sum ${Q100} | cut -d' ' -f1)"
echo "  Q100 fai:    ${Q100}.fai ($(wc -l < ${Q100}.fai) contigs)"
echo

for HAP in PATERNAL MATERNAL; do
    OUT=${WORK}/q100_${HAP}.fa
    if [[ ! -s "${OUT}.gz" ]]; then
        echo "[$(date -Iseconds)] Splitting ${HAP} hap..."
        REGIONS=$(awk -v hap=_${HAP} '$1 ~ hap {print $1}' ${Q100}.fai)
        samtools faidx ${Q100} ${REGIONS} > ${OUT}
        bgzip -@ 8 -f ${OUT}
        samtools faidx ${OUT}.gz
        echo "[$(date -Iseconds)] ${HAP} hap written: $(ls -la ${OUT}.gz)"
    else
        echo "[$(date -Iseconds)] ${OUT}.gz already exists, skip split."
    fi
done

for HAP in PATERNAL MATERNAL; do
    BAM=${WORK}/q100_${HAP}_to_chm13.bam
    if [[ ! -s "${BAM}.csi" ]]; then
        echo "[$(date -Iseconds)] minimap2 align ${HAP}..."
        minimap2 -a -x asm5 --cs -r2k -t 16 ${REF} ${WORK}/q100_${HAP}.fa.gz \
            | samtools sort -@ 4 -m 4G -o ${BAM} --write-index
        echo "[$(date -Iseconds)] ${HAP} BAM: $(ls -la ${BAM})"
    else
        echo "[$(date -Iseconds)] ${BAM} already aligned, skip."
    fi
done

for HAP in PATERNAL MATERNAL; do
    OUTDIR=${WORK}/svim_asm_${HAP}
    OUTVCF=${WORK}/q100_${HAP}_mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200.vcf.gz
    if [[ ! -s "${OUTVCF}" ]]; then
        echo "[$(date -Iseconds)] svim-asm call ${HAP}..."
        rm -rf ${OUTDIR}
        mkdir -p ${OUTDIR}
        python3 ${SVIM} haploid ${OUTDIR} ${WORK}/q100_${HAP}_to_chm13.bam ${REF} \
            --sample q100_${HAP} \
            --min_mapq 20 \
            --min_sv_size 10 \
            --max_partition_size 10 \
            --partition_max_distance 1000 \
            --max_edit_distance 200 \
            --max_sv_size 500000 \
            --query_gap_tolerance 200 \
            --reference_gap_tolerance 200
        bcftools sort ${OUTDIR}/variants.vcf -O z -o ${OUTVCF}
        tabix -p vcf ${OUTVCF}
        echo "[$(date -Iseconds)] ${HAP} VCF written: $(ls -la ${OUTVCF})"
    else
        echo "[$(date -Iseconds)] ${OUTVCF} already exists, skip."
    fi
done

ACTIVE=${PROJECT}/resources/chm13v2.0_censat_v2.1.active.bed
echo
echo "=== Q100 SV counts (genome-wide + active HOR) ==="
for HAP in PATERNAL MATERNAL; do
    VCF=${WORK}/q100_${HAP}_mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200.vcf.gz
    GW=$(bcftools view -H ${VCF} | wc -l)
    HOR=$(bcftools view -H -R ${ACTIVE} ${VCF} | wc -l)
    echo "  ${HAP}: genome-wide=${GW} active-HOR=${HOR}"
done

echo
echo "[$(date -Iseconds)] DONE"
BATCH

sed -i "s|DISCOVERY_LOGS|${LOGS}|g" "${SBATCH_SCRIPT}"

echo "Submitting sbatch..."
JOBID=$(sbatch --parsable "${SBATCH_SCRIPT}")
echo "[OK] q100_align_call sbatch JOBID=${JOBID}"
echo "  log: ${LOGS}/q100_align_call_${JOBID}.log"
