#!/bin/bash
#SBATCH --job-name=canon_decomp
#SBATCH --output=%x_%A_%a.log
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
ROOT=${CENSV_ROOT}
SUF=mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200
VCF=$ROOT/variant_calling/svim_asm/cohort_svim_asm_${SUF}.vcf.gz
OVL=$ROOT/resources/overlap_samples.txt
OUT=$ROOT/analyses/cdr_qtl/canonical_geno/decompose
PIXI=$HOME/.pixi/bin/pixi
CLUST=$ROOT/workflow/scripts/haplogroup_clustering.py

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)
if [ -n "${1:-}" ]; then CHROM=$1; else CHROM=${CHROMS[$((SLURM_ARRAY_TASK_ID-1))]}; fi

mkdir -p "$OUT"
echo "[$(date)] canonical decompose $CHROM  vcf=$VCF  overlap=$OVL -> $OUT"
cd "$ROOT"
$PIXI run python "$CLUST" \
  --vcf "$VCF" \
  --chrom "$CHROM" \
  --method cosine_nj \
  --decompose \
  --overlap --overlap-file "$OVL" \
  --output "$OUT"
echo "[$(date)] DONE $CHROM"
ls -la "$OUT/$CHROM/" | grep -E "sv_genotypes|sv_variants|sv_samples" || echo "WARN: expected outputs missing for $CHROM"
