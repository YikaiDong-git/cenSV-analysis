#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
PROJ=${CENSV_ROOT}
cd $PROJ
SUFFIX="mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
STRATS=(fksize100k fksize200k fksize500k edge500 edge1000 edge2000 midpoint500 midpoint1000)
CHUNKED_CHROMS=(17 1 21 19 18)
declare -A CHUNKS=( [17]=30 [1]=44 [21]=26 [19]=46 [18]=30 )
SINGLE_CHROMS=(15 4 6 9 16 12 3 8 7 20 2 10 13 14 5 22 11 X)

mkdir -p logs/en_chunk logs/en_single

strat_to_snpdir() {
  case "$1" in
    fksize100k) echo "snp_extraction/fksize100000";;
    fksize200k) echo "snp_extraction/fksize200000";;
    fksize500k) echo "snp_extraction/fksize500000";;
    edge500|edge1000|edge2000|midpoint500|midpoint1000) echo "snp_extraction/$1";;
  esac
}

N_ARR=0; N_MRG=0; N_SINGLE=0

for chr in "${CHUNKED_CHROMS[@]}"; do
  N=${CHUNKS[$chr]}
  LAST=$((N-1))
  SV_DIR="sv_extraction/decompose/$SUFFIX/chr${chr}"
  for strat in "${STRATS[@]}"; do
    SNP_DIR=$(strat_to_snpdir $strat)
    EN_DIR="results/elastic_net/decompose/cosine_nj/${strat}_${SUFFIX}/chr${chr}"
    mkdir -p "$EN_DIR/chunks"

    ARR_JID=$(sbatch --parsable --array=0-${LAST}%50 --mem=16G -c 12 --time=24:00:00 \
      --job-name=en_chr${chr}_${strat} \
      --output=logs/en_chunk/chr${chr}_${strat}_%A_%a.log --open-mode=append \
      --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/elastic_net_train_chunk.py \
        --snp-dir $SNP_DIR \
        --sv-dir  $SV_DIR \
        --output  $EN_DIR/chunks/chunk_\$SLURM_ARRAY_TASK_ID \
        --chromosome chr${chr} \
        --sv-start  \$((SLURM_ARRAY_TASK_ID * 50)) \
        --sv-end    \$(((SLURM_ARRAY_TASK_ID + 1) * 50)) \
        --n-folds 4 --n-jobs 12 \
        --min-snp-freq 0.05 --min-sv-freq 0.05 --max-sv-freq 0.99 \
        --sv-filter all")
    echo "array chr${chr}_${strat} jobid=$ARR_JID tasks=$N"
    N_ARR=$((N_ARR+1))

    MRG_JID=$(sbatch --parsable --dependency=afterok:$ARR_JID --mem=4G -c 2 --time=1:00:00 \
      --job-name=merge_chr${chr}_${strat} \
      --output=logs/en_chunk/merge_chr${chr}_${strat}_%A.log \
      --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/merge_en_chunks.py \
        --chunks-dir $EN_DIR/chunks \
        --output     $EN_DIR")
    echo "merge chr${chr}_${strat} jobid=$MRG_JID (dep=$ARR_JID)"
    N_MRG=$((N_MRG+1))
  done
done

for chr in "${SINGLE_CHROMS[@]}"; do
  SV_DIR="sv_extraction/decompose/$SUFFIX/chr${chr}"
  for strat in "${STRATS[@]}"; do
    SNP_DIR=$(strat_to_snpdir $strat)
    EN_DIR="results/elastic_net/decompose/cosine_nj/${strat}_${SUFFIX}/chr${chr}"
    mkdir -p "$EN_DIR"

    S_JID=$(sbatch --parsable --mem=16G -c 12 --time=24:00:00 \
      --job-name=en_chr${chr}_${strat} \
      --output=logs/en_single/chr${chr}_${strat}_%A.log \
      --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/elastic_net_train.py \
        --snp-dir $SNP_DIR \
        --sv-dir  $SV_DIR \
        --output  $EN_DIR \
        --chromosome chr${chr} \
        --n-folds 4 --n-jobs 12 \
        --min-snp-freq 0.05 --min-sv-freq 0.05 --max-sv-freq 0.99 \
        --sv-filter all")
    echo "single chr${chr}_${strat} jobid=$S_JID"
    N_SINGLE=$((N_SINGLE+1))
  done
done

echo "dispatched: arrays=$N_ARR merges=$N_MRG singles=$N_SINGLE"
