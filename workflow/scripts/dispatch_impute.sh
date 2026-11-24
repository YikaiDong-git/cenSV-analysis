#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
PROJ=${CENSV_ROOT}
cd $PROJ
SUFFIX="mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
STRATS=(fksize100k fksize200k fksize500k edge500 edge1000 edge2000 midpoint500 midpoint1000)
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
R2=0.1
SAMPLES=resources/imputation_target_samples.txt

mkdir -p logs/impute
JIDS=""

for strat in "${STRATS[@]}"; do
  case "$strat" in
    fksize100k) SNP_SUFFIX=fksize100000;;
    fksize200k) SNP_SUFFIX=fksize200000;;
    fksize500k) SNP_SUFFIX=fksize500000;;
    *) SNP_SUFFIX=$strat;;
  esac
  SH=workflow/scripts/impute_${strat}.sbatch
  cat > $SH <<EOF
#SBATCH --array=0-22%23
#SBATCH --mem=16G
#SBATCH -c 1
#SBATCH --time=2:00:00
#SBATCH --job-name=impute_${strat}
#SBATCH --output=logs/impute/${strat}_%A_%a.log
set -euo pipefail
cd $PROJ
CHROMS=(${CHROMS[@]})
CH=\${CHROMS[\$SLURM_ARRAY_TASK_ID]}
MODEL_DIR=results/elastic_net/decompose/cosine_nj/${strat}_${SUFFIX}/chr\${CH}
OUT_DIR=results/imputation/decompose/cosine_nj/${strat}_${SUFFIX}/chr\${CH}
mkdir -p \$OUT_DIR
~/.pixi/bin/pixi run python3 workflow/scripts/impute_1kg.py \\
  --model-dir \$MODEL_DIR \\
  --vcf resources/1000G/1kGP_high_coverage_Illumina.chr\${CH}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \\
  --samples $SAMPLES \\
  --output \$OUT_DIR \\
  --chromosome chr\${CH} \\
  --r2-threshold $R2
EOF
  JID=$(sbatch --parsable $SH)
  echo "impute $strat jobid=$JID"
  JIDS="${JIDS}:${JID}"
done

JIDS_CLEAN=${JIDS#:}
echo "DEPS=afterok:$JIDS_CLEAN"
echo "dispatched 8 arrays × 23 tasks = 184 imputation jobs"
echo "$JIDS_CLEAN" > logs/impute/.dep_jids.txt
