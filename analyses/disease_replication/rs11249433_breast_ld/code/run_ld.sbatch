#!/bin/bash
#SBATCH -J ld_rs11249433
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH --nodes=1
#SBATCH -o ${CENSV_ROOT}/analyses/disease_replication/rs11249433_breast_ld/logs/%x.%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
export OMP_NUM_THREADS=1
BASE=${CENSV_WORKSPACE}
TRI=$BASE/Project1_Centromere/analyses/disease_replication
WD=$TRI/rs11249433_breast_ld
HITS=$HOME/AUO_handoff${CENSV_ROOT}/analyses/disease_replication/aou_replication/results_vds/step30_phewide_ext2/phewide_disease_hits_enriched_flagged.tsv
WTS=$TRI/phewide_all_handoff/full_11129_ext2/sv_weights_long_ext2_chr1.tsv
EUR=$BASE/_shared/datasets/1KGP_PLINK_b38_EUR_20260530/chr1_eur
PLINK2=$HOME/.pixi/bin/plink2
PY=$BASE/Project1_Centromere/.pixi/envs/default/bin/python3
SCR=${CENSV_SCRATCH}/p1_rs11249433_ld_${SLURM_JOB_ID}
mkdir -p "$SCR" "$WD/output"

$PLINK2 --bfile "$EUR" --chr 1 --from-bp 120000000 --to-bp 147000000 \
  --set-all-var-ids @:# --rm-dup exclude-all --mac 1 --threads 1 \
  --make-pgen --out "$SCR/win" > /dev/null

$PY "$WD/code/make_score.py" "$HITS" "$WTS" "$SCR/win.pvar" "$SCR/score.tsv" "$WD/output/sv_coverage.tsv"
NCOL=$(head -1 "$SCR/score.tsv" | awk -F'\t' '{print NF}')

$PLINK2 --pfile "$SCR/win" --threads 1 \
  --score "$SCR/score.tsv" 1 2 header-read variance-standardize cols=+scoresums \
  --score-col-nums 3-"$NCOL" --out "$SCR/sv" > /dev/null
printf '1:121538815\n1:143272445\n' > "$SCR/snps.txt"
$PLINK2 --pfile "$SCR/win" --threads 1 --extract "$SCR/snps.txt" --export A --out "$SCR/snps" > /dev/null

$PY "$WD/code/ld_table.py" "$SCR/score.tsv" "$SCR/sv.sscore" "$SCR/snps.raw" "$WD/output/sv_coverage.tsv" "$WD/output/ld_breast_svs_vs_rs11249433.tsv"
cp "$SCR/sv.log" "$WD/logs/plink_score.${SLURM_JOB_ID}.log"
rm -rf "$SCR"
echo DONE
