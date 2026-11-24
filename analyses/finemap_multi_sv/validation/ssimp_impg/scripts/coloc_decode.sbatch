#!/bin/bash
#SBATCH --job-name=coloc_decode
#SBATCH --mem=24G
#SBATCH -c 4
#SBATCH --nodes=1
#SBATCH --exclude=node01,node02,node03,node04,node09,node10,node11,node12,node13,node14
#SBATCH --time=00:40:00
#SBATCH --output=${CENSV_ROOT}/analyses/finemap_multi_sv/validation/ssimp_impg/work/coloc_pros1/coloc_decode.%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
P=${CENSV_WORKSPACE}
BASE=$P/Project1_Centromere/analyses/finemap_multi_sv
V=$BASE/validation/ssimp_impg
WD=$V/work/coloc_pros1
SCR=${CENSV_SCRATCH}/pqtl_pros1
DECODE=$SCR/Proteomics_SMP_PC0_2700_56_PROS1_Protein_S_10032022.txt.gz
PIXI=~/.pixi/bin/pixi
RSTATS_MANIFEST=$P/_shared/envs/workspace/v1/pixi.toml
BT=$P/_shared/envs/workspace/v1/.pixi/envs/biotools/bin
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
mkdir -p "$WD"
echo "[$(date)] slice deCODE Protein-S cis (chr3:93.79-94.92 Mb)"
zcat "$DECODE" | awk -F'\t' 'NR==1 || ($1=="chr3" && $2>=93790000 && $2<=94920000)' > "$WD/decode_proteinS_PROS1.b38.slice.tsv"
wc -l "$WD/decode_proteinS_PROS1.b38.slice.tsv"
echo "[$(date)] prep (harmonize vs FinnGen+UKBB + build FinnGen LD)"
"$BT/python" "$V/scripts/coloc_decode_prep.py"
echo "[$(date)] coloc R (coloc.abf + coloc.signals)"
"$PIXI" run --manifest-path "$RSTATS_MANIFEST" -e rstats Rscript "$V/scripts/coloc_decode.R"
echo "[$(date)] DONE_DECODE_COLOC_JOB"
