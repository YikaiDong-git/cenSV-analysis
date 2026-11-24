#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
#SBATCH --job-name=cdr_rescan
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
set -euo pipefail

PROJ=${CENSV_ROOT}
RD=$PROJ/analyses/cdr_phenotypes
WU=$PROJ/analyses/cdr_qtl_full_panel
CANON=$PROJ/analyses/cdr_qtl

NPERM=${NPERM:-100000}
TASKFILE=${TASKFILE:-$RD/output/tasks/rescan_tasks.txt}

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKFILE")
CHROM=$(echo "$LINE" | awk '{print $1}')
PHENO=$(echo "$LINE" | awk '{print $2}')
FRAME=$(echo "$LINE" | awk '{print $3}')
[ -z "$CHROM" ] && { echo "no task for index $SLURM_ARRAY_TASK_ID"; exit 1; }
[ -z "$FRAME" ] && { echo "task line has no frame: '$LINE'"; exit 1; }

COMPAT=$RD/output/compat/compat_$FRAME
if [ "$PHENO" = "cdr_pos_offset_deepest" ]; then
  PHENOFILE=$COMPAT/per_hap_chrom_phenotypes_deepest.tsv
else
  PHENOFILE=$COMPAT/per_hap_chrom_phenotypes.tsv
fi
[ -s "$PHENOFILE" ] || { echo "missing phenotype file: $PHENOFILE"; exit 1; }

OUTBASE=${OUTBASE:-$RD/output/scan/scan_1e5}
OUTDIR=$OUTBASE/$FRAME
mkdir -p "$OUTDIR"
echo "[$(date)] rescan $CHROM / $PHENO / frame=$FRAME  nperm=$NPERM"
echo "[pheno-file] $PHENOFILE"

cd "$PROJ"
~/.pixi/bin/pixi run python "$CANON/code/scan/b6_cdr_pos_qtl_perchrom.py" \
  --chrom "$CHROM" --phenotype "$PHENO" \
  --sv-source "$WU/geno/perchrom" \
  --pheno-file "$PHENOFILE" \
  --n-perm "$NPERM" \
  --extra-cols sv_locus,geno_group,status,cv_r2,strategy,arm \
  --out-dir "$OUTDIR"

echo "[$(date)] DONE $CHROM / $PHENO / $FRAME"
