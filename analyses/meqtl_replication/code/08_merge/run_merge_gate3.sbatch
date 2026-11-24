#!/bin/bash
#SBATCH --job-name=mergegate3
#SBATCH --output=${CENSV_ROOT}/analyses/meqtl_replication/logs/v2_lift/mergegate3_%j.log
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
P=${CENSV_ROOT}
U=$P/analyses/meqtl_replication
echo "HOST $(hostname)  START $(date -u +%Y-%m-%dT%H:%M:%SZ)"

N=$(ls $U/output/v2_lift/per_hap/*.tsv.gz 2>/dev/null | wc -l)
W=$(wc -l < $U/output/a1_lift/all_tags.txt)
echo "haplotype files: $N / $W"

cd "$P" || exit 1
~/.pixi/bin/pixi run python "$U/code/08_merge/merge_and_gate3.py"
rc=$?
echo "EXIT $rc  END $(date -u +%Y-%m-%dT%H:%M:%SZ)"
[ $rc -eq 0 ] && echo PASS || echo "FAIL rc=$rc"
exit $rc
