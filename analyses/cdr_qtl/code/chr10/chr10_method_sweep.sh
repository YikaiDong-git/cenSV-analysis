#!/bin/bash
#SBATCH --job-name=chr10_methsweep
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=1:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_qtl/logs/chr10_methsweep_%j.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
set -euo pipefail

P=${CENSV_ROOT}
VC=$P/variant_calling/svim_asm
PARAMS=mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200
DEC=$P/sv_extraction/decompose/$PARAMS
PIXI=~/.pixi/bin/pixi
WORK=${CENSV_SCRATCH}/chr10_recover_3447741   # reuse cohort built by JID 3447741
COHORT="$WORK/cohort_svim_asm_${PARAMS}.vcf.gz"
cd "$P"

if [ ! -f "$COHORT" ]; then
  echo "[$(date +%T)] cohort missing in scratch -> rebuild (plain merge)"
  mkdir -p "$WORK"
  ls "$VC"/*_hap*_${PARAMS}.vcf.gz | sort > "$WORK/vcf_list.txt"
  $PIXI run bcftools merge --file-list "$WORK/vcf_list.txt" -O z -o "$COHORT"
  $PIXI run tabix -f -p vcf "$COHORT"
fi
echo "[$(date +%T)] cohort: $(du -h "$COHORT" | cut -f1)  chr10_input_records: $($PIXI run bcftools view -H "$COHORT" chr10 2>/dev/null | wc -l)"

OLD_NPZ="$DEC/chr10/chr10_sv_genotypes.npz"
OLD_DECOMP="$DEC/chr10/chr10_decomposition.tsv"
echo "[$(date +%T)] TARGET: npz $($PIXI run python3 -c "import numpy as np; d=np.load('$OLD_NPZ'); a=d[list(d.keys())[0]]; print(a.shape, int(a.sum()))")  decomp_rows=$(wc -l < "$OLD_DECOMP")"

for M in cosine_nj correlation_nj jaccard_ward cosine_ward; do
  OUT="$WORK/dec_$M"
  echo "===================================================================="
  echo "[$(date +%T)] METHOD=$M"
  set +e
  $PIXI run python3 workflow/scripts/haplogroup_clustering.py \
    --vcf "$COHORT" --output "$OUT" --chrom chr10 \
    --min-samples 2 --k-range 2,30 --decompose --method "$M" \
    > "$WORK/log_$M.txt" 2>&1
  rc=$?
  set -e
  if [ $rc -ne 0 ]; then echo "  decompose rc=$rc; tail:"; tail -4 "$WORK/log_$M.txt"; continue; fi
  NEW_NPZ="$OUT/chr10/chr10_sv_genotypes.npz"
  NEW_DECOMP="$OUT/chr10/chr10_decomposition.tsv"
  if [ ! -f "$NEW_NPZ" ]; then echo "  no npz (method makes no nj_tree -> no decompose)"; continue; fi
  $PIXI run python3 - "$NEW_NPZ" "$OLD_NPZ" "$NEW_DECOMP" "$OLD_DECOMP" <<'PY'
import sys, numpy as np
nn,on,nd,od = sys.argv[1:5]
A=np.load(nn); B=np.load(on); A=A[list(A.keys())[0]]; B=B[list(B.keys())[0]]
import subprocess
nrows=sum(1 for _ in open(nd)); orows=sum(1 for _ in open(od))
print(f"  regen npz {A.shape} sum {int(A.sum())} | decomp_rows {nrows}")
print(f"  exist npz {B.shape} sum {int(B.sum())} | decomp_rows {orows}")
exact = A.shape==B.shape and np.array_equal(A,B)
print(f"  EXACT_NPZ_MATCH={exact}")
if A.shape==B.shape and not exact:
    print(f"  n_diff_cells={int((A!=B).sum())}")
PY
done
echo "[$(date +%T)] === SWEEP DONE ==="
echo "WORK=$WORK"
