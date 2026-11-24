#!/bin/bash
#SBATCH --job-name=chr10_recover
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=1:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_qtl/logs/chr10_recover_%j.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
set -euo pipefail

P=${CENSV_ROOT}
VC=$P/variant_calling/svim_asm
PARAMS=mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200
DEC=$P/sv_extraction/decompose/$PARAMS
PIXI=~/.pixi/bin/pixi
WORK=${CENSV_SCRATCH}/chr10_recover_${SLURM_JOB_ID:-manual}
mkdir -p "$WORK/decompose_out"
cd "$P"

echo "[$(date +%T)] === STEP 1: rebuild cohort VCF (exact original: plain bcftools merge, no flags) ==="
ls "$VC"/*_hap*_${PARAMS}.vcf.gz | sort > "$WORK/vcf_list.txt"
echo "  n input per-hap VCFs: $(wc -l < "$WORK/vcf_list.txt")"
COHORT="$WORK/cohort_svim_asm_${PARAMS}.vcf.gz"

set +e
$PIXI run bcftools merge --file-list "$WORK/vcf_list.txt" -O z -o "$COHORT" 2> "$WORK/merge.err"
rc=$?
set -e
if [ $rc -ne 0 ]; then
  echo "  PLAIN MERGE FAILED (rc=$rc). stderr tail:"
  tail -5 "$WORK/merge.err"
  echo "  This means plain merge is not byte-reproducible in this bcftools; STOP and report."
  exit 3
fi
$PIXI run tabix -f -p vcf "$COHORT"
echo "  cohort built: $(du -h "$COHORT" | cut -f1)"
echo "  chr10 records in cohort: $($PIXI run bcftools view -H "$COHORT" chr10 2>/dev/null | wc -l)"

echo "[$(date +%T)] === STEP 2: re-decompose chr10 (committed params, jaccard_nj default) ==="
$PIXI run python3 workflow/scripts/haplogroup_clustering.py \
  --vcf "$COHORT" --output "$WORK/decompose_out" --chrom chr10 \
  --min-samples 2 --k-range 2,30 --decompose

NEW="$WORK/decompose_out/chr10/chr10_sv_genotypes.npz"
NEWVAR="$WORK/decompose_out/chr10/chr10_sv_variants.tsv"
OLD="$DEC/chr10/chr10_sv_genotypes.npz"

echo "[$(date +%T)] === STEP 3: GATE — regenerated npz must array-match existing ==="
$PIXI run python3 - "$NEW" "$OLD" "$NEWVAR" <<'PY'
import sys, numpy as np, pandas as pd
new_p, old_p, newvar_p = sys.argv[1:4]
new = np.load(new_p); old = np.load(old_p)
nk = list(new.keys())[0]; ok = list(old.keys())[0]
A = new[nk]; B = old[ok]
print("regenerated:", A.shape, A.dtype, "sum", int(A.sum()))
print("existing:   ", B.shape, B.dtype, "sum", int(B.sum()))
exact = A.shape == B.shape and np.array_equal(A, B)
print("EXACT_MATCH:", exact)
if not exact and A.shape == B.shape:
    print("n_differing_cells:", int((A != B).sum()))
    As = A[np.lexsort(A.T[::-1])]; Bs = B[np.lexsort(B.T[::-1])]
    print("ROW_PERMUTATION_MATCH:", np.array_equal(As, Bs))
elif not exact:
    print("SHAPE_MISMATCH: cannot be the same decompose; merge/method differ")
nv = pd.read_csv(newvar_p, sep="\t")
print("regenerated sv_variants rows:", len(nv), "cols:", list(nv.columns)[:6])
PY
echo "[$(date +%T)] === DONE. If EXACT_MATCH True -> copy chr10_sv_variants.tsv into canonical dir (manual, gated). ==="
echo "WORK=$WORK"
