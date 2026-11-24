#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${HOME:?set HOME (see README)}"
set -uo pipefail
SUF=mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200
P1=${CENSV_ROOT}
BC=${HOME}/.pixi/envs/bcftools/bin/bcftools
PY=$P1/.pixi/envs/default/bin/python
VC=$P1/variant_calling/svim_asm
CANON=$P1/sv_extraction/decompose/$SUF/chr10
WORK=${CENSV_SCRATCH}/chr10_regen2
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$P1"

echo "[1] vcf list"; ls "$VC"/*_hap?_"$SUF".vcf.gz > "$WORK/vcflist.txt"; echo "    n=$(wc -l < "$WORK/vcflist.txt")"
echo "[2] merge -m none -r chr10"
"$BC" merge -m none -r chr10 -l "$WORK/vcflist.txt" -Oz -o "$WORK/cohort.vcf.gz" 2> "$WORK/merge.err"
"$BC" index -t "$WORK/cohort.vcf.gz"
echo "    records: $("$BC" view -H "$WORK/cohort.vcf.gz" 2>/dev/null | wc -l)  samples: $("$BC" query -l "$WORK/cohort.vcf.gz" 2>/dev/null | wc -l)"
echo "[3] decompose chr10 cosine_nj"
"$PY" workflow/scripts/haplogroup_clustering.py --vcf "$WORK/cohort.vcf.gz" --output "$WORK/decomp_out" \
    --chrom chr10 --method cosine_nj --decompose --min-samples 2 --k-range 2,30 > "$WORK/decomp.log" 2>&1
NEWDIR="$WORK/decomp_out/chr10"; echo "    outputs: $(ls "$NEWDIR" 2>/dev/null | tr '\n' ' ')"
EXIST="$CANON/chr10_sv_genotypes.npz"; NEW="$NEWDIR/chr10_sv_genotypes.npz"
echo "[4] verify"
if [ -f "$NEW" ] && cmp -s "$EXIST" "$NEW"; then
  echo "BYTEMATCH_PASS"; cp "$NEWDIR/chr10_sv_variants.tsv" "$CANON/chr10_sv_variants.tsv"; echo "RECOVERED (byte) -> $CANON"
elif [ -f "$NEW" ]; then
  echo "BYTEMATCH_FAIL -> try genotype-vector alignment"
  "$PY" - "$EXIST" "$NEW" "$NEWDIR/chr10_sv_variants.tsv" "$WORK/chr10_sv_variants_ALIGNED.tsv" <<'PYEOF'
import sys, numpy as np, pandas as pd
E=np.load(sys.argv[1]); R=np.load(sys.argv[2])
Eg=E[list(E.files)[0]]; Rg=R[list(R.files)[0]]
rv=pd.read_csv(sys.argv[3], sep="\t")
print(f"  shapes E={Eg.shape} R={Rg.shape} regen_var_rows={len(rv)}")
if Eg.shape!=Rg.shape:
    print("  shape mismatch -> not equivalent; STOP"); sys.exit(0)
if Eg.shape[0]<Eg.shape[1]: Eg=Eg  # rows=SVs assumed
from collections import defaultdict
idx=defaultdict(list)
for j in range(Rg.shape[0]): idx[Rg[j].tobytes()].append(j)
aligned=[None]*Eg.shape[0]; used=set(); n_ok=0; n_amb=0; n_miss=0
for i in range(Eg.shape[0]):
    cands=[j for j in idx.get(Eg[i].tobytes(),[]) if j not in used]
    if len(cands)==1: aligned[i]=cands[0]; used.add(cands[0]); n_ok+=1
    elif len(cands)>1: aligned[i]=cands[0]; used.add(cands[0]); n_amb+=1
    else: n_miss+=1
print(f"  genotype-vector alignment: unique={n_ok} ambiguous={n_amb} unmatched={n_miss} / {Eg.shape[0]}")
if n_miss==0 and n_amb < 0.05*Eg.shape[0]:
    out=rv.iloc[[aligned[i] for i in range(Eg.shape[0])]].reset_index(drop=True)
    out.to_csv(sys.argv[4], sep="\t", index=False)
    print(f"  WROTE aligned variant table ({len(out)} rows) -> {sys.argv[4]}")
    print("  ALIGN_PASS")
else:
    print("  ALIGN_FAIL (too many unmatched/ambiguous -> not faithfully recoverable)")
PYEOF
else
  echo "BYTEMATCH_FAIL (no regen npz produced); decomp.log tail:"; tail -5 "$WORK/decomp.log"
fi
echo "CHR10_REGEN2_DONE"
