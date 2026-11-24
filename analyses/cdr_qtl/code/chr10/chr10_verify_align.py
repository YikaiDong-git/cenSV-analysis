import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
P=CENSV_ROOT + ""
PARAMS="mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
DEC=f"{P}/sv_extraction/decompose/{PARAMS}/chr10"
WORK=CENSV_SCRATCH + "/chr10_recover_3447749_NA"
import glob
cands = glob.glob(CENSV_SCRATCH + "/chr10_recover_*/dec_cosine_nj/chr10")
assert cands, "no cosine_nj output found"
RG = sorted(cands)[-1]
print("regen dir:", RG)

E = np.load(f"{DEC}/chr10_sv_genotypes.npz"); E = E[list(E.keys())[0]]
R = np.load(f"{RG}/chr10_sv_genotypes.npz");   R = R[list(R.keys())[0]]
print("E", E.shape, "sum", int(E.sum()), "| R", R.shape, "sum", int(R.sum()))
assert E.shape == R.shape

dcol = (E != R).sum(axis=1)
naff = int((dcol > 0).sum())
print(f"columns with ANY diff: {naff}/{E.shape[0]}   total diff cells: {int(dcol.sum())}   max per-col: {int(dcol.max())}")
print("affected column indices:", np.where(dcol>0)[0].tolist())

csE = E.sum(axis=1); csR = R.sum(axis=1)
print("column-sum identical:", int((csE==csR).sum()), "/", E.shape[0])
print("columns where colsum differs:", np.where(csE!=csR)[0].tolist())

sv = pd.read_csv(f"{RG}/chr10_sv_variants.tsv", sep="\t")
print("regen sv_variants rows:", len(sv), "== npz cols:", len(sv)==R.shape[0])
print("cols:", list(sv.columns))
print("sv.ac == R colsum (all):", bool((sv['ac'].values == csR).all()))
aff = np.where(dcol>0)[0]
if len(aff):
    print("\n-- drifted columns: coords + ac(E) vs ac(R) --")
    show = sv.iloc[aff][['chrom','pos','end','svtype','svlen','id','sv_category']].copy()
    show['ac_existing']=csE[aff]; show['ac_regen']=csR[aff]; show['ndiff']=dcol[aff]
    print(show.to_string(index=False, max_colwidth=40))
