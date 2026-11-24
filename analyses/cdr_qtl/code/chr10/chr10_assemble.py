"""Assemble the certified chr10_sv_variants.tsv, row-aligned to the AUTHORITATIVE npz."""
import numpy as np, pandas as pd, glob, sys
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
P=CENSV_ROOT + ""
PARAMS="mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
DEC=f"{P}/sv_extraction/decompose/{PARAMS}/chr10"
RG=sorted(glob.glob(CENSV_SCRATCH + "/chr10_recover_*/dec_cosine_nj/chr10"))[-1]

E=np.load(f"{DEC}/chr10_sv_genotypes.npz"); E=E[list(E.keys())[0]]
R=np.load(f"{RG}/chr10_sv_genotypes.npz");   R=R[list(R.keys())[0]]
SV=pd.read_csv(f"{RG}/chr10_sv_variants.tsv", sep="\t")
assert E.shape==R.shape==(len(SV),196), (E.shape,R.shape,len(SV))

diff_cols=np.where((E!=R).any(axis=1))[0]
print("affected cols:", diff_cols.tolist())

import itertools
from collections import defaultdict
src=np.arange(E.shape[0])
pos=SV['pos'].values
groups=defaultdict(list)
for i in diff_cols: groups[int(pos[i])].append(int(i))
unresolved=[]
for p,cols in groups.items():
    best=None
    for perm in itertools.permutations(cols):
        if all(np.array_equal(E[cols[k],:], R[perm[k],:]) for k in range(len(cols))):
            best=perm; break
    if best is None:
        unresolved+=cols
        print(f"  UNRESOLVED group pos={p} cols={cols}")
    else:
        for k in range(len(cols)): src[cols[k]]=best[k]
        if tuple(cols)!=best: print(f"  resolved pos={p}: {cols} <- regen {list(best)}")

ok_perm = sorted(src.tolist())==list(range(E.shape[0]))
mism=[i for i in range(E.shape[0]) if not np.array_equal(E[i,:], R[src[i],:])]
print("mapping is permutation:", ok_perm, "| residual mismatches:", len(mism), "| unresolved:", unresolved)

final=SV.iloc[src].reset_index(drop=True).copy()
an=2*196
ac=E.sum(axis=1).astype(int)
final['ac']=ac; final['an']=an
final['af']=np.round(ac/an,6); final['maf']=np.round(np.minimum(ac/an,1-ac/an),6)
final['chrom']='chr10'; final['n_miss']=0

assert bool((final['ac'].values==E.sum(axis=1)).all()), "ac mismatch vs authoritative npz"
assert len(final)==7304 and final['pos'].notna().all()
print("final rows:", len(final), "| ac==E.colsum:", bool((final['ac'].values==E.sum(axis=1)).all()))
print("svtype counts:", final['svtype'].value_counts().to_dict())

if len(mism)==0 and ok_perm and not unresolved:
    out=f"{DEC}/chr10_sv_variants.tsv"
    final.to_csv(out, sep="\t", index=False)
    print("WROTE", out)
    pv=f"{RG}/chr10_variants.tsv"
    import os
    if os.path.exists(pv):
        import shutil; shutil.copy(pv, f"{DEC}/chr10_variants.tsv"); print("copied parent chr10_variants.tsv")
    print("CERTIFIED: chr10_sv_variants.tsv row-aligned to authoritative npz (exact).")
else:
    print("NOT WRITING — alignment not fully certified; investigate.")
