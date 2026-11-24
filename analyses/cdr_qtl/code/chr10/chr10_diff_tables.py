import pandas as pd, glob, os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
WU=CENSV_ROOT + "/analyses/cdr_qtl"
OLD=f"{WU}/results/_pre_chr10_tables_20260627"
NEW=f"{WU}/results/tables"
for t in sorted(glob.glob(f"{OLD}/T*.tsv")):
    name=os.path.basename(t)
    new=f"{NEW}/{name}"
    if not os.path.exists(new): print(f"{name}: NEW MISSING"); continue
    a=pd.read_csv(t, sep="\t"); b=pd.read_csv(new, sep="\t")
    print(f"\n=== {name}: old {a.shape} -> new {b.shape} ===")
    ccol=next((c for c in b.columns if c.lower() in ("chrom","chromosome","chrom_only","chr")), None)
    if ccol and ccol in a.columns:
        oldc=set(a[ccol].astype(str)); newc=set(b[ccol].astype(str))
        added=sorted(newc-oldc); removed=sorted(oldc-newc)
        print(f"  chrom col '{ccol}': +{added or 'none'}  -{removed or 'none'}")
        if ccol:
            common=sorted(oldc & newc)
            am=a[a[ccol].astype(str).isin(common)].sort_values(ccol).reset_index(drop=True)
            bm=b[b[ccol].astype(str).isin(common)].sort_values(ccol).reset_index(drop=True)
            shared_cols=[c for c in a.columns if c in b.columns]
            try:
                same = am[shared_cols].round(6).equals(bm[shared_cols].round(6)) if am.shape==bm.shape else "shape-diff"
            except Exception as e:
                same=f"compare-err:{e}"
            print(f"  non-chr10 rows identical: {same}")
    else:
        try: print(f"  identical: {a.round(6).equals(b.round(6))}")
        except Exception as e: print(f"  compare-err: {e}")
    if ccol and ccol in b.columns:
        c10=b[b[ccol].astype(str).str.contains('chr10')]
        if len(c10): print(f"  chr10 row(s) now present:\n{c10.to_string(index=False, max_colwidth=22)}")
