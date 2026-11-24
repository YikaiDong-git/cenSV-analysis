#!/usr/bin/env python3
import os, csv
from collections import defaultdict

HOME = os.path.expanduser("~")
XCSV = f"{HOME}/.local/lib/python3.11/site-packages/phetk/phecode/phecodeX.csv"
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
CAT = f"{BASE}/aou_replication/results_vds/step21_catalog"
OUT = f"{CAT}/phecodeX"; os.makedirs(OUT, exist_ok=True)
FLAG2VOCAB = {"9": "ICD9CM", "10": "ICD10CM"}


def main():
    npeople = {}
    for r in csv.DictReader(open(f"{CAT}/icd_code_catalog.tsv"), delimiter="\t"):
        npeople[(r["vocabulary"], r["code"])] = int(r["n_people"])

    codes = defaultdict(list)
    info = {}
    n_unmapped_flag = 0
    with open(XCSV) as f:
        for r in csv.DictReader(f):
            vocab = FLAG2VOCAB.get(str(r["flag"]).strip())
            if vocab is None:
                n_unmapped_flag += 1
                continue
            pc = r["phecode"]
            codes[pc].append((vocab, r["ICD"]))
            if pc not in info:
                info[pc] = (r["phecode_string"], r.get("phecode_category", ""), (r.get("sex") or "Both").strip())

    with open(f"{OUT}/phecodeX_map_expanded.tsv", "w") as fo:
        fo.write("vocabulary_id\tconcept_code\tphecode\n")
        for pc, cc in codes.items():
            for vocab, code in cc:
                fo.write(f"{vocab}\t{code}\t{pc}\n")
    with open(f"{OUT}/phecodeX_info.tsv", "w") as fo:
        fo.write("phecode\tlabel\tcategory\tsex\tn_icd_codes\n")
        for pc, cc in codes.items():
            lab, cat, sx = info[pc]
            fo.write(f"{pc}\t{lab}\t{cat}\t{sx}\t{len(cc)}\n")

    rows = []
    for pc, cc in codes.items():
        ppl = [npeople.get(k, 0) for k in cc]
        floor = max(ppl) if ppl else 0
        upper = sum(ppl)
        n_present = sum(1 for x in ppl if x > 0)
        lab, cat, sx = info[pc]
        rows.append((pc, lab, cat, sx, len(cc), n_present, floor, upper))
    rows.sort(key=lambda r: -r[6])
    with open(f"{OUT}/phecodeX_aou_profile.tsv", "w") as fo:
        fo.write("phecode\tlabel\tcategory\tsex\tn_icd_codes\tn_codes_in_aou\tcases_floor\tcases_upper\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")

    total = len(rows)
    print(f"=== phecodeX built: {total} phecodes, {sum(len(v) for v in codes.values()):,} ICD rows "
          f"({n_unmapped_flag} non-ICD9/10 rows skipped) ===\n", flush=True)
    print("SURVIVAL by AoU cases_floor (single-code MCC=1 proxy — MCC=2 will be lower):")
    for thr in (0, 20, 50, 100, 200, 500, 1000):
        n = sum(1 for r in rows if r[6] >= thr)
        print(f"  cases_floor >= {thr:6d}:  {n:5d} phecodes", flush=True)
    zero = sum(1 for r in rows if r[6] == 0)
    print(f"\nuntestable (0 AoU-present codes): {zero};  testable: {total - zero}")
    sexsp = sum(1 for r in rows if r[3] in ("Male", "Female"))
    sexsp50 = sum(1 for r in rows if r[3] in ("Male", "Female") and r[6] >= 50)
    print(f"sex-specific phecodes: {sexsp} ({sexsp50} with >=50 cases_floor)\n")
    print("CATEGORY (phecodes with cases_floor>=50, by category):")
    bycat = defaultdict(int); bycat_all = defaultdict(int)
    for r in rows:
        bycat_all[r[2]] += 1
        if r[6] >= 50:
            bycat[r[2]] += 1
    for g in sorted(bycat_all, key=lambda x: -bycat[x]):
        print(f"  {g[:32]:32}  {bycat[g]:4d} / {bycat_all[g]:4d}", flush=True)
    print("\nTOP 10 best-powered phecodeX:")
    for r in rows[:10]:
        print(f"  {r[0]:8} {r[1][:40]:40} {r[2][:16]:16} floor={r[6]:7d}", flush=True)
    print(f"\nwrote {OUT}/{{phecodeX_map_expanded.tsv, phecodeX_info.tsv, phecodeX_aou_profile.tsv}}")


if __name__ == "__main__":
    main()
