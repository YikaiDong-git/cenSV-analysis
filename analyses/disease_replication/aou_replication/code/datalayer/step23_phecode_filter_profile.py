#!/usr/bin/env python3
import os, csv
from collections import defaultdict

HOME = os.path.expanduser("~")
PRS = f"{HOME}/workspace/working/PRS/resources/phewas"
CAT = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/"
                         "disease_replication/aou_replication/results_vds/step21_catalog")
OUT = CAT

SEX_GROUPS = {"pregnancy complications", "genitourinary"}


def main():
    npeople = {}
    for r in csv.DictReader(open(f"{CAT}/icd_code_catalog.tsv"), delimiter="\t"):
        npeople[(r["vocabulary"], r["code"])] = int(r["n_people"])
    codes = defaultdict(list)
    for r in csv.DictReader(open(f"{PRS}/phecode_map_expanded.tsv"), delimiter="\t"):
        codes[r["phecode"]].append((r["vocabulary_id"], r["concept_code"]))
    info = {}
    for r in csv.DictReader(open(f"{PRS}/phecode_info.tsv"), delimiter="\t"):
        info[r["phecode"]] = (r["description"], r.get("group", ""))

    rows = []
    for pc, cc in codes.items():
        ppl = [npeople.get(k, 0) for k in cc]
        floor = max(ppl) if ppl else 0
        upper = sum(ppl)
        n_codes = len(cc); n_codes_present = sum(1 for x in ppl if x > 0)
        desc, group = info.get(pc, ("", ""))
        rows.append((pc, desc, group, n_codes, n_codes_present, floor, upper))

    rows.sort(key=lambda r: -r[5])
    with open(f"{OUT}/phecode_aou_profile.tsv", "w") as fo:
        fo.write("phecode\tdescription\tgroup\tn_icd_codes\tn_codes_in_aou\tcases_floor\tcases_upper\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")

    total = len(rows)
    print(f"=== {total} phecodes profiled against AoU (cases_floor = max single-code person-count) ===\n")

    print("PERSPECTIVE 1 — POWER (survival by case-count floor):")
    for thr in (0, 20, 50, 100, 200, 500, 1000, 5000, 10000):
        n = sum(1 for r in rows if r[5] >= thr)
        print(f"  cases_floor >= {thr:6d}:  {n:5d} phecodes")
    print()

    print("PERSPECTIVE 2 — COVERAGE (phecodes with NO AoU code at all = untestable):")
    zero = sum(1 for r in rows if r[5] == 0)
    print(f"  {zero} phecodes have 0 AoU-present ICD codes (drop outright); {total-zero} have >=1.\n")

    print("PERSPECTIVE 3 — CATEGORY (phecodes with cases_floor>=200, by group):")
    bycat = defaultdict(int); bycat_all = defaultdict(int)
    for r in rows:
        bycat_all[r[2]] += 1
        if r[5] >= 200:
            bycat[r[2]] += 1
    for g in sorted(bycat_all, key=lambda x: -bycat[x]):
        print(f"  {g[:34]:34}  {bycat[g]:4d} / {bycat_all[g]:4d}  (>=200 cases / total)")
    print()

    print("PERSPECTIVE 4 — SEX-SPECIFIC groups (need per-sex handling):")
    for g in sorted(SEX_GROUPS):
        n = sum(1 for r in rows if r[2] == g and r[5] >= 200)
        print(f"  {g}: {n} phecodes >=200 cases")
    print()

    print("TOP 12 best-powered phecodes:")
    for r in rows[:12]:
        print(f"  {r[0]:8} {r[1][:38]:38} {r[2][:20]:20} floor={r[5]:7d}")
    print(f"\nwrote {OUT}/phecode_aou_profile.tsv")


if __name__ == "__main__":
    main()
