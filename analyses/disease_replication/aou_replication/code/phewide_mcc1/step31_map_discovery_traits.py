#!/usr/bin/env python3
import os, csv, re
from collections import Counter

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
PH = f"{BASE}/phewide_all_handoff"
CATL = f"{BASE}/aou_replication/results_vds/step21_catalog/loinc_catalog.tsv"
XINFO = f"{BASE}/aou_replication/results_vds/step21_catalog/phecodeX/phecodeX_info.tsv"
OUT = f"{BASE}/aou_replication/results_vds/step31_phewide"; os.makedirs(OUT, exist_ok=True)

CURATED_LOINC = {
    "50": "8302-2",
    "23104": "39156-5",
    "21001": "39156-5",
    "21002": "29463-7",
    "30080": "777-3",
    "30000": "6690-2",
    "30010": "789-8",
    "30120": "731-0",
    "30140": "751-8",
    "30130": "742-7",
    "30760": "2085-9",
    "30780": "2089-1",
    "30690": "2093-3",
    "30870": "2571-8",
    "30740": "2345-7",
    "30750": "4548-4",
    "30700": "2160-0",
    "30610": "1742-6",
    "30620": "1920-8",
    "30840": "1975-2",
    "30600": "1751-7",
    "100024": "17861-6",
    "30890": "1989-3",
    "HEIGHT_IRN": "8302-2",
    "BMI_IRN": "39156-5",
}
STOP = set("in of and the a to blood serum plasma level count number volume mass fraction "
           "measurement quantile left right total by no nos".split())


def toks(s):
    return set(w for w in re.split(r"[^a-z0-9]+", s.lower()) if w and w not in STOP)


def best_match(name, cand_tokens):
    q = toks(name)
    if not q:
        return None, 0.0
    best, bs = None, 0.0
    for cid, ct in cand_tokens:
        inter = len(q & ct)
        if not inter:
            continue
        j = inter / len(q | ct)
        if j > bs:
            bs, best = j, cid
    return best, bs


def main():
    tl = {r["trait_code"]: (r["trait_name"], r["source_code"])
          for r in csv.DictReader(open(f"{PH}/trait_list.tsv"), delimiter="\t")}
    pair_traits = sorted({(r["study"], r["trait"]) for r in csv.DictReader(open(f"{PH}/sv_trait_pairs.tsv"), delimiter="\t")})

    loinc = [(r["loinc"], r["name"]) for r in csv.DictReader(open(CATL), delimiter="\t")]
    loinc_tok = [(lid, toks(nm)) for lid, nm in loinc]
    loinc_name = dict(loinc)
    phe = [(r["phecode"], r["label"]) for r in csv.DictReader(open(XINFO), delimiter="\t")]
    phe_tok = [(pid, toks(lab)) for pid, lab in phe]
    phe_name = dict(phe)

    rows = []
    for study, trait in pair_traits:
        name = tl.get(trait, (trait, study))[0]
        field = re.sub(r"_irnt$", "", trait)
        aou_type = aou_id = aou_name = ""; score = 0.0; how = ""
        fg_quant = study == "finngen" and re.search(r"height|body-mass|\bbmi\b", name, re.I)
        if study == "neale" or fg_quant:
            if field in CURATED_LOINC:
                aou_id = CURATED_LOINC[field]; score = 1.0; how = "curated"
            else:
                aou_id, score = best_match(name, loinc_tok)
            if aou_id:
                aou_type = "loinc"; aou_name = loinc_name.get(aou_id, "")
        else:
            aou_id, score = best_match(name, phe_tok)
            if aou_id:
                aou_type = "phecode"; aou_name = phe_name.get(aou_id, "")
        mapped = bool(aou_id) and (how == "curated" or score >= 0.34)
        rows.append([study, trait, name, aou_type if mapped else "", aou_id if mapped else "",
                     aou_name if mapped else "", f"{score:.2f}", how or ("name" if mapped else "unmapped")])

    with open(f"{OUT}/trait_crosswalk.tsv", "w") as f:
        f.write("study\ttrait_code\ttrait_name\taou_type\taou_id\taou_name\tmatch_score\tmethod\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    cov = Counter()
    for r in rows:
        key = (r[0], "mapped" if r[3] else "unmapped")
        cov[key] += 1
    print("=== trait crosswalk coverage (distinct pair-traits) ===")
    for study in ("neale", "finngen"):
        m = cov[(study, "mapped")]; u = cov[(study, "unmapped")]
        print(f"  {study:8}: mapped {m} / {m+u}", flush=True)
    print("\n  sample mapped:")
    for r in rows:
        if r[3] and r[7] != "curated":
            print(f"    {r[1]:14} {r[2][:34]:34} -> {r[4]:9} {r[5][:34]} (j={r[6]})")
    print("\n  sample UNMAPPED (no AoU equivalent / low score):")
    for r in [x for x in rows if not x[3]][:12]:
        print(f"    {r[0]:8} {r[1]:14} {r[2][:44]}")
    print(f"\nwrote {OUT}/trait_crosswalk.tsv", flush=True)


if __name__ == "__main__":
    main()
