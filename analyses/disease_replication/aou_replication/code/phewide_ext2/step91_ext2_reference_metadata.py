#!/usr/bin/env python3
"""Step 91 (ext2 arm) — ARM-NEUTRAL REFERENCE METADATA, derived inside the ext2 tree."""
import os
import re
import sys
import csv
import json
import hashlib
import argparse
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from frozen_tree_sig import tree_sig                     # noqa: E402

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
PH = f"{BASE}/phewide_all_handoff"
RV = f"{BASE}/aou_replication/results_vds"
CARRY = f"{RV}/phewide_carryover/inputs"
CATL = f"{RV}/step21_catalog/loinc_catalog.tsv"
XINFO = f"{RV}/step21_catalog/phecodeX/phecodeX_info.tsv"
OUT = f"{RV}/step30_phewide_ext2/reference"

PUB_PSUM = f"{RV}/step30_phewide/phecodeX_pheno_summary.tsv"
PUB_QSUM = f"{RV}/step30_phewide/quant_pheno_summary.tsv"
PUB_XW = f"{RV}/step31_phewide/trait_crosswalk.tsv"

FROZEN = ["step30_phewide", "step28_phewide", "phewide_carryover", "step30_phewide_mcc1", "step31_phewide"]
CHUNK = 40_000

LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def guard_write(path):
    assert not os.path.islink(path), f"REFUSING to write through symlink: {path}"
    d = os.path.dirname(path)
    assert not os.path.islink(d), f"REFUSING to write into symlinked dir: {d}"


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def build_phecode_summary():
    info = {r["phecode"]: r for r in csv.DictReader(open(XINFO), delimiter="\t")}
    hdr = open(f"{CARRY}/disease.pheno").readline().rstrip("\n").split("\t")
    phecodes = [c for c in hdr if c not in ("FID", "IID")]
    say(f"  disease.pheno: {len(phecodes)} phecode columns; phecodeX_info: {len(info)} phecodes")

    cases = Counter()
    nrow = 0
    for ch in pd.read_csv(f"{CARRY}/disease.pheno", sep="\t", usecols=phecodes,
                          dtype=np.float32, na_values=["NA", "-9", ""], chunksize=CHUNK):
        nrow += len(ch)
        s = (ch == 1).sum()
        for k, v in s.items():
            cases[k] += int(v)
    say(f"  scanned {nrow:,} cohort rows")

    rows = []
    missing = []
    for ph in phecodes:
        rec = info.get(ph)
        if rec is None:
            missing.append(ph)
            rows.append([ph, "", "", str(cases[ph])])
        else:
            rows.append([ph, rec["label"], rec["sex"], str(cases[ph])])
    if missing:
        say(f"  ⚠️ {len(missing)} phecode(s) absent from phecodeX_info: {missing[:5]}")
    return ["phecode", "label", "sex", "cases_analysed"], rows


def build_quant_summary():
    cat = {r["loinc"]: r["name"] for r in csv.DictReader(open(CATL), delimiter="\t")}
    hdr = open(f"{CARRY}/quant.pheno").readline().rstrip("\n").split("\t")
    loincs = [c for c in hdr if c not in ("FID", "IID")]
    say(f"  quant.pheno: {len(loincs)} LOINC columns; loinc_catalog: {len(cat)} LOINCs")

    npeople = Counter()
    nrow = 0
    for ch in pd.read_csv(f"{CARRY}/quant.pheno", sep="\t", usecols=loincs,
                          dtype=np.float32, na_values=["NA", "-9", ""], chunksize=CHUNK):
        nrow += len(ch)
        s = ch.notna().sum()
        for k, v in s.items():
            npeople[k] += int(v)
    say(f"  scanned {nrow:,} cohort rows")

    rows = [[lo, cat.get(lo, ""), str(npeople[lo])] for lo in loincs]
    absent = [lo for lo in loincs if lo not in cat]
    if absent:
        say(f"  ⚠️ {len(absent)} LOINC(s) absent from loinc_catalog: {absent[:5]}")
    return ["loinc", "name", "n_people"], rows


CURATED_LOINC = {
    "50": "8302-2", "23104": "39156-5", "21001": "39156-5", "21002": "29463-7",
    "30080": "777-3", "30000": "6690-2", "30010": "789-8", "30120": "731-0",
    "30140": "751-8", "30130": "742-7", "30760": "2085-9", "30780": "2089-1",
    "30690": "2093-3", "30870": "2571-8", "30740": "2345-7", "30750": "4548-4",
    "30700": "2160-0", "30610": "1742-6", "30620": "1920-8", "30840": "1975-2",
    "30600": "1751-7", "100024": "17861-6", "30890": "1989-3",
    "HEIGHT_IRN": "8302-2", "BMI_IRN": "39156-5",
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


def build_crosswalk():
    tl = {r["trait_code"]: (r["trait_name"], r["source_code"])
          for r in csv.DictReader(open(f"{PH}/trait_list.tsv"), delimiter="\t")}
    pair_traits = sorted({(r["study"], r["trait"])
                          for r in csv.DictReader(open(f"{PH}/sv_trait_pairs.tsv"), delimiter="\t")})
    loinc = [(r["loinc"], r["name"]) for r in csv.DictReader(open(CATL), delimiter="\t")]
    loinc_tok = [(lid, toks(nm)) for lid, nm in loinc]
    loinc_name = dict(loinc)
    phe = [(r["phecode"], r["label"]) for r in csv.DictReader(open(XINFO), delimiter="\t")]
    phe_tok = [(pid, toks(lab)) for pid, lab in phe]
    phe_name = dict(phe)
    say(f"  {len(pair_traits)} distinct (study, trait) pairs from sv_trait_pairs.tsv")

    rows = []
    for study, trait in pair_traits:
        name = tl.get(trait, (trait, study))[0]
        field = re.sub(r"_irnt$", "", trait)
        aou_type = aou_id = aou_name = ""
        score = 0.0
        how = ""
        fg_quant = study == "finngen" and re.search(r"height|body-mass|\bbmi\b", name, re.I)
        if study == "neale" or fg_quant:
            if field in CURATED_LOINC:
                aou_id = CURATED_LOINC[field]
                score = 1.0
                how = "curated"
            else:
                aou_id, score = best_match(name, loinc_tok)
            if aou_id:
                aou_type = "loinc"
                aou_name = loinc_name.get(aou_id, "")
        else:
            aou_id, score = best_match(name, phe_tok)
            if aou_id:
                aou_type = "phecode"
                aou_name = phe_name.get(aou_id, "")
        mapped = bool(aou_id) and (how == "curated" or score >= 0.34)
        rows.append([study, trait, name, aou_type if mapped else "", aou_id if mapped else "",
                     aou_name if mapped else "", f"{score:.2f}", how or ("name" if mapped else "unmapped")])
    return ["study", "trait_code", "trait_name", "aou_type", "aou_id", "aou_name",
            "match_score", "method"], rows


def p3_coverage():
    pair = {}
    for r in csv.DictReader(open(f"{PH}/sv_trait_pairs.tsv"), delimiter="\t"):
        pair.setdefault(r["sv_locus"], set()).add((r["study"], r["trait"]))
    order = f"{RV}/step30_phewide_ext2/disease_matrices/disease_svorder.txt"
    scored = {ln.split("\t")[0] for ln in open(order) if ln.strip()}
    testable = sorted(set(pair) & scored)
    rows = [["n_pair_svs", str(len(pair))],
            ["n_pair_trait_combinations", str(sum(len(v) for v in pair.values()))],
            ["n_scored_svs_ext2", str(len(scored))],
            ["n_pair_svs_testable_ext2", str(len(testable))],
            ["n_pair_svs_not_scored_ext2", str(len(set(pair) - scored))]]
    say(f"  discovery pair SVs {len(pair)} · testable in this arm {len(testable)} "
        f"({100.0 * len(testable) / max(1, len(pair)):.1f}%)")
    return ["metric", "value"], rows


def write_tsv(path, hdr, rows):
    guard_write(path)
    with open(path, "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def gate_phecode_summary(built):
    B = list(csv.DictReader(open(built), delimiter="\t"))
    ok = True

    if os.path.exists(PUB_PSUM):
        P = list(csv.DictReader(open(PUB_PSUM), delimiter="\t"))
        b3 = [(r["phecode"], r["label"], r["sex"]) for r in B]
        p3 = [(r["phecode"], r["label"], r["sex"]) for r in P]
        g1 = b3 == p3
        say(f"    G1 (phecode,label,sex) vs published : {'PASS — identical' if g1 else 'FAIL'}")
        if not g1:
            ok = False
            diff = [(x, y) for x, y in zip(b3, p3) if x != y][:5]
            for x, y in diff:
                say(f"       rebuilt {x}  vs published {y}")
    else:
        say("    G1 SKIPPED — no published copy")
        P = []

    meta = f"{RV}/step30_phewide_ext2/disease_matrices/disease_meta.tsv"
    if os.path.exists(meta):
        E = {r["phecode"]: int(r["cases"]) for r in csv.DictReader(open(meta), delimiter="\t")}
        bad = [(r["phecode"], int(r["cases_analysed"]), E[r["phecode"]])
               for r in B if r["phecode"] in E and int(r["cases_analysed"]) != E[r["phecode"]]]
        cov = sum(1 for r in B if r["phecode"] in E)
        g2 = not bad and cov == len(B)
        say(f"    G2 cases_analysed vs ext2 engine    : "
            f"{'PASS — ' + str(cov) + '/' + str(len(B)) + ' exact' if g2 else 'FAIL (' + str(len(bad)) + ' mismatched)'}")
        if not g2:
            ok = False
            for ph, x, y in bad[:5]:
                say(f"       {ph}: rebuilt {x} vs engine {y}")
    else:
        say("    G2 SKIPPED — no ext2 disease_meta.tsv")
        E = {}

    if P and E:
        rows = []
        ratios = []
        for r in P:
            ph = r["phecode"]
            if ph not in E:
                continue
            pub_c, ana = int(r["cases_mcc2"]), E[ph]
            rows.append([ph, str(pub_c), str(ana), f"{pub_c / ana:.4f}" if ana else "NA"])
            if ana:
                ratios.append(pub_c / ana)
        rp = f"{OUT}/cases_reconciliation.tsv"
        write_tsv(rp, ["phecode", "published_cases_mcc2", "cases_analysed", "ratio"], rows)
        a = np.array(ratios)
        say(f"    reconciliation: published/analysed  median {np.median(a):.3f}  "
            f"p05 {np.percentile(a, 5):.3f}  p95 {np.percentile(a, 95):.3f}  "
            f"| overstated >20% in {int((a > 1.2).sum())}/{len(a)} phecodes")
        say(f"    -> {os.path.basename(rp)}")
    return ok


def gate(name, built, published):
    if not os.path.exists(published):
        say(f"  {name:<28} GATE SKIPPED — no published copy at {published}")
        return None
    a = md5_file(built)
    b = md5_file(published)
    if a == b:
        say(f"  {name:<28} GATE PASS   md5 {a}  (identical to the published copy)")
        return True
    la = open(built).read().splitlines()
    lb = open(published).read().splitlines()
    say(f"  {name:<28} GATE FAIL   rebuilt {len(la)} lines / published {len(lb)} lines")
    say(f"      rebuilt md5 {a}")
    say(f"      published   {b}")
    shown = 0
    for i in range(max(len(la), len(lb))):
        x = la[i] if i < len(la) else "<absent>"
        y = lb[i] if i < len(lb) else "<absent>"
        if x != y:
            say(f"      line {i + 1}: rebuilt   {x[:120]}")
            say(f"      line {i + 1}: published {y[:120]}")
            shown += 1
            if shown >= 10:
                say("      ... (further differences suppressed)")
                break
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-counts", action="store_true",
                    help="crosswalk + coverage only; skip the two multi-GB .pheno scans")
    a = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    say("=" * 100)
    say("STEP 91 — arm-neutral reference metadata, derived inside the ext2 tree")
    say("=" * 100)
    base = {n: tree_sig(f"{RV}/{n}") for n in FROZEN}
    say(f"frozen baseline captured for {len(FROZEN)} trees")
    say()

    results = {}
    verdicts = {}

    if not a.skip_counts:
        say("[1/4] phecodeX_pheno_summary.tsv  <- carryover disease.pheno + phecodeX_info")
        hdr, rows = build_phecode_summary()
        p = f"{OUT}/phecodeX_pheno_summary.tsv"
        write_tsv(p, hdr, rows)
        verdicts["phecodeX_pheno_summary"] = gate_phecode_summary(p)
        results["phecodeX_pheno_summary"] = dict(path=p, n_rows=len(rows), md5=md5_file(p))
        say()

        say("[2/4] quant_pheno_summary.tsv  <- carryover quant.pheno + loinc_catalog")
        hdr, rows = build_quant_summary()
        p = f"{OUT}/quant_pheno_summary.tsv"
        write_tsv(p, hdr, rows)
        verdicts["quant_pheno_summary"] = gate("quant_pheno_summary", p, PUB_QSUM)
        results["quant_pheno_summary"] = dict(path=p, n_rows=len(rows), md5=md5_file(p))
        say()
    else:
        say("[1-2/4] SKIPPED (--skip-counts)")
        say()

    say("[3/4] trait_crosswalk.tsv  <- trait_list + sv_trait_pairs + both catalogs")
    hdr, rows = build_crosswalk()
    p = f"{OUT}/trait_crosswalk.tsv"
    write_tsv(p, hdr, rows)
    verdicts["trait_crosswalk"] = gate("trait_crosswalk", p, PUB_XW)
    cov = Counter()
    for r in rows:
        cov[(r[0], "mapped" if r[3] else "unmapped")] += 1
    for study in ("neale", "finngen"):
        m, u = cov[(study, "mapped")], cov[(study, "unmapped")]
        say(f"    {study:8} mapped {m} / {m + u}")
    results["trait_crosswalk"] = dict(path=p, n_rows=len(rows), md5=md5_file(p),
                                      mapped=sum(1 for r in rows if r[3]))
    say()

    say("[4/4] P3 pair-SV coverage in THIS arm")
    hdr, rows = p3_coverage()
    p = f"{OUT}/p3_pair_coverage.tsv"
    write_tsv(p, hdr, rows)
    results["p3_pair_coverage"] = dict(path=p, **{k: v for k, v in rows})
    say()

    say("=" * 100)
    say("VERDICTS")
    how = {"phecodeX_pheno_summary": "PASS (G1 label/sex identical + G2 counts exact vs engine)"}
    for k, v in verdicts.items():
        say(f"  {k:<28} {how.get(k, 'PASS (byte-identical)') if v else ('SKIPPED' if v is None else 'FAIL')}")
    say("=" * 100)

    changed = [n for n in FROZEN if tree_sig(f"{RV}/{n}")[0] != base[n][0]]
    say(f"frozen-tree integrity: {'0 trees changed — OK' if not changed else 'CHANGED: ' + str(changed)}")
    assert not changed, f"frozen trees modified: {changed}"

    with open(f"{OUT}/step91_summary.json", "w") as f:
        json.dump(dict(tables=results, gates={k: v for k, v in verdicts.items()},
                       frozen_unchanged=True), f, indent=2)
    with open(f"{OUT}/step91_run.log", "w") as f:
        f.write("\n".join(LOG) + "\n")

    failed = [k for k, v in verdicts.items() if v is False]
    if failed:
        say(f"\nSTEP91_GATE_FAILED: {failed}")
        say("A rebuild that does not reproduce the published table means the table encodes something")
        say("about the 6,513 panel we have not accounted for. Reconcile before repointing anything.")
        sys.exit(1)
    say("\nSTEP91_DONE")


if __name__ == "__main__":
    main()
