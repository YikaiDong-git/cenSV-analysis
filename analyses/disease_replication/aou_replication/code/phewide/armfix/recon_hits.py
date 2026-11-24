#!/usr/bin/env python3
"""Reconcile the disease hit table old (pre-armfix) vs new (post-armfix), keyed by (sv_locus, phecode)."""
import os, csv
from collections import OrderedDict

RV = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds")
OLD = f"{RV}/archive_2026-07-11_pre_armfix/step30_phewide/phewide_disease_hits_enriched.tsv"
NEW = f"{RV}/step30_phewide/phewide_disease_hits_enriched.tsv"
AFF = {"chr1", "chr9", "chr16"}


def load(p):
    d = OrderedDict()
    for r in csv.DictReader(open(p), delimiter="\t"):
        d[(r["sv_locus"], r["phecode"])] = r
    return d


def fnum(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return float("nan")


def main():
    old, new = load(OLD), load(NEW)
    ok, nk = set(old), set(new)
    retained = sorted(ok & nk); dropped = sorted(ok - nk); added = sorted(nk - ok)
    print(f"OLD hits: {len(old)} ({len({k[0] for k in old})} SVs) | NEW hits: {len(new)} ({len({k[0] for k in new})} SVs)")
    print(f"  retained: {len(retained)} | dropped: {len(dropped)} | added: {len(added)}\n")

    def chrom(k): return k[0].split("_")[0]

    print("=== UNAFFECTED-chrom hits (chr10/11/17 — must be STABLE) ===")
    stable_bad = []
    for k in retained:
        if chrom(k) in AFF:
            continue
        po, pn = fnum(old[k]["p_score"]), fnum(new[k]["p_score"])
        rel = abs(po - pn) / max(abs(po), 1e-300)
        flag = "  <-!! MOVED" if rel > 1e-2 else ""
        if flag:
            stable_bad.append(k)
        print(f"  {k[0][:30]:30} {k[1]:9} p {po:.2e}->{pn:.2e} q {old[k]['q_bh_tier']}->{new[k]['q_bh_tier']}{flag}")
    for k in dropped + added:
        if chrom(k) not in AFF:
            print(f"  {'DROPPED' if k in dropped else 'ADDED':7} {k[0][:30]:30} {k[1]:9}  <-!! unexpected off-target")
    print(f"  unaffected stability: {'OK' if not stable_bad else 'CHECK '+str(stable_bad)}\n")

    print("=== AFFECTED-chrom hits (chr1/9/16 — expected to CHANGE) ===")
    print(f"  retained (with delta):")
    for k in retained:
        if chrom(k) not in AFF:
            continue
        po, pn = fnum(old[k]["p_score"]), fnum(new[k]["p_score"])
        print(f"    {k[0][:30]:30} {k[1]:9} p {po:.2e}->{pn:.2e} q {old[k]['q_bh_tier']}->{new[k]['q_bh_tier']} "
              f"firth_conf {old[k]['firth_confirmed']}->{new[k]['firth_confirmed']}")
    print(f"  DROPPED (were hits, no longer): {[f'{k[0].split(chr(95))[1]}/{k[1]}' for k in dropped if chrom(k) in AFF]}")
    print(f"  ADDED (new hits): {[f'{k[0].split(chr(95))[1]}/{k[1]}' for k in added if chrom(k) in AFF]}")


if __name__ == "__main__":
    main()
