#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

SRC = (CENSV_ROOT + "/analyses/"
       "cdr_calls/output/calls")
OUT = (CENSV_ROOT + "/analyses/"
       "cdr_definition/output/callsets")
GAP = 150_000


def merge(iv, gap, strict):
    out = []
    for st, en in iv:
        if out:
            d = st - out[-1][1]
            if (d < gap) if strict else (d <= gap):
                out[-1][1] = max(out[-1][1], en)
                continue
        out.append([st, en])
    return out


def load(tagdir):
    by = {}
    p = os.path.join(tagdir, "cdrs.bed")
    if not os.path.exists(p):
        return by
    with open(p) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            by.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    for k in by:
        by[k].sort()
    return by


def main():
    tags = sorted(d for d in glob.glob(SRC + "/*") if os.path.isdir(d))
    assert tags, "no canonical dip directories under %s" % SRC

    os.makedirs(OUT + "/gao_cdrs", exist_ok=True)
    os.makedirs(OUT + "/gao_cdrs_lt", exist_ok=True)

    n_dips = n_le = n_lt = 0
    n_arrays = 0
    ties = 0
    overlapping_pairs = 0
    rows = []

    for tagdir in tags:
        tag = os.path.basename(tagdir)
        by = load(tagdir)
        if not by:
            continue
        os.makedirs(OUT + "/gao_cdrs/%s" % tag, exist_ok=True)
        os.makedirs(OUT + "/gao_cdrs_lt/%s" % tag, exist_ok=True)
        fle = open(OUT + "/gao_cdrs/%s/cdrs.bed" % tag, "w")
        flt = open(OUT + "/gao_cdrs_lt/%s/cdrs.bed" % tag, "w")
        for ctg, iv in sorted(by.items()):
            n_arrays += 1
            for a, b in zip(iv, iv[1:]):
                d = b[0] - a[1]
                if d == GAP:
                    ties += 1
                if d < 0:
                    overlapping_pairs += 1

            mle = merge(iv, GAP, strict=False)
            mlt = merge(iv, GAP, strict=True)

            ident = merge(iv, -1, strict=False)
            assert [list(x) for x in iv] == ident, \
                "%s %s: gap=-1 merge altered the input" % (tag, ctg)
            sat = merge(iv, 10 ** 12, strict=False)
            assert sat == [[iv[0][0], max(e for _, e in iv)]], \
                "%s %s: unbounded merge did not collapse to the full span" % (tag, ctg)
            assert len(mle) <= len(mlt) <= len(iv), \
                "%s %s: merge is not monotone in the rule (%d/%d/%d)" \
                % (tag, ctg, len(mle), len(mlt), len(iv))
            starts = set(s for s, _ in iv)
            ends = set(e for _, e in iv)
            for s, e in mle:
                assert s in starts and e in ends, \
                    "%s %s: merged CDR [%d,%d] has an invented boundary" % (tag, ctg, s, e)
            assert sum(e - s for s, e in mle) >= sum(e - s for s, e in iv), \
                "%s %s: merged CDRs cover fewer bases than the dips they contain" % (tag, ctg)

            for s, e in mle:
                fle.write("%s\t%d\t%d\n" % (ctg, s, e))
            for s, e in mlt:
                flt.write("%s\t%d\t%d\n" % (ctg, s, e))

            n_dips += len(iv)
            n_le += len(mle)
            n_lt += len(mlt)
            rows.append((tag, ctg, len(iv), len(mle), len(mlt),
                         sum(e - s for s, e in iv), sum(e - s for s, e in mle)))
        fle.close()
        flt.close()

    with open(OUT + "/gao_unit_summary.tsv", "w") as f:
        f.write("tag\tcontig\tn_dips\tn_cdr_le\tn_cdr_lt\tdip_bp\tcdr_bp\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")

    print("haplotypes            : %d" % len(tags))
    print("active arrays          : %d" % n_arrays)
    print("dips                   : %d" % n_dips)
    print("Gao CDRs (gap <= 150kb): %d" % n_le)
    print("Gao CDRs (gap <  150kb): %d" % n_lt)
    print("pairs separated by EXACTLY 150000 bp: %d" % ties)
    print("overlapping dip pairs in the input   : %d" % overlapping_pairs)
    assert overlapping_pairs == 0, \
        "the canonical dip set contains overlapping calls; the merge assumes it does not"
    if ties == 0:
        print("[gate] no boundary ties -> the le and lt rules are the SAME partition here")
        assert n_le == n_lt, "no ties but the two rules disagree -- merge() is inconsistent"
    else:
        print("[note] %d boundary ties: le and lt differ by %d CDRs" % (ties, n_lt - n_le))
    print("[gate] all per-array merge assertions PASSED")


if __name__ == "__main__":
    sys.exit(main())
