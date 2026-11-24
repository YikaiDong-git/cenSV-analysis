#!/usr/bin/env python3
from __future__ import print_function

import argparse
import gzip
import io
import os
import subprocess
import sys
import tempfile

URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"

ADJUSTMENTS = {
    "chr1": ((121700000, 125100000), (121700000, 128100000),
             "q edge moved 3.0 Mb outward, past the start of the 1q12 block"),
    "chr2": ((91800000, 96000000), (91800000, 94200000),
             "q edge moved 1.8 Mb inward, off the SNP-poor band"),
    "chr6": ((58500000, 62600000), (57800000, 63300000),
             "both edges moved 0.7 Mb outward"),
    "chrX": ((58100000, 63800000), (56600000, 62100000),
             "both edges moved inward, 1.5 Mb on p and 1.7 Mb on q"),
}

DROP = {"chrY"}

ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX"]


def fetch(path=None):
    if path:
        raw = open(path, "rb").read()
    else:
        tmp = tempfile.NamedTemporaryFile(suffix=".gz", delete=False)
        tmp.close()
        subprocess.check_call(["curl", "-sSL", "--fail", URL, "-o", tmp.name])
        raw = open(tmp.name, "rb").read()
        os.unlink(tmp.name)
    return gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode("ascii", "replace")


def acen_merge(text):
    out = {}
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 5 or f[4] != "acen":
            continue
        c, s, e = f[0], int(f[1]), int(f[2])
        if c in out:
            out[c] = (min(out[c][0], s), max(out[c][1], e))
        else:
            out[c] = (s, e)
    return out


def build(cytoband_path=None):
    merged = acen_merge(fetch(cytoband_path))
    rows = []
    for c in ORDER:
        if c in DROP or c not in merged:
            continue
        got = merged[c]
        if c in ADJUSTMENTS:
            expect, used, why = ADJUSTMENTS[c]
            if got != expect:
                raise SystemExit(
                    "FAIL %s: UCSC cytoBand now merges to %s, but the adjustment recorded here\n"
                    "was derived from %s. Re-derive the adjustment before trusting the override.\n"
                    "  (%s)" % (c, got, expect, why))
            rows.append((c, used[0], used[1]))
        else:
            rows.append((c, got[0], got[1]))
    n_adj = sum(1 for c, _, _ in rows if c in ADJUSTMENTS)
    if n_adj != len(ADJUSTMENTS):
        raise SystemExit("FAIL: %d of %d adjustments applied" % (n_adj, len(ADJUSTMENTS)))
    return rows


def render(rows):
    return "".join("%s\t%d\t%d\n" % r for r in rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    ap.add_argument("--check", metavar="BED")
    ap.add_argument("--cytoband")
    a = ap.parse_args()

    rows = build(a.cytoband)
    text = render(rows)

    if a.check:
        have = open(a.check).read()
        if have == text:
            print("OK  %s reproduces from cytoBand plus the %d recorded adjustments "
                  "(%d intervals)" % (a.check, len(ADJUSTMENTS), len(rows)))
            return 0
        print("FAIL %s does not match what this script builds" % a.check)
        h = dict((l.split("\t")[0], l) for l in have.splitlines() if l.strip())
        w = dict((l.split("\t")[0], l) for l in text.splitlines() if l.strip())
        for c in sorted(set(h) | set(w), key=lambda x: ORDER.index(x) if x in ORDER else 99):
            if h.get(c) != w.get(c):
                print("  %-6s committed=%-28s built=%s"
                      % (c, h.get(c, "(absent)"), w.get(c, "(absent)")))
        return 1

    if not a.out:
        sys.stdout.write(text)
        return 0
    open(a.out, "w").write(text)
    print("wrote %s: %d intervals, %d adjusted (%s)"
          % (a.out, len(rows), len(ADJUSTMENTS), ", ".join(sorted(ADJUSTMENTS))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
