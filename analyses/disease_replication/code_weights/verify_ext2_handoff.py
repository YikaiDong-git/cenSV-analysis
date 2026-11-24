#!/usr/bin/env python3
"""Full-scale verification of the ext2 handoff before it is shipped."""
import collections
import csv
import glob
import os
import sys
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

B = CENSV_ROOT + ""
HAND = B + "/analyses/disease_replication/phewide_all_handoff"
EXT2 = HAND + "/full_11129_ext2"
W = B + "/analyses/sv_catalog/output/shared/weights_decompose_two_arm"

fail = []


def check(label, ok, detail=""):
    print("  %-62s %s%s" % (label, "OK" if ok else "FAIL", ("  " + detail) if detail else ""))
    if not ok:
        fail.append(label)


ext2 = collections.defaultdict(dict)
ext2_txt = collections.defaultdict(dict)
dup = 0
a1_bad = 0
nrows = 0
for fp in sorted(glob.glob(EXT2 + "/sv_weights_long_ext2_chr*.tsv*")):
    with open(fp) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            nrows += 1
            k = (int(r["pos"]), r["ref"], r["alt"])
            if r["a1_effect"] != r["alt"]:
                a1_bad += 1
            if k in ext2[r["sv_locus"]]:
                dup += 1
            ext2[r["sv_locus"]][k] = float(r["beta_std"])
            ext2_txt[r["sv_locus"]][k] = r["beta_std"]
print("ext2: %d predictors, %d weight rows" % (len(ext2), nrows))

shipped = collections.defaultdict(dict)
shipped_txt = collections.defaultdict(dict)
for pat in (HAND + "/sv_weights_long_chr*.tsv",
            HAND + "/discovery_null_ext/sv_weights_long_ext_chr*.tsv"):
    for fp in sorted(glob.glob(pat)):
        with open(fp) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                kk = (int(r["pos"]), r["ref"], r["alt"])
                shipped[r["sv_locus"]][kk] = float(r["beta_std"])
                shipped_txt[r["sv_locus"]][kk] = r["beta_std"]
print("shipped: %d predictors" % len(shipped))

print("\n=== assertions ===")
check("ext2 covers all 11,129 candidate predictors", len(ext2) == 11129,
      "got %d" % len(ext2))
check("a1_effect == alt on every row", a1_bad == 0, "%d violations" % a1_bad)
check("no duplicate (sv_locus,pos,ref,alt) rows", dup == 0, "%d dupes" % dup)

missing_sv = [s for s in shipped if s not in ext2]
check("every shipped predictor present in ext2", not missing_sv,
      "%d absent" % len(missing_sv))

beta_bad = snp_dropped = 0
worst_rel = 0.0
worst_ex = None
n_cmp = 0
for sv, sh in shipped_txt.items():
    e = ext2_txt.get(sv, {})
    for k, v in sh.items():
        if k not in e:
            snp_dropped += 1
            continue
        n_cmp += 1
        if e[k] != v:
            beta_bad += 1
            a, b = float(v), float(e[k])
            rel = abs(a - b) / max(abs(a), 1e-300)
            if rel > worst_rel:
                worst_rel, worst_ex = rel, (sv, k, v, e[k])
check("the byte-identity comparison actually ran (not vacuous)", n_cmp > 100000,
      "compared %d SNPs" % n_cmp)
check("no shipped SNP dropped by ext2 (superset property)", snp_dropped == 0,
      "%d dropped" % snp_dropped)
check("beta agreement: differences confined to the last digit, <1e-9 relative",
      worst_rel < 1e-9,
      "%d of %d differ (%.4f%%), max relative %.2e%s"
      % (beta_bad, n_cmp, 100.0 * beta_bad / max(n_cmp, 1), worst_rel,
         ("  e.g. %s old=%s new=%s" % (worst_ex[0], worst_ex[2], worst_ex[3]))
         if worst_ex else ""))

crlf = []
for fp in sorted(glob.glob(EXT2 + "/*.tsv")):
    with open(fp, "rb") as fh:
        if bytes((13, 10)) in fh.read(200000):
            crlf.append(os.path.basename(fp))
check("LF line endings, matching the existing bundles", not crlf,
      "CRLF in %d file(s) e.g. %s" % (len(crlf), crlf[:2]))

meta = list(csv.DictReader(open(EXT2 + "/sv_metadata_ext2.tsv"), delimiter="\t"))
check("metadata rows == 11,129", len(meta) == 11129, "got %d" % len(meta))
check("metadata sv_locus set == weight-file sv_locus set",
      {m["sv_locus"] for m in meta} == set(ext2))
check("metadata n_weight_snps matches actual row counts",
      all(int(m["n_weight_snps"]) == len(ext2[m["sv_locus"]]) for m in meta))
t = collections.Counter(m["tier"] for m in meta)
check("tier counts (79 / 1,684 / 2,847 / 6,519)",
      [t["A_chrX_reference_only"], t["B_ref_snps_present"], t["C_zero_ref_match"],
       t["admitted_already_shipped"]] == [79, 1684, 2847, 6519], str(dict(t)))
newn = sum(1 for m in meta if m["already_shipped"] == "no")
check("new-to-AoU predictors == 4,616", newn == 4616, "got %d" % newn)

c14 = [m for m in meta if m["chrom"] == "chr14"]
cx = [m for m in meta if m["chrom"] == "chrX"]
check("chr14 now has weight vectors (was absent from AoU)", len(c14) == 239,
      "%d predictors" % len(c14))
check("chrX now has weight vectors (was absent from AoU)", len(cx) == 79,
      "%d predictors" % len(cx))

disk = {os.path.basename(f) for f in glob.glob(W + "/chr*/sv*/*_weights.tsv")}
refd = set()
for m in sorted(glob.glob(W + "/chr*/chr*_MANIFEST.tsv")):
    with open(m) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            refd.add(os.path.basename(r["weight_path"]))
orphan = sorted(disk - refd)
print("\n=== weight files on disk not referenced by any manifest: %d ===" % len(orphan))
for o in orphan[:10]:
    print("  %s" % o)
check("orphan files are accounted for (expected 4, none are candidates)",
      len(orphan) == 4, "%d orphans" % len(orphan))

print("\n%s" % ("ALL CHECKS PASSED - safe to ship" if not fail
                else "FAILED: %d\n  - %s" % (len(fail), "\n  - ".join(fail))))
sys.exit(1 if fail else 0)
