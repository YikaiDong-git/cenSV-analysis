#!/usr/bin/env python3
"""Emit AoU-scoreable weight vectors for ALL 11,129 candidate predictors."""
import collections
import csv
import glob
import os
import sys
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BASE = CENSV_ROOT + ""
W = BASE + "/analyses/sv_catalog/output/shared/weights_decompose_two_arm"
ADMIT = (BASE + "/analyses/sv_catalog/output/ref_ukbb/admit/"
         "admit_decompose_cv01.tsv")
HAND = BASE + "/analyses/disease_replication/phewide_all_handoff"
OUT = HAND + "/full_11129_ext2"
MAXROWS = 0

os.makedirs(OUT, exist_ok=True)


def log(m):
    print(m, flush=True)


admit = {}
for r in csv.DictReader(open(ADMIT), delimiter="\t"):
    admit[r["sv_locus"]] = r
log("admit table rows: %d (ADMIT %d)"
    % (len(admit), sum(1 for v in admit.values() if v["status"] == "ADMIT")))


def tier_of(r):
    if r is None:
        return "unknown"
    if r["status"] == "ADMIT":
        return "admitted_already_shipped"
    if r["status"] == "FAIL_chrX_no_LD_ref":
        return "A_chrX_reference_only"
    return ("B_ref_snps_present" if int(r["n_matched"]) >= 1 else "C_zero_ref_match")


shipped_n = collections.Counter()
for pat in (HAND + "/sv_weights_long_chr*.tsv",
            HAND + "/discovery_null_ext/sv_weights_long_ext_chr*.tsv"):
    for fp in sorted(glob.glob(pat)):
        with open(fp) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            col = "sv_locus" if "sv_locus" in rd.fieldnames else rd.fieldnames[0]
            for r in rd:
                shipped_n[r[col]] += 1
log("already-shipped predictors: %d" % len(shipped_n))

disk = {}
for fp in glob.glob(W + "/chr*/sv*/*_weights.tsv"):
    disk.setdefault(os.path.basename(fp), []).append(fp)
dupes = {k: v for k, v in disk.items() if len(v) > 1}
log("weight files on disk: %d (%d basenames, %d duplicated basenames)"
    % (sum(len(v) for v in disk.values()), len(disk), len(dupes)))
for k in list(dupes)[:5]:
    log("  DUPLICATE basename: %s -> %s" % (k, dupes[k]))

manifest = []
for m in sorted(glob.glob(W + "/chr*/chr*_MANIFEST.tsv")):
    with open(m) as fh:
        manifest.extend(list(csv.DictReader(fh, delimiter="\t")))
log("manifest rows: %d" % len(manifest))


fmt = collections.Counter()


def parse_alleles(snp_id, snp_pos):
    for sep in (":", "_"):
        p = snp_id.split(sep)
        if len(p) == 4 and p[1].lstrip("-").isdigit():
            if int(p[1]) != int(snp_pos):
                return None
            fmt[sep] += 1
            return p[2].upper(), p[3].upper()
    return None


skipped = []

by_chrom = collections.defaultdict(list)
meta = []
unresolved = empty = 0
seen_locus = set()

for i, r in enumerate(manifest):
    if i % 2000 == 0:
        log("  ...%d/%d" % (i, len(manifest)))
    hits = disk.get(os.path.basename(r["weight_path"]), [])
    if not hits:
        unresolved += 1
        continue
    with open(hits[0]) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        empty += 1
        continue
    loc = rows[0]["locus_id"]
    if loc in seen_locus:
        continue
    seen_locus.add(loc)

    n = 0
    for x in rows:
        al = parse_alleles(x["snp_id"], x["snp_pos"])
        if al is None:
            skipped.append((loc, x["snp_id"], x["snp_pos"]))
            continue
        ref, alt = al
        cc = x["snp_chrom"]
        cc = cc if cc.startswith("chr") else "chr" + cc
        by_chrom[cc].append((loc, cc, int(x["snp_pos"]), ref, alt,
                             alt, "%.10g" % float(x["weight"])))
        n += 1

    a = admit.get(loc)
    sn = shipped_n.get(loc, 0)
    meta.append(dict(
        sv_locus=loc, chrom=r["chrom"], sv_index=r["sv_index"],
        sv_category=r["sv_category"], sv_strategy=r["strategy"], sv_arm=r["arm"],
        cv_r2=r["cv_r2"], n_weight_snps=n,
        admit_status=(a["status"] if a else "not_in_admit_table"),
        tier=tier_of(a),
        already_shipped=("yes" if sn else "no"),
        n_snps_previously_shipped=sn,
        n_snps_added_vs_shipped=(n - sn) if sn else "",
    ))

log("\npredictors written : %d" % len(meta))
log("unresolved paths   : %d" % unresolved)
log("empty weight files : %d" % empty)
log("snp_id formats seen: %s" % dict(fmt))
log("SNP rows skipped   : %d" % len(skipped))
for s in skipped[:10]:
    log("  SKIPPED %s snp_id=%r snp_pos=%s" % s)

zero = [m["sv_locus"] for m in meta if m["n_weight_snps"] == 0]
if zero or skipped:
    log("\nFATAL: %d predictors have zero weight rows, %d SNP rows unparsed"
        % (len(zero), len(skipped)))
    for z in zero[:10]:
        log("  ZERO-WEIGHT %s" % z)
    sys.exit(2)

def ckey(c):
    c = c.replace("chr", "")
    return (0, int(c)) if c.isdigit() else (1, c)


total_rows = 0
for cc in sorted(by_chrom, key=ckey):
    rows = by_chrom[cc]
    total_rows += len(rows)
    parts = ([rows[i:i + MAXROWS] for i in range(0, len(rows), MAXROWS)] or [[]]
             if MAXROWS else [rows])
    for k, chunk in enumerate(parts):
        suffix = "" if len(parts) == 1 else "_part%d" % (k + 1)
        fp = "%s/sv_weights_long_ext2_%s.tsv%s" % (OUT, cc, suffix)
        with open(fp, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sv_locus", "chrom", "pos", "ref", "alt", "a1_effect",
                        "beta_std"])
            w.writerows(chunk)
        log("  %-46s %7d rows  %5.2f MB"
            % (os.path.basename(fp), len(chunk), os.path.getsize(fp) / 1e6))

mcols = ["sv_locus", "chrom", "sv_index", "sv_category", "sv_strategy", "sv_arm",
         "cv_r2", "n_weight_snps", "admit_status", "tier", "already_shipped",
         "n_snps_previously_shipped", "n_snps_added_vs_shipped"]
mp = OUT + "/sv_metadata_ext2.tsv"
with open(mp, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=mcols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(sorted(meta, key=lambda m: (ckey(m["chrom"]), int(m["sv_index"]))))
log("\nwrote %s (%d rows, %.2f MB)" % (os.path.basename(mp), len(meta),
                                       os.path.getsize(mp) / 1e6))

log("\ntotal weight rows: %d" % total_rows)
t = collections.Counter(m["tier"] for m in meta)
log("\nby tier:")
for k, v in sorted(t.items()):
    log("  %-28s %6d predictors" % (k, v))

new = [m for m in meta if m["already_shipped"] == "no"]
old = [m for m in meta if m["already_shipped"] == "yes"]
log("\nNEW to All of Us      : %6d predictors" % len(new))
log("already shipped       : %6d predictors" % len(old))
if old:
    d = sorted(int(m["n_snps_added_vs_shipped"]) for m in old)
    log("  SNPs added vs what they hold: median %d, p90 %d, max %d, unchanged %d"
        % (d[len(d) // 2], d[int(len(d) * 0.9)], d[-1], sum(1 for x in d if x == 0)))

sys.exit(0 if unresolved == 0 else 1)
