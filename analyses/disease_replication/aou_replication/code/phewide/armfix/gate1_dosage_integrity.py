#!/usr/bin/env python3
import os, glob, hashlib
import numpy as np
import pgenlib

HOME = os.path.expanduser("~")
RV = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds"
NEW = f"{RV}/step28_phewide/pgen"
OLD = f"{RV}/archive_2026-07-11_pre_armfix/step28_phewide"
AFFECTED = {"chr1", "chr9", "chr16"}


def sv_hashes(pgdir):
    prefixes = sorted(glob.glob(f"{pgdir}/svpgen300_*.pgen"),
                      key=lambda p: int(p.split("svpgen300_")[1].split(".pgen")[0]))
    out = {}
    for pg in prefixes:
        ids = [ln.split("\t")[2] for ln in open(pg[:-5] + ".pvar") if not ln.startswith("#")]
        r = pgenlib.PgenReader(pg.encode())
        n = r.get_raw_sample_ct()
        buf = np.empty(n, np.float32)
        for v in range(r.get_variant_ct()):
            r.read_dosages(v, buf)
            out[ids[v]] = hashlib.md5(buf.tobytes()).hexdigest()
        r.close()
    return out


def main():
    print("hashing NEW pgen dosages ...", flush=True); new = sv_hashes(NEW)
    print("hashing ARCHIVED pgen dosages ...", flush=True); old = sv_hashes(OLD)
    common = sorted(set(new) & set(old))
    print(f"new SVs={len(new)} old SVs={len(old)} common={len(common)}", flush=True)

    unaff_match = unaff_diff = aff_match = aff_diff = 0
    bad_unaff, unchanged_aff = [], []
    for sv in common:
        chrom = sv.split("_")[0]
        same = new[sv] == old[sv]
        if chrom in AFFECTED:
            if same:
                aff_match += 1; unchanged_aff.append(sv)
            else:
                aff_diff += 1
        else:
            if same:
                unaff_match += 1
            else:
                unaff_diff += 1; bad_unaff.append(sv)

    print("\n=== GATE 1 RESULT ===")
    print(f"  UNAFFECTED (18 chroms): {unaff_match} bit-identical, {unaff_diff} DIFFER (must be 0)")
    print(f"  AFFECTED  (chr1/9/16):  {aff_diff} changed, {aff_match} unchanged (expect ~all changed)")
    if bad_unaff:
        print(f"  !! FAIL: {len(bad_unaff)} unaffected-chrom SVs changed, e.g. {bad_unaff[:8]}")
    if unchanged_aff:
        print(f"  note: {len(unchanged_aff)} affected-chrom SVs unchanged "
              f"(OK if their predictors were all inside the old window), e.g. {unchanged_aff[:5]}")
    ok = (unaff_diff == 0) and (aff_diff > 0)
    print(f"\n  GATE 1: {'PASS ✓' if ok else 'FAIL ✗'}")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
