#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path


SUBSTRATE_COLS = [
    "sv_locus", "sv_index", "sv_id", "sv_category", "strategy", "arm",
    "cv_r2", "primary_chrom", "primary_block_start", "primary_block_end",
    "n_weight_snps", "n_matched", "ld_match_frac", "n_in_block", "status",
]
OUTPUT_COLS = [
    "pair_id", "sv_locus", "sv_id_builder", "sv_chrom", "svpos",
    "strategy", "arm", "win_start", "win_end", "bundle_dir",
    "sv_category", "cv_r2",
]


def parse_sv_locus(sv_locus: str, sv_index_expected: str) -> str:
    parts = sv_locus.split("_")
    if len(parts) < 3:
        raise ValueError(f"sv_locus too short to parse: {sv_locus!r}")
    if not parts[1].startswith("sv"):
        raise ValueError(f"sv_locus[1] does not start with 'sv': {sv_locus!r}")
    locus_idx = parts[1][2:]
    if locus_idx != sv_index_expected:
        raise ValueError(
            f"sv_index mismatch sv_locus={sv_locus!r} idx_token={locus_idx!r} "
            f"expected={sv_index_expected!r}"
        )
    svtype = parts[2]
    if svtype not in {"BND", "DEL", "DUP:TANDEM", "INS"}:
        raise ValueError(f"unexpected svtype in sv_locus={sv_locus!r}: {svtype!r}")
    return f"sv_{sv_index_expected}_{svtype}"


def normalize_chrom(primary_chrom: str) -> str:
    c = primary_chrom.strip()
    if c.startswith("chr"):
        return c
    return f"chr{c}"


def atomic_write_tsv(rows, out_path: Path) -> None:
    tmp = out_path.with_suffix(out_path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        fh.write("\t".join(OUTPUT_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in OUTPUT_COLS) + "\n")
    os.replace(tmp, out_path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--substrate", type=Path, required=True,
                    help="admit_decompose_cv01.tsv")
    ap.add_argument("--out-manifest", type=Path, required=True,
                    help="output MASTER manifest TSV")
    args = ap.parse_args()

    rows = []
    pair_id = 0
    with open(args.substrate) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header != SUBSTRATE_COLS:
            raise ValueError(
                f"substrate header mismatch\n"
                f"  expected: {SUBSTRATE_COLS}\n"
                f"  got     : {header}"
            )
        for line_no, line in enumerate(fh, start=2):
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(SUBSTRATE_COLS):
                raise ValueError(
                    f"row width mismatch line {line_no}: "
                    f"{len(parts)} vs expected {len(SUBSTRATE_COLS)}"
                )
            rec = dict(zip(SUBSTRATE_COLS, parts))
            if rec["status"] != "ADMIT":
                continue
            sv_chrom = normalize_chrom(rec["primary_chrom"])
            sv_id_builder = parse_sv_locus(rec["sv_locus"], rec["sv_index"])
            win_start = int(rec["primary_block_start"])
            win_end = int(rec["primary_block_end"])
            svpos = (win_start + win_end) // 2
            arm = rec["arm"]
            if arm not in {"best8", "best6_alt"}:
                raise ValueError(f"unexpected arm={arm!r} at line {line_no}")
            bundle_dir = f"{sv_chrom}/{rec['sv_locus']}__{arm}"
            pair_id += 1
            rows.append({
                "pair_id": pair_id,
                "sv_locus": rec["sv_locus"],
                "sv_id_builder": sv_id_builder,
                "sv_chrom": sv_chrom,
                "svpos": svpos,
                "strategy": rec["strategy"],
                "arm": arm,
                "win_start": win_start,
                "win_end": win_end,
                "bundle_dir": bundle_dir,
                "sv_category": rec["sv_category"],
                "cv_r2": rec["cv_r2"],
            })

    if pair_id != 6519:
        raise ValueError(
            f"expected 6,519 ADMIT rows, got pair_id={pair_id}; "
            "substrate may be stale or wrong file"
        )

    args.out_manifest.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_tsv(rows, args.out_manifest)
    sys.stderr.write(
        f"[manifest] wrote {pair_id} rows to {args.out_manifest}\n"
    )
    best8 = sum(1 for r in rows if r["arm"] == "best8")
    best6 = sum(1 for r in rows if r["arm"] == "best6_alt")
    chroms = {}
    for r in rows:
        chroms[r["sv_chrom"]] = chroms.get(r["sv_chrom"], 0) + 1
    sys.stderr.write(
        f"[breakdown] best8={best8}  best6_alt={best6}  chroms={len(chroms)}\n"
    )
    sys.stderr.write(
        f"[chrom-top5] " +
        "  ".join(f"{c}={n}" for c, n in
                 sorted(chroms.items(), key=lambda kv: -kv[1])[:5]) + "\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
