#!/usr/bin/env python3
"""Prepare nonredundant observed-SV dosage profiles for the focused trans scan."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROBUST_STATUSES = {"cross_panel_robust", "single_panel_robust"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consensus", type=Path, required=True)
    parser.add_argument("--candidate-evidence", type=Path, required=True)
    parser.add_argument("--output-dosages", type=Path, required=True)
    parser.add_argument("--output-metadata", type=Path, required=True)
    parser.add_argument("--output-candidates", type=Path, required=True)
    return parser.parse_args()


def dosage_key(values: np.ndarray) -> bytes:
    normalized = np.nan_to_num(values.astype(np.float64), nan=-9.0)
    return normalized.tobytes()


def main() -> int:
    args = parse_args()
    consensus = pd.read_csv(args.consensus, sep="\t")
    evidence = pd.read_csv(args.candidate_evidence, sep="\t")
    selected = consensus[consensus["conditioning_consensus"].isin(ROBUST_STATUSES)].copy()
    if selected.empty:
        raise ValueError("No robust candidates found in conditioning consensus")

    evidence_cols = [
        "sv_id",
        "sv_index",
        "sv_dosages",
        "sv_chrom",
        "sv_pos",
        "sv_end",
        "svtype",
        "svlen",
        "cv_r2",
        "svfreq",
        "n_observed_carriers",
    ]
    evidence_subset = evidence[evidence_cols].drop_duplicates("sv_id").rename(
        columns={column: f"evidence_{column}" for column in evidence_cols if column != "sv_id"}
    )
    selected = selected.merge(
        evidence_subset, on="sv_id", how="left", validate="one_to_one"
    )
    if selected[["evidence_sv_index", "evidence_sv_dosages"]].isna().any().any():
        missing = selected.loc[selected["evidence_sv_dosages"].isna(), "sv_id"].tolist()
        raise ValueError(f"Missing dosage provenance for candidates: {missing}")

    dosage_rows: list[dict] = []
    metadata_rows: list[dict] = []
    sample_order: list[str] | None = None
    seen_profiles: dict[bytes, int] = {}

    priority = {"cross_panel_robust": 0, "single_panel_robust": 1}
    selected["_priority"] = selected["conditioning_consensus"].map(priority)
    selected = selected.sort_values(
        ["_priority", "physical_locus_id", "evidence_cv_r2"], ascending=[True, True, False]
    )

    for row in selected.itertuples(index=False):
        dosage = pd.read_csv(row.evidence_sv_dosages, sep="\t")
        source_index = int(row.evidence_sv_index)
        match = dosage[dosage["sv_index"] == source_index]
        if len(match) != 1:
            raise ValueError(f"Expected one dosage row for {row.sv_id}; found {len(match)}")
        samples = [column for column in dosage.columns if column != "sv_index"]
        if sample_order is None:
            sample_order = samples
        elif samples != sample_order:
            raise ValueError(f"Sample order differs in {row.evidence_sv_dosages}")

        values = match.iloc[0][samples].to_numpy(dtype=float)
        key = dosage_key(values)
        if key in seen_profiles:
            output_index = seen_profiles[key]
            meta = metadata_rows[output_index]
            meta["collapsed_sv_ids"] = f"{meta['collapsed_sv_ids']};{row.sv_id}"
            meta["collapsed_source_indices"] = f"{meta['collapsed_source_indices']};{source_index}"
            meta["n_collapsed_records"] += 1
            continue

        output_index = len(dosage_rows)
        seen_profiles[key] = output_index
        dosage_rows.append({"sv_index": output_index, **dict(zip(samples, values, strict=True))})
        metadata_rows.append(
            {
                "sv_index": output_index,
                "id": row.sv_id,
                "chrom": row.sv_chrom,
                "pos": int(row.sv_pos),
                "end": int(row.evidence_sv_end),
                "svtype": row.evidence_svtype,
                "svlen": row.evidence_svlen,
                "cv_r2": row.evidence_cv_r2,
                "svfreq": row.evidence_svfreq,
                "physical_locus_id": row.physical_locus_id,
                "conditioning_consensus": row.conditioning_consensus,
                "n_robust_panels": int(row.n_robust_panels),
                "n_observed_carriers": row.evidence_n_observed_carriers,
                "source_sv_index": source_index,
                "source_dosages": row.evidence_sv_dosages,
                "collapsed_sv_ids": row.sv_id,
                "collapsed_source_indices": str(source_index),
                "n_collapsed_records": 1,
            }
        )

    for path in [args.output_dosages, args.output_metadata, args.output_candidates]:
        path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dosage_rows).to_csv(args.output_dosages, sep="\t", index=False)
    metadata = pd.DataFrame(metadata_rows)
    metadata.to_csv(args.output_metadata, sep="\t", index=False)
    metadata.to_csv(args.output_candidates, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
