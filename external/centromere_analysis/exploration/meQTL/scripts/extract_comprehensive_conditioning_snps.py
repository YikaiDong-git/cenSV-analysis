#!/usr/bin/env python3
"""Extract one shared phenotype-blind SNP matrix for an SV/window/callset task."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from comprehensive_conditioning import interval_covered, ld_prune, phenotype_independent_snp_selection


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--matrix-output", type=Path, required=True)
    parser.add_argument("--variants-output", type=Path, required=True)
    parser.add_argument("--status-output", type=Path, required=True)
    return parser.parse_args()


def write_outputs(args: argparse.Namespace, status: dict, samples: list[str], metadata: pd.DataFrame | None = None,
                  genotypes: np.ndarray | None = None, sv: np.ndarray | None = None,
                  condition_indices: np.ndarray | None = None) -> None:
    args.matrix_output.parent.mkdir(parents=True, exist_ok=True)
    metadata = metadata if metadata is not None else pd.DataFrame(columns=["variant_id", "chrom", "pos", "ref", "alt", "maf", "call_rate", "sv_abs_correlation", "condition_selected"])
    genotypes = genotypes if genotypes is not None else np.empty((len(samples), 0), dtype=np.float32)
    sv = sv if sv is not None else np.full(len(samples), np.nan)
    condition_indices = condition_indices if condition_indices is not None else np.array([], dtype=np.int32)
    compression = {"method": "gzip", "mtime": 0} if args.variants_output.suffix == ".gz" else None
    metadata.to_csv(args.variants_output, sep="\t", index=False, compression=compression)
    np.savez_compressed(args.matrix_output, sample_ids=np.asarray(samples), genotypes=np.asarray(genotypes, dtype=np.float32),
                        sv_dosage=np.asarray(sv, dtype=np.float64), variant_ids=metadata.get("variant_id", pd.Series(dtype=str)).astype(str).to_numpy(),
                        condition_indices=np.asarray(condition_indices, dtype=np.int32))
    args.status_output.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")


def parse_call(value: str, min_dp: float, min_gq: float) -> float:
    fields = (value.split(":") + [".", "."])[:3]
    gt, dp, gq = fields
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or "." in alleles:
        return np.nan
    if min_dp > 0 and (dp == "." or float(dp) < min_dp):
        return np.nan
    if min_gq > 0 and (gq == "." or float(gq) < min_gq):
        return np.nan
    return float(int(alleles[0]) + int(alleles[1]))


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text())
    root = Path(config["workflow_output"])
    tasks = pd.read_csv(root / "manifests/snp_tasks.tsv", sep="\t").set_index("snp_task_id")
    svs = pd.read_csv(root / "manifests/svs.tsv", sep="\t").set_index("sv_id")
    task = tasks.loc[args.task_id]
    svrow = svs.loc[str(task.sv_id)]
    panel_name = str(task.panel)
    panel = config["panels"][panel_name]
    design = np.load(config["canonical_design"], allow_pickle=True)
    samples = design["sample_ids"].astype(str).tolist()
    if "pc_counts" in design.files:
        cov_counts = dict(zip(design["pc_counts"].astype(int), design["covariate_counts"].astype(int)))
        covariates = np.asarray(design["covariates"][:, :cov_counts[10]], float)
    else:
        covariates = np.asarray(design["covariates"], float)
    start = max(0, int(svrow.pos) - int(task.window_bp))
    end = int(svrow.end) + int(task.window_bp)
    chrom = str(svrow.chrom)
    vcf = Path(str(panel["vcf_template"]).format(chrom=chrom))
    coverage_value = panel.get("coverage_template")
    coverage = Path(str(coverage_value).format(chrom=chrom)) if coverage_value else None
    base = {"task_id": args.task_id, "sv_id": str(task.sv_id), "panel": panel_name,
            "panel_interpretation": "same-participant sensitivity callset; not independent replication",
            "window_bp": int(task.window_bp), "chrom": chrom, "start": start, "end": end,
            "vcf": str(vcf), "coverage": str(coverage) if coverage else None,
            "selection_uses_methylation": False}
    available_windows = {int(value) for value in panel.get("available_windows_bp", [])}
    if int(task.window_bp) not in available_windows:
        write_outputs(args, {**base, "status": "callset_window_unavailable",
                             "unavailable_reason": str(panel.get("unavailable_reason", "not configured as available"))}, samples)
        return
    if not vcf.exists() or (coverage is not None and not coverage.exists()):
        write_outputs(args, {**base, "status": "missing_panel_input"}, samples)
        return
    if not Path(str(vcf) + ".tbi").exists():
        write_outputs(args, {**base, "status": "missing_panel_index"}, samples)
        return
    if coverage is not None and not interval_covered(coverage, chrom, start, end):
        write_outputs(args, {**base, "status": "incomplete_interval_coverage"}, samples)
        return

    dosage_table = pd.read_csv(svrow.sv_dosages, sep="\t")
    dosage_hit = dosage_table.loc[pd.to_numeric(dosage_table.sv_index).eq(int(svrow.sv_index))]
    if len(dosage_hit) != 1:
        write_outputs(args, {**base, "status": "sv_dosage_lookup_failed"}, samples)
        return
    sv = pd.to_numeric(dosage_hit.iloc[0].reindex(samples), errors="coerce").to_numpy(float)
    if not np.isfinite(sv).all():
        write_outputs(args, {**base, "status": "sv_dosage_missing"}, samples, sv=sv)
        return

    max_raw = int(config["max_raw_snps"])
    metadata: list[dict[str, object]] = []
    genotype_rows: list[np.ndarray] = []
    tmp_root = os.environ.get("TMPDIR")
    with tempfile.TemporaryDirectory(prefix="meqtl_snp_", dir=tmp_root) as temp_dir:
        sample_file = Path(temp_dir) / "samples.txt"
        sample_file.write_text("\n".join(samples) + "\n")
        command = ["bcftools", "query", "-S", str(sample_file), "-r", f"{chrom}:{start + 1}-{end}",
                   "-i", 'TYPE="snp" && N_ALT=1', "-f", "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT:%DP:%GQ]\\n", str(vcf)]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        assert process.stdout is not None
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4 + len(samples):
                process.kill()
                raise RuntimeError(f"Unexpected bcftools field count for {args.task_id}")
            calls = np.asarray([parse_call(value, float(config["min_dp"]), float(config["min_gq"]))
                                for value in fields[4:]], float)
            call_rate = float(np.isfinite(calls).mean())
            if call_rate < float(config["min_site_call_rate"]):
                continue
            calls[np.isnan(calls)] = np.nanmean(calls)
            allele_frequency = float(calls.mean() / 2)
            maf = min(allele_frequency, 1 - allele_frequency)
            if maf < float(config["maf_min"]) or calls.std() == 0:
                continue
            metadata.append({"variant_id": f"{fields[0]}:{fields[1]}:{fields[2]}:{fields[3]}", "chrom": fields[0],
                             "pos": int(fields[1]), "ref": fields[2], "alt": fields[3], "maf": maf, "call_rate": call_rate})
            genotype_rows.append(calls)
            if len(genotype_rows) > max_raw:
                process.kill()
                write_outputs(args, {**base, "status": "oversized_raw_snp_count", "max_raw_snps": max_raw}, samples, sv=sv)
                return
        stderr = process.stderr.read() if process.stderr else ""
        if process.wait() != 0:
            raise RuntimeError(f"bcftools query failed: {stderr[-2000:]}")
    if not genotype_rows:
        write_outputs(args, {**base, "status": "no_eligible_snps", "n_raw_eligible": 0}, samples, sv=sv)
        return
    raw = np.vstack(genotype_rows).T
    ids = [str(row["variant_id"]) for row in metadata]
    selected, scores = phenotype_independent_snp_selection(raw, sv, covariates, int(config["max_susie_snps"]), ids)
    selected_raw = raw[:, selected]
    selected_meta = pd.DataFrame([metadata[index] for index in selected]).reset_index(drop=True)
    selected_meta["sv_abs_correlation"] = scores[selected]
    condition_relative = ld_prune(np.arange(len(selected)), selected_raw, int(config["max_condition_snps"]), float(config["condition_ld_r2_max"]))
    selected_meta["condition_selected"] = False
    selected_meta.loc[condition_relative, "condition_selected"] = True
    status = {**base, "status": "ok", "n_raw_eligible": len(raw.T), "n_susie_snps": len(selected),
              "n_condition_snps": len(condition_relative), "phenotype_independent_selection": "absolute partial correlation with SV after PC10 covariates; stable variant-ID tie break; LD prune for conditional subset"}
    write_outputs(args, status, samples, selected_meta, selected_raw, sv, condition_relative)


if __name__ == "__main__":
    main()
