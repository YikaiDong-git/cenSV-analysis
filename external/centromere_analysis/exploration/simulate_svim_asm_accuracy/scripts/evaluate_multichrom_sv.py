"""Score one multi-chromosome simulated assembly against complete sequence truth."""
import json
import copy
import hashlib
import shutil
import os
from pathlib import Path
import tempfile

from classify_sv_call_errors import classify_call_errors
from complete_truth_reference import CompositeReference
from multichrom_truth_scoring import score_complete_truth
from parse_complete_truth_calls import parse_complete_calls

SCOPES = ("raw_full", "production_active", "target_active")


def _write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write("\n")
        os.replace(name, path)
    except BaseException:
        if os.path.exists(name):
            os.unlink(name)
        raise


def _validate_template(reference, target_chrom):
    if target_chrom not in reference.template.references:
        raise ValueError(f"Mutated template missing target contig {target_chrom}")
    if target_chrom not in reference.original.references:
        raise ValueError(f"Original reference missing target contig {target_chrom}")
    if reference.template.get_reference_length(target_chrom) != reference.original.get_reference_length(target_chrom):
        raise ValueError("Substitution template and reference lengths differ")


def validate_truth(template_path, reference_path, truth, windows):
    reference = CompositeReference(template_path, reference_path)
    try:
        for event in truth:
            _validate_template(reference, event["chrom"])
        for window in windows:
            chrom, start, end = window["chrom"], window["start"], window["end"]
            if not 0 <= start < end <= reference.get_reference_length(chrom):
                raise ValueError(f"Invalid negative window: {window.get('id')}")
        raw_type = {"TANDEM": "DUP:TANDEM", "DISPERSED": "DUP:INT"}
        calls = [dict(id="oracle_" + event["id"], chrom=event["chrom"],
                      kind=raw_type.get(event["kind"], event["kind"]),
                      edit=(event["start"], event["end"], event["alt"])) for event in truth]
        result = score_complete_truth(reference, truth, calls, windows)
        metrics = result["metrics"]
        checks = dict(all_truth_recovered=metrics["tp_events"] == len(truth),
                      no_false_or_unresolved_calls=metrics["fp_calls"] == 0 and metrics["unresolved_calls"] == 0,
                      no_unresolved_events=metrics["unresolved_events"] == 0,
                      whole_sequence_reconstructed=metrics["full_sequence_reconstruction"] is True,
                      no_false_negative_windows=metrics["fp_windows"] == 0 and
                                                metrics["unresolved_negative_windows"] == 0)
        return dict(passed=all(checks.values()), checks=checks, metrics=metrics)
    finally:
        reference.close()


def score_simulation(case_dir, template_path, reference_path, target_chrom, *, diagnostics=True):
    case_dir = Path(case_dir)
    case = json.loads((case_dir / "truth.json").read_text())
    truth, windows = case["truth"], case["negative_windows"]
    if case.get("chrom", target_chrom) != target_chrom:
        raise ValueError("Case chromosome disagrees with target chromosome")
    if any(t["chrom"] != target_chrom for t in truth):
        raise ValueError("Truth contains events on a non-target chromosome")
    if any(w["chrom"] != target_chrom for w in windows):
        raise ValueError("Negative windows contain a non-target chromosome")
    trace = case_dir / "svim" / "candidate_trace.json"
    reference = CompositeReference(template_path, reference_path)
    try:
        _validate_template(reference, target_chrom)
        raw = parse_complete_calls(case_dir / "svim" / "variants.vcf", reference, trace)
        active = parse_complete_calls(case_dir / "active.vcf.gz", reference, trace)
        raw_ids = [c["id"] for c in raw]
        active_ids = [c["id"] for c in active]
        if None in raw_ids or len(set(raw_ids)) != len(raw_ids):
            raise ValueError("Raw VCF requires unique, nonempty call IDs")
        if len(set(active_ids)) != len(active_ids) or not set(active_ids) <= set(raw_ids):
            raise ValueError("Filtered VCF call IDs must be unique and present in raw VCF")
        subsets = {"raw_full": raw, "production_active": active,
                   "target_active": [c for c in active if c["chrom"] == target_chrom]}
        out = case_dir / "scoring"
        results = {}
        archive = case_dir / "scoring_initial"
        archived_old = False
        for scope in SCOPES:
            old_path = out / f"{scope}.json"
            if old_path.exists() and not json.loads(old_path.read_text()).get("conventions", {}).get("exact_search_version"):
                archive.mkdir(parents=True, exist_ok=True)
                for source in (old_path, out / f"parsed_calls_{scope}.json"):
                    if source.exists() and not (archive / source.name).exists():
                        shutil.copy2(source, archive / source.name)
                archived_old = True
        if archived_old and (out / "call_errors.json").exists() and not (archive / "call_errors.json").exists():
            shutil.copy2(out / "call_errors.json", archive / "call_errors.json")
        classifier_sha = hashlib.sha256((Path(__file__).parent / "classify_sv_call_errors.py").read_bytes()).hexdigest()
        cache = []
        for scope, selected in subsets.items():
            cached = next((result for prior, result in cache if prior == selected), None)
            if cached is None:
                cached = score_complete_truth(reference, truth, selected, windows)
                cache.append((selected, cached))
            result = copy.deepcopy(cached)
            result["metadata"] = dict(scope=scope, case_id=case.get("id"), target_chrom=target_chrom,
                                      scenario=case.get("scenario"), array_rate=case.get("array_rate"),
                                      target_call_count=sum(c["chrom"] == target_chrom for c in selected),
                                      off_target_call_count=sum(c["chrom"] != target_chrom for c in selected),
                                      diagnostics_requested=bool(diagnostics),
                                      strictscorer_version=result["conventions"].get("exact_search_version"),
                                      classifier_sha256=classifier_sha if diagnostics else None)
            results[scope] = result
            _write_json(out / f"{scope}.json", result)
            _write_json(out / f"parsed_calls_{scope}.json", selected)
        errors_path = out / "call_errors.json"
        if diagnostics:
            errors = classify_call_errors(reference, truth, active, results["production_active"])
            errors["provenance"] = dict(classifier_sha256=classifier_sha,
                                        strictscorer_version=results["production_active"]["conventions"].get("exact_search_version"))
            _write_json(errors_path, errors)
        elif errors_path.exists():
            errors_path.unlink()
        return dict(case_id=case.get("id"), target_chrom=target_chrom,
                    scoring_dir=str(out), diagnostics=bool(diagnostics),
                    scopes={name: result["metadata"] for name, result in results.items()})
    finally:
        reference.close()
