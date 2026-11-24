#!/usr/bin/env python3
"""Aggregate completed multi-chromosome complete-truth scoring outputs."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

SCOPES = ("raw_full", "production_active", "target_active")


def _hash(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _load(path):
    with path.open() as fh:
        return json.load(fh)


def _table(path, rows, columns=()):
    fields = list(columns)
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, sort_keys=True, separators=(",", ":"))
                             if isinstance(v, (dict, list, tuple)) else v for k, v in row.items()})


def _group_rates(rows, keys, numerator, total, unresolved):
    grouped = defaultdict(lambda: [0, 0, 0])
    for row in rows:
        value = grouped[tuple(row.get(k) for k in keys)]
        value[0] += numerator(row)
        value[1] += total(row)
        value[2] += unresolved(row)
    return [dict(zip(keys, key), numerator=n, denominator=d, unresolved=u,
                 lower=n / d if d else None, upper=(n + u) / d if d else None)
            for key, (n, d, u) in sorted(grouped.items(), key=lambda x: tuple(str(z) for z in x[0]))]


def _size_bin(size):
    return "50-1000" if size <= 1000 else "1001-10000" if size <= 10000 else "10001-50000" if size <= 50000 else ">50000"


def summarize(manifest_path, *, allow_partial=False):
    manifest_path = Path(manifest_path)
    m = _load(manifest_path)
    root = manifest_path.parent
    bundle_ids = [b["id"] for b in m["bundles"]]
    if len(bundle_ids) != len(set(bundle_ids)):
        raise ValueError("Duplicate bundle IDs")
    if any(set(b["simulations"]) != {"control", "alpha_mixed", "TE_mixed"} or
           len(b["simulations"]) != 3 for b in m["bundles"]):
        raise ValueError("Every bundle must specify the three expected simulations")
    simulations, events, calls, errors, windows, paired = [], [], [], [], [], []
    missing, hashes, loaded = [], {str(manifest_path): _hash(manifest_path)}, {}
    strict_versions, classifier_hashes = set(), set()
    classifier_scripts, diagnostic_versions, diagnostic_settings = set(), set(), set()
    for bundle in m["bundles"]:
        for scenario in bundle["simulations"]:
            folder = root / "bundles" / bundle["id"] / scenario
            files = [folder / "truth.json", folder / "scoring/call_errors.json"]
            files += [folder / f"scoring/{scope}.json" for scope in SCOPES]
            files += [folder / f"scoring/parsed_calls_{scope}.json" for scope in SCOPES]
            absent = [str(p) for p in files if not p.is_file()]
            if absent:
                missing.append(dict(bundle=bundle["id"], scenario=scenario, missing=absent))
                continue
            for f in files:
                hashes[str(f)] = _hash(f)
            case = _load(files[0])
            score = {scope: _load(folder / f"scoring/{scope}.json") for scope in SCOPES}
            parsed = {scope: _load(folder / f"scoring/parsed_calls_{scope}.json") for scope in SCOPES}
            diag = _load(files[1]); production = score["production_active"]
            versions = {score[scope].get("conventions", {}).get("exact_search_version") for scope in SCOPES}
            strict_versions.update(versions)
            provenance = diag.get("provenance", {})
            classifier_hashes.add(provenance.get("classifier_sha256"))
            classifier_scripts.add(provenance.get('classifier_script', 'classify_sv_call_errors.py'))
            diagnostic_versions.add(provenance.get('algorithm_version', 'original_bounded_diagnostics'))
            diagnostic_settings.add(json.dumps({key: diag['summary'].get(key) for key in
                ('candidate_flank_bp', 'max_span_bp', 'max_calls', 'distance_cap',
                 'max_comparisons_per_call', 'max_alignment_work', 'max_enumeration_states_per_call')}, sort_keys=True))
            diag_version = provenance.get("strictscorer_version")
            if diag_version is not None and diag_version != production.get("conventions", {}).get("exact_search_version"):
                raise ValueError(f"Diagnostic scorer version mismatch in {folder}")
            meta_hash = production.get("metadata", {}).get("classifier_sha256")
            if meta_hash is not None and meta_hash != provenance.get("classifier_sha256"):
                raise ValueError(f"Diagnostic classifier hash mismatch in {folder}")
            if len(diag["rows"]) != len(production["calls"]) or len(parsed["production_active"]) != len(production["calls"]):
                raise ValueError(f"Call row count mismatch in {folder}")
            for i, (d, c, p) in enumerate(zip(diag["rows"], production["calls"], parsed["production_active"])):
                if (d["call_index"] != i or c["index"] != i or d["call_id"] != c["id"] or c["id"] != p["id"] or d["strict_status"] != c["status"] or
                    (d["category"] == "strict_match") != (c["status"] == "tp")):
                    raise ValueError(f"Call alignment mismatch in {folder} at index {i}")
            if len(case["truth"]) != len(production["events"]):
                raise ValueError(f"Truth event count mismatch in {folder}")
            meta = dict(bundle=bundle["id"], case_id=case["id"], chrom=bundle["chrom"],
                        array_rate=bundle["array_rate"], scenario=scenario)
            loaded[(bundle["id"], scenario)] = (case, score)
            for scope in SCOPES:
                s = score[scope]
                if len(s["calls"]) != len(parsed[scope]) or s["metrics"]["n_calls"] != len(s["calls"]):
                    raise ValueError(f"Strict scope mismatch in {folder}: {scope}")
                sims = dict(meta, scope=scope, **s["metrics"])
                sims.update(target_call_count=s["metadata"]["target_call_count"], off_target_call_count=s["metadata"]["off_target_call_count"])
                simulations.append(sims)
            metric = production["metrics"]
            observed = Counter(c["status"] for c in production["calls"])
            if any(observed[status] != metric[key] for status, key in (("tp", "tp_calls"), ("fp", "fp_calls"), ("unresolved", "unresolved_calls"))):
                raise ValueError(f"Strict call accounting mismatch in {folder}")
            if sum(diag["summary"]["category_counts"].values()) != len(diag["rows"]) or Counter(d["category"] for d in diag["rows"]) != Counter(diag["summary"]["category_counts"]):
                raise ValueError(f"Diagnostic category accounting mismatch in {folder}")
            truth_by_id = {t["id"]: t for t in case["truth"]}
            if len(truth_by_id) != len(case["truth"]):
                raise ValueError(f"Duplicate truth IDs in {folder}")
            for ev in production["events"]:
                t = truth_by_id[ev["id"]]
                size = max(t["end"] - t["start"], len(t["alt"]))
                events.append(dict(meta, event_id=ev["id"], true_type=t["kind"], size_bp=size,
                    size_bin=_size_bin(size), pair_id=t.get("pair_id"), sequence_class=t.get("sequence_class"),
                    start=t["start"], end=t["end"], status=ev["status"], recovered=ev["recovered"],
                    matched_call_ids=ev.get("matched_call_ids", [])))
            for c, p, d in zip(production["calls"], parsed["production_active"], diag["rows"]):
                calls.append(dict(meta, call_index=c["index"], call_id=c["id"], call_chrom=c["chrom"],
                                  called_type=c["kind"], status=c["status"], reason=c.get("reason"),
                                  start=c["start"], end=c["end"], edit=p.get("edit"),
                                  matched_truth_ids=c.get("matched_truth_ids", [])))
                errors.append(dict(meta, **d))
            for w in score["target_active"]["negative_windows"]:
                windows.append(dict(meta, **w))
    version_inconsistent = (len(strict_versions) > 1 or len(classifier_hashes) > 1 or
                            len(classifier_scripts) > 1 or len(diagnostic_versions) > 1)
    if version_inconsistent and not allow_partial:
        raise RuntimeError("Mixed strict scorer versions or classifier hashes; rescore before final aggregation")
    if missing and not allow_partial:
        raise RuntimeError(f"Incomplete scoring: {len(missing)} simulations missing; use --allow-partial for inspection")
    primary_events = [r for r in events if r["scenario"] == "alpha_mixed"]
    primary_calls = [r for r in calls if r["scenario"] in ("alpha_mixed", "control")]
    primary_errors = [r for r in errors if r["scenario"] in ("alpha_mixed", "control")]
    sensitivity = lambda rows, keys: _group_rates(rows, keys, lambda r: int(r["recovered"]),
                                                  lambda r: 1, lambda r: int(r["status"] == "unresolved"))
    sens_type = sensitivity(primary_events, ["true_type"])
    sens_chr = sensitivity(primary_events, ["chrom", "true_type", "array_rate"])
    sens_size = sensitivity(primary_events, ["true_type", "size_bin", "array_rate"])
    precision_rows = []
    for condition, scenarios in (("primary", ("alpha_mixed", "control")), ("TE_comparison", ("TE_mixed", "control"))):
        subset = [dict(r, condition=condition) for r in calls if r["scenario"] in scenarios]
        precision_rows.extend(_group_rates(subset, ["condition", "called_type"],
            lambda r: int(r["status"] == "tp"), lambda r: 1, lambda r: int(r["status"] == "unresolved")))
    for condition in ("primary", "TE_comparison"):
        for typ in ("INS", "DEL", "DUP:TANDEM", "DUP:INT", "INV", "BND"):
            if not any(r["condition"] == condition and r["called_type"] == typ for r in precision_rows):
                precision_rows.append(dict(condition=condition, called_type=typ, numerator=0, denominator=0,
                                           unresolved=0, lower=None, upper=None))
    control_rows = []
    controls = [r for r in simulations if r["scope"] == "production_active" and r["scenario"] == "control"]
    types = sorted({r["called_type"] for r in calls} | {"INS", "DEL", "DUP:TANDEM", "DUP:INT", "INV", "BND"})
    for control in controls:
        key = control["bundle"]
        selected = [r for r in calls if r["bundle"] == key and r["scenario"] == "control"]
        full_windows = [r for r in windows if r["bundle"] == key and r["scenario"] == "control" and r["end"]-r["start"] == 50000 and r["eligible_negative"]]
        for typ in ["ALL"] + types:
            global_calls = [r for r in selected if typ == "ALL" or r["called_type"] == typ]
            target_calls = [r for r in global_calls if r["call_chrom"] == control["chrom"]]
            ids = {r["call_id"] for r in target_calls if r["status"] == "fp"}
            unresolved_ids = {r["call_id"] for r in target_calls if r["status"] == "unresolved"}
            affected = sum(bool(ids.intersection(w["fp_call_ids"])) for w in full_windows)
            uncertain = sum(not ids.intersection(w["fp_call_ids"]) and bool(unresolved_ids.intersection(w["unresolved_call_ids"])) for w in full_windows)
            control_rows.append(dict(bundle=key, chrom=control["chrom"], array_rate=control["array_rate"], called_type=typ,
                global_records=len(global_calls), target_records=len(target_calls), off_target_records=len(global_calls)-len(target_calls),
                global_empty=int(not global_calls), target_empty=int(not target_calls),
                full_windows=len(full_windows), affected_windows=affected, unresolved_windows=uncertain,
                clean_windows=len(full_windows)-affected-uncertain))
    for bundle in m["bundles"]:
        left = loaded.get((bundle["id"], "alpha_mixed")); right = loaded.get((bundle["id"], "TE_mixed"))
        if not left or not right:
            continue
        if left[0].get("template_sha256") != right[0].get("template_sha256"):
            raise ValueError(f"Paired payload background mismatch: {bundle['id']}")
        lt = {x["id"]: x for x in left[0]["truth"]}; rt = {x["id"]: x for x in right[0]["truth"]}
        le = {x["id"]: x for x in left[1]["production_active"]["events"]}; re = {x["id"]: x for x in right[1]["production_active"]["events"]}
        for name in ("paired_299", "paired_6087"):
            if name not in lt or name not in rt:
                continue
            a, b = lt[name], rt[name]
            if (a["chrom"], a["start"], a["end"], len(a["alt"]), a.get("pair_id")) != (b["chrom"], b["start"], b["end"], len(b["alt"]), b.get("pair_id")):
                raise ValueError(f"Paired payload geometry mismatch: {bundle['id']} {name}")
            paired.append(dict(bundle=bundle["id"], chrom=bundle["chrom"], array_rate=bundle["array_rate"],
                event_id=name, pair_id=a.get("pair_id"), size_bp=len(a["alt"]), start=a["start"],
                alpha_status=le[name]["status"], alpha_recovered=le[name]["recovered"],
                TE_status=re[name]["status"], TE_recovered=re[name]["recovered"]))
    categories = Counter(r["category"] for r in primary_errors)
    primary_counts = Counter(r["status"] for r in primary_calls)
    exact = [r["dna_edit_distance"] for r in primary_errors if r.get("dna_edit_distance") is not None]
    absolute = {label: sum((lo is None or x > lo) and (hi is None or x <= hi) for x in exact)
                for label, lo, hi in (("0", None, 0), ("1", 0, 1), ("2-10", 1, 10), ("11-100", 10, 100), ("101-1000", 100, 1000), (">1000", 1000, None))}
    relative = [r["dna_error_fraction_of_event"] for r in primary_errors if r.get("dna_error_fraction_of_event") is not None]
    relative_bins = {label: sum((lo is None or x > lo) and (hi is None or x <= hi) for x in relative)
                     for label, lo, hi in (("0", None, 0), ("(0,.001]", 0, .001), ("(.001,.01]", .001, .01), ("(.01,.05]", .01, .05), ("(.05,.5]", .05, .5), (">.5", .5, None))}
    expected = sum(len(b["simulations"]) for b in m["bundles"])
    if "array_substitution_rates" in m["config"]:
        design_expected = len(m["config"]["chromosomes"]) * len(m["config"]["array_substitution_rates"]) * 3
        if expected != design_expected:
            raise ValueError(f"Manifest has {expected} simulations, expected {design_expected} from configuration")
    summary = dict(complete=not missing and not version_inconsistent, partial=bool(missing) or version_inconsistent,
        version_inconsistent=version_inconsistent,
        strict_scorer_versions=sorted(strict_versions, key=str), classifier_hashes=sorted(classifier_hashes, key=str),
        classifier_script=next(iter(classifier_scripts)) if len(classifier_scripts) == 1 else None,
        diagnostic_algorithm_versions=sorted(diagnostic_versions),
        diagnostic_search_settings=[json.loads(x) for x in sorted(diagnostic_settings)],
        diagnostic_search_incomplete_records=sum(r['category'] == 'unsupported_or_bounded' for r in primary_errors),
        diagnostic_bounded_records=sum(r.get('search_bounded') is True for r in primary_errors),
        diagnostic_resolution_counts=dict(Counter(r.get('explanation_resolution', 'not_required') for r in primary_errors)),
        expected_simulations=expected,
        completed_simulations=len(loaded), missing_simulations=missing,
        chromosomes_expected=m["config"]["chromosomes"], chromosomes_completed=sorted({k[0].split('_')[0] for k in loaded}),
        primary_scope="production_active", primary_conditions=["alpha_mixed", "control"],
        primary_strict_counts=dict(primary_counts), primary_records=len(primary_calls),
        primary_diagnostic_categories=dict(categories), exact_dna_error_absolute_bands=absolute,
        exact_dna_error_relative_bands=relative_bins,
        dna_error_censored_count=sum(r.get("dna_edit_distance_censored") is True for r in primary_errors),
        diagnostic_rows_with_lower_bound=sum(r.get("dna_edit_distance_lower_bound") is not None for r in primary_errors),
        censored_lower_bound_bands={label: sum(r.get("dna_edit_distance_censored") is True and
            r.get("dna_edit_distance_lower_bound") is not None and
            (lo is None or r["dna_edit_distance_lower_bound"] > lo) and
            (hi is None or r["dna_edit_distance_lower_bound"] <= hi) for r in primary_errors)
            for label, lo, hi in (("0", None, 0), ("1", 0, 1), ("2-10", 1, 10),
                                  ("11-100", 10, 100), ("101-1000", 100, 1000), (">1000", 1000, None))},
        group_records=sum(len(r.get("group_call_ids", [])) > 1 for r in primary_errors),
        improvement_proven_records=sum(r.get("improvement_proven") is True for r in primary_errors),
        control_count=len(controls), paired_payload_count=len(paired), source_sha256=hashes,
        taxonomy="Strict match, strict duplicate, near-correct, distorted, off destination, near unattributed, unsupported/bounded, or null control false call; per-record local groups can overlap and cannot estimate relaxed precision.",
        definitions="Primary alpha_mixed plus controls; TE_mixed is separate. Numerator/denominator bounds preserve unresolved; no independence-based intervals.")
    error_categories_by_rate = [dict(array_rate=k[0], category=k[1], count=v)
                                for k, v in sorted(Counter((r["array_rate"], r["category"]) for r in primary_errors).items())]
    error_magnitude = [dict(bundle=r["bundle"], case_id=r["case_id"], chrom=r["chrom"],
                            array_rate=r["array_rate"], scenario=r["scenario"], call_index=r["call_index"],
                            call_id=r["call_id"], category=r["category"],
                            dna_edit_distance=r.get("dna_edit_distance"),
                            dna_edit_distance_lower_bound=r.get("dna_edit_distance_lower_bound"),
                            dna_edit_distance_censored=r.get("dna_edit_distance_censored"),
                            dna_error_fraction_of_event=r.get("dna_error_fraction_of_event"),
                            improvement_proven=r.get("improvement_proven"),
                            marginal_improvement_proven=r.get("marginal_improvement_proven"),
                            group_call_ids=r.get("group_call_ids", [])) for r in primary_errors]
    out = root / "summary"; out.mkdir(parents=True, exist_ok=True)
    for name, rows in [("per_simulation", simulations), ("events", events), ("calls", calls),
        ("call_errors_flat", errors), ("sensitivity_by_type", sens_type), ("by_chr_type_rate", sens_chr),
        ("by_size_type_rate", sens_size), ("precision_by_reported_type", precision_rows),
        ("control_specificity", control_rows), ("paired_payload_recovery", paired),
        ("error_categories_by_rate", error_categories_by_rate), ("error_magnitude", error_magnitude)]:
        _table(out / f"{name}.tsv", rows, columns={
            "per_simulation": ("bundle", "case_id", "chrom", "array_rate", "scenario", "scope"),
            "events": ("bundle", "case_id", "chrom", "array_rate", "scenario", "event_id"),
            "calls": ("bundle", "case_id", "chrom", "array_rate", "scenario", "call_index", "call_id"),
            "call_errors_flat": ("bundle", "case_id", "chrom", "array_rate", "scenario", "call_index", "call_id", "category"),
            "sensitivity_by_type": ("true_type", "numerator", "denominator", "unresolved", "lower", "upper"),
            "by_chr_type_rate": ("chrom", "true_type", "array_rate", "numerator", "denominator", "unresolved", "lower", "upper"),
            "by_size_type_rate": ("true_type", "size_bin", "array_rate", "numerator", "denominator", "unresolved", "lower", "upper"),
            "precision_by_reported_type": ("condition", "called_type", "numerator", "denominator", "unresolved", "lower", "upper"),
            "control_specificity": ("bundle", "chrom", "array_rate", "called_type", "global_empty", "target_empty", "full_windows", "clean_windows", "affected_windows", "unresolved_windows"),
            "paired_payload_recovery": ("bundle", "chrom", "array_rate", "event_id", "pair_id", "alpha_status", "alpha_recovered", "TE_status", "TE_recovered"),
            "error_categories_by_rate": ("array_rate", "category", "count"),
            "error_magnitude": ("bundle", "case_id", "chrom", "array_rate", "scenario", "call_index", "call_id", "category", "dna_edit_distance", "dna_edit_distance_lower_bound", "dna_edit_distance_censored")
        }.get(name, ()))
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()
    result = summarize(args.manifest, allow_partial=args.allow_partial)
    print(f"Completed {result['completed_simulations']}/{result['expected_simulations']} simulations; partial={result['partial']}")
