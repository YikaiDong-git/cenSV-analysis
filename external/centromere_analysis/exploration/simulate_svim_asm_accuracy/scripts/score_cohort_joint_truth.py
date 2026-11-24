"""Exact cohort scoring with separately certified false-junction exclusions."""

from collections import Counter

from cohort_bnd_junction_validation import validate_false_junctions
from multichrom_truth_scoring import _overlap, _span, score_complete_truth


def _ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def score_joint_truth(reference, truth, calls, negative_windows, **strictparams):
    receipts = validate_false_junctions(reference, truth, calls,
        max_representation_shift=strictparams.get("max_representation_shift", 1_000_000))
    excluded = {key for key, value in receipts.items()
                if value["definite_fp"] and value["complete_endpoint_chromosome_absence"]}
    kept = [call for call in calls if call["id"] not in excluded]
    result = score_complete_truth(reference, truth, kept, negative_windows, **strictparams)
    result["junction_validation"] = dict(
        version="exact_bnd_junction_contradiction_v1",
        n_bnd_records=len(receipts), n_proven_false_records=len(excluded),
        n_distinct_adjacencies=len({r["adjacency_id"] for r in receipts.values() if "adjacency_id" in r}),
        n_distinct_proven_false_adjacencies=len({receipts[key]["adjacency_id"] for key in excluded}),
        excluded_from_recovery_search=sorted(excluded), receipts=receipts,
        interpretation="Absent required exact junction in both complete edited endpoint chromosomes disproves a record; "
                       "a present junction never certifies full-allele recovery.")
    result["conventions"]["junction_contradiction_version"] = "exact_bnd_junction_contradiction_v1"
    if not excluded:
        return result

    by_id = {row["id"]: row for row in result["calls"]}
    rows = []
    for index, call in enumerate(calls):
        if call["id"] in excluded:
            start, end = _span(call)
            row = dict(id=call["id"], index=index, chrom=call["chrom"], kind=call["kind"],
                       start=start, end=end, status="fp", reason="proven_false_exact_junction",
                       matched_truth_ids=[], junction_validation=receipts[call["id"]])
        else:
            row = dict(by_id[call["id"]], index=index)
        rows.append(row)
    result["calls"] = rows

    for window in result["negative_windows"]:
        interval = window["start"], window["end"]
        fp, uncertain = [], []
        for call, row in zip(calls, rows):
            if call["chrom"] == window["chrom"] and _overlap(_span(call), interval):
                if row["status"] == "fp":
                    fp.append(row["id"])
                elif row["status"] == "unresolved":
                    uncertain.append(row["id"])
        window.update(fp_call_ids=fp, unresolved_call_ids=uncertain,
                      false_positive_window=bool(fp), unresolved_window=bool(uncertain))
    counts = Counter(row["status"] for row in rows)
    negatives = [w for w in result["negative_windows"] if w["eligible_negative"]]
    fp_windows = sum(w["false_positive_window"] for w in negatives)
    uncertain_only = sum(w["unresolved_window"] and not w["false_positive_window"] for w in negatives)
    result["metrics"].update(
        n_calls=len(calls), tp_calls=counts["tp"], fp_calls=counts["fp"],
        unresolved_calls=counts["unresolved"], precision_lower=_ratio(counts["tp"], len(calls)),
        precision_upper=_ratio(counts["tp"] + counts["unresolved"], len(calls)),
        n_negative_windows=len(negatives), fp_windows=fp_windows,
        unresolved_negative_windows=uncertain_only,
        specificity_lower=_ratio(len(negatives) - fp_windows - uncertain_only, len(negatives)),
        specificity_upper=_ratio(len(negatives) - fp_windows, len(negatives)),
        full_sequence_reconstruction=False,
        full_sequence_reconstruction_status="not_equivalent_proven_false_junction")
    return result
