"""Call-led diagnostics for strict complete-truth scoring."""
from collections import Counter
from itertools import combinations

import edlib

from complete_truth_scoring import _apply, _distance, _edit, _intervals, _net, _overlap, _size, _span

_MAX_COMPARISONS_PER_CALL = 64
_MAX_ALIGNMENT_WORK = 100_000_000


def _distance_bounded(a, b, cap):
    common = min(len(a), len(b))
    prefix = 0
    while prefix < common and a[prefix] == b[prefix]:
        prefix += 1
    a, b = a[prefix:], b[prefix:]
    common = min(len(a), len(b))
    suffix = 0
    while suffix < common and a[len(a) - suffix - 1] == b[len(b) - suffix - 1]:
        suffix += 1
    if suffix:
        a, b = a[:-suffix], b[:-suffix]
    if not a or not b:
        exact = max(len(a), len(b))
        return exact, exact, False, False
    lower = abs(len(a) - len(b))
    if lower > cap:
        return None, lower, True, False
    bands = sorted({lower, cap, *(k for k in (1, 2, 5, 10, 20, 50, 100, 200, 500,
                                                 1000, 2000, 5000, 10000)
                                   if lower <= k <= cap)})
    for band in bands:
        if max(len(a), len(b)) * min(band + 1, max(len(a), len(b))) > _MAX_ALIGNMENT_WORK:
            return None, lower, True, True
        distance = edlib.align(a, b, mode="NW", task="distance", k=band)["editDistance"]
        if distance >= 0:
            return distance, distance, False, False
        lower = max(lower, band + 1)
    return None, lower, True, False


def _reconstruct(template, truth, event_index, calls, indices, max_span, distance_cap):
    event = truth[event_index]
    edits = [_edit(event)] + [calls[i]["edit"] for i in indices]
    chrom = event["chrom"]
    length = template.get_reference_length(chrom)
    left = max(0, min(e[0] for e in edits) - 1)
    right = min(length, max(e[1] for e in edits) + 1)
    neighbors = {event_index}
    while True:
        old = (left, right, len(neighbors))
        for i, other in enumerate(truth):
            if other["chrom"] == chrom and _overlap(_span(other, True), (left, right)):
                neighbors.add(i)
                left = max(0, min(left, int(other["start"]) - 1))
                right = min(length, max(right, int(other["end"]) + 1))
        if (left, right, len(neighbors)) == old:
            break
    if right - left > max_span:
        return None, "reconstruction_span_bound"
    reference = template.fetch(chrom, left, right).upper()
    expected = _apply(reference, left, [_edit(truth[i]) for i in neighbors])
    background = [_edit(truth[i]) for i in neighbors if i != event_index]
    baseline = _apply(reference, left, background)
    proposed = _apply(reference, left, background + [calls[i]["edit"] for i in indices])
    if expected is None or baseline is None or proposed is None:
        return None, "overlapping_or_unordered_edits"
    actual, actual_lower, actual_censored, actual_work = _distance_bounded(expected, proposed, distance_cap)
    base, base_lower, base_censored, base_work = _distance_bounded(expected, baseline, distance_cap)
    size = _size(_edit(event))
    if size <= 0:
        return None, "zero_size_truth"
    net_error = sum(_net(calls[i]["edit"]) for i in indices) - _net(_edit(event))
    improvement = actual is not None and ((base is not None and actual < base) or
                                          (base is None and actual < base_lower))
    effect_fraction_lower = (1 - actual / base_lower
                             if actual is not None and base_lower > 0 else None)
    supported_effect = effect_fraction_lower is not None and effect_fraction_lower >= .5
    return dict(
        baseline_effect_fraction_recovered_lower_bound=effect_fraction_lower,
        baseline_effect_support_at_least_half=supported_effect,
        weak_positive_improvement=bool(improvement and not supported_effect),
        candidate_event_id=event["id"], group_call_ids=[calls[i].get("id", str(i)) for i in indices],
        dna_edit_distance=actual, dna_edit_distance_lower_bound=actual_lower,
        dna_edit_distance_censored=actual_censored,
        baseline_dna_edit_distance=base, baseline_dna_edit_distance_lower_bound=base_lower,
        baseline_dna_edit_distance_censored=base_censored,
        dna_error_fraction_of_event=None if actual is None else actual / size,
        dna_error_fraction_of_local_alt=None if actual is None else actual / max(len(expected), len(proposed), 1),
        event_edit_size_bp=size, truth_local_alt_length=len(expected), proposed_local_alt_length=len(proposed),
        start_displacement_bp=min(calls[i]["edit"][0] for i in indices) - int(event["start"]),
        end_displacement_bp=max(calls[i]["edit"][1] for i in indices) - int(event["end"]),
        net_length_error_bp=net_error, relative_net_length_error=abs(net_error) / size,
        improvement_proven=improvement,
        baseline_improvement_exact_proven=bool(improvement and base is not None),
        reconstruction_start=left, reconstruction_end=right,
        comparison_work_bounded=actual_work or base_work,
    ), None


def classify_call_errors(template, truth, calls, strict_result, *, candidate_flank=10000,
                         max_span=2000000, max_calls=8, distance_cap=20000):
    if candidate_flank < 0 or max_span < 1 or max_calls < 1 or distance_cap < 0:
        raise ValueError("Invalid diagnostic bounds")
    strict = strict_result["calls"]
    if len(strict) != len(calls):
        raise ValueError("Strict call rows must correspond one-to-one with calls")
    by_index = {row.get("index", i): row for i, row in enumerate(strict)}
    if set(by_index) != set(range(len(calls))):
        raise ValueError("Strict call indices are missing or duplicated")
    recovered_rows = {row["id"]: row for row in strict_result["events"] if row.get("recovered")}
    recovered = set(recovered_rows)
    rows = []
    for index, call in enumerate(calls):
        sr = by_index[index]
        row = dict(call_index=index, call_id=call.get("id", str(index)), chrom=call["chrom"],
                   called_type=call.get("kind"), strict_status=sr["status"],
                   category=None, representation_unassessed=False, search_bounded=False,
                   candidate_event_id=None, group_call_ids=[], comparison_count=0,
                   nearest_truth_distance_bp=None, nearest_truth_event_id=None,
                   nearest_truth_support_distance_bp=None, nearest_truth_support_event_id=None)
        if sr["status"] == "tp":
            row["category"] = "strict_match"
            rows.append(row)
            continue
        if not truth:
            row.update(category="null_control_false_call", representation_unassessed=call.get("edit") is None)
            rows.append(row)
            continue
        if sr["status"] not in ("fp", "unresolved") or call.get("edit") is None:
            row.update(category="unsupported_or_bounded", representation_unassessed=True)
            rows.append(row)
            continue
        span = _span(call)
        candidates = sorted(((_distance(span, _span(t, True)), i) for i, t in enumerate(truth)
                             if t["chrom"] == call["chrom"]), key=lambda x: (x[0], x[1]))
        if candidates:
            row["nearest_truth_distance_bp"], ti = candidates[0]
            row["nearest_truth_event_id"] = truth[ti]["id"]
        support = sorted((min(_distance(span, interval) for interval in _intervals(t)), i)
                         for i, t in enumerate(truth) if t["chrom"] == call["chrom"])
        if support:
            row["nearest_truth_support_distance_bp"], ti = support[0]
            row["nearest_truth_support_event_id"] = truth[ti]["id"]
        if sr.get("reason") == "duplicate_equivalent_call":
            row.update(category="redundant_exact_representation",
                       attribution_source="strict_scorer_duplicate_equivalent_call",
                       strict_duplicate_truth_ids=sr.get("matched_truth_ids", []))
            rows.append(row)
            continue
        candidates = [(d, i) for d, i in candidates if d <= candidate_flank]
        if not candidates:
            row["category"] = "off_truth_no_local_target"
            rows.append(row)
            continue
        choices = []
        bounded = False
        for _, event_index in candidates:
            nearby = [j for j, c in enumerate(calls) if c["chrom"] == call["chrom"] and
                      c.get("edit") is not None and by_index[j]["status"] in ("fp", "unresolved") and
                      _distance(_span(c), _span(truth[event_index], True)) <= candidate_flank]
            nearby.sort(key=lambda j: (_distance(_span(calls[j]), span), j))
            if len(nearby) > max_calls:
                bounded = True
            selected = nearby[:max_calls]
            if index not in selected:
                selected[-1:] = [index]
            selected.sort(key=lambda j: (_span(calls[j]), j))
            groups = [group for n in range(1, len(selected) + 1)
                      for group in combinations(selected, n) if index in group]
            for group in groups:
                if row["comparison_count"] >= _MAX_COMPARISONS_PER_CALL:
                    bounded = True
                    break
                row["comparison_count"] += 1
                result, failure = _reconstruct(template, truth, event_index, calls, group, max_span, distance_cap)
                if failure:
                    bounded = True
                    continue
                if result["comparison_work_bounded"] or result["dna_edit_distance_censored"]:
                    bounded = True
                distance = result["dna_edit_distance"]
                marginal_distance = result["baseline_dna_edit_distance"]
                marginal_lower = result["baseline_dna_edit_distance_lower_bound"]
                marginal_censored = result["baseline_dna_edit_distance_censored"]
                marginal_proven = result["improvement_proven"]
                marginal_exact = result["baseline_improvement_exact_proven"]
                if len(group) > 1:
                    without, failure = _reconstruct(template, truth, event_index, calls,
                        tuple(j for j in group if j != index), max_span, distance_cap)
                    if failure:
                        marginal_proven = marginal_exact = False
                        marginal_distance = marginal_lower = None
                        marginal_censored = True
                        bounded = True
                    else:
                        marginal_distance = without["dna_edit_distance"]
                        marginal_lower = without["dna_edit_distance_lower_bound"]
                        marginal_censored = without["dna_edit_distance_censored"]
                        marginal_proven = distance is not None and (
                            (marginal_distance is not None and distance < marginal_distance) or
                            (marginal_distance is None and distance < marginal_lower))
                        marginal_exact = bool(marginal_proven and marginal_distance is not None)
                        if without["comparison_work_bounded"]:
                            bounded = True
                result.update(marginal_improvement_proven=marginal_proven,
                              marginal_improvement_exact_proven=marginal_exact,
                              without_focal_dna_edit_distance=marginal_distance,
                              without_focal_dna_edit_distance_lower_bound=marginal_lower,
                              without_focal_dna_edit_distance_censored=marginal_censored)
                if len(group) == 1 and distance == 0 and truth[event_index]["id"] in recovered:
                    category, rank = "redundant_exact_representation", 0
                elif (distance is not None and result["improvement_proven"] and marginal_proven and result["baseline_effect_support_at_least_half"] and
                      result["dna_error_fraction_of_event"] <= .01 and
                      result["relative_net_length_error"] <= .10):
                    category, rank = "near_correct_true_event", 1
                elif (distance is not None and result["improvement_proven"] and marginal_proven and result["baseline_effect_support_at_least_half"] and
                       .01 < result["dna_error_fraction_of_event"] <= .50):
                    category, rank = "distorted_true_event", 2
                else:
                    category, rank = "near_truth_unattributed", 3
                choices.append((rank, result["dna_edit_distance_lower_bound"], len(group), event_index,
                                tuple(group), category, result))
            if row["comparison_count"] >= _MAX_COMPARISONS_PER_CALL:
                bounded = True
                break
        row["search_bounded"] = bounded
        row["group_attribution_scope"] = "focal record only; other members require their own marginal checks"
        if choices:
            best = min(choices, key=lambda x: x[:5])
            row.update(best[-1]); row["category"] = best[-2]
            matched_event = recovered_rows.get(row["candidate_event_id"])
            row["candidate_event_already_recovered"] = matched_event is not None
            row["candidate_event_strict_call_ids"] = (matched_event.get("matched_call_ids", [])
                                                       if matched_event is not None else [])
            if row["category"] == "near_truth_unattributed" and bounded:
                row["category"] = "unsupported_or_bounded"
                row["representation_unassessed"] = True
        elif bounded:
            row.update(category="unsupported_or_bounded", representation_unassessed=True)
        else:
            row["category"] = "near_truth_unattributed"
        rows.append(row)
    counts = Counter(r["category"] for r in rows)
    return {"rows": rows, "summary": {"n_calls": len(calls), "category_counts": dict(counts),
            "search_bounded_records": sum(r["search_bounded"] for r in rows),
            "candidate_flank_bp": candidate_flank, "max_span_bp": max_span, "max_calls": max_calls,
            "distance_cap": distance_cap, "max_comparisons_per_call": _MAX_COMPARISONS_PER_CALL,
            "max_alignment_work": _MAX_ALIGNMENT_WORK,
            "interpretation": "Non-disjoint focal-record diagnostics; groups can reuse records and events. Off-local means no event destination within the specified flank; source support can be closer. Strict labels and precision remain unchanged."}}
