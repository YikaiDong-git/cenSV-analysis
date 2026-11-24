"""Call-led diagnostics for strict complete-truth scoring."""
from collections import Counter
from functools import lru_cache

import edlib

from complete_truth_scoring import _apply, _distance, _edit, _intervals, _net, _overlap, _size, _span

ALGORITHM_VERSION = 'net_pruned_cached_diagnostics_v2'
_MAX_COMPARISONS_PER_CALL = 100_000
_MAX_ENUMERATION_STATES = 2_000_000
_MAX_ALIGNMENT_WORK = 6_400_000_000


class DiagnosticSearchLimit(RuntimeError):
    pass


def _common_prefix(a, b):
    size = min(len(a), len(b)); pos = 0; block = 8192
    while pos + block <= size and a[pos:pos + block] == b[pos:pos + block]:
        pos += block
    lo, hi = 0, min(block, size - pos)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if a[pos:pos + mid] == b[pos:pos + mid]:
            lo = mid
        else:
            hi = mid - 1
    return pos + lo


def _compatible(edits):
    ordered = sorted(edits, key=lambda e: (e[0], e[1]))
    cursor = -1; insertions = set()
    for start, end, _ in ordered:
        if start < cursor or end < start:
            return False
        if start == end:
            if start in insertions:
                return False
            insertions.add(start)
        cursor = end
    return True


def _net_range(calls, indices, focal):
    optional = [_net(calls[i]['edit']) for i in indices if i != focal]
    fixed = _net(calls[focal]['edit'])
    return fixed + sum(min(x, 0) for x in optional), fixed + sum(max(x, 0) for x in optional)


def _groups_in_net_range(calls, indices, focal, target, allowance, states,
                         max_states=_MAX_ENUMERATION_STATES):
    optional = [i for i in indices if i != focal]
    nets = [_net(calls[i]['edit']) for i in optional]
    n = len(optional); low = [0] * (n + 1); high = [0] * (n + 1)
    for pos in range(n - 1, -1, -1):
        low[pos] = low[pos + 1] + min(nets[pos], 0)
        high[pos] = high[pos + 1] + max(nets[pos], 0)
    def visit(pos, net, chosen):
        states[0] += 1
        if states[0] > max_states:
            raise DiagnosticSearchLimit('Diagnostic subset-state limit')
        if net + high[pos] < target - allowance or net + low[pos] > target + allowance:
            return
        if pos == n:
            yield tuple(sorted(chosen, key=lambda i: (_span(calls[i]), i)))
            return
        j = optional[pos]
        if _compatible([calls[i]['edit'] for i in (*chosen, j)]):
            yield from visit(pos + 1, net + nets[pos], (*chosen, j))
        yield from visit(pos + 1, net, chosen)
    yield from visit(0, _net(calls[focal]['edit']), (focal,))


def _distance_bounded(a, b, cap):
    prefix = _common_prefix(a, b)
    a, b = a[prefix:], b[prefix:]
    suffix = _common_prefix(a[::-1], b[::-1])
    if suffix:
        a, b = a[:-suffix], b[:-suffix]
    return _distance_trimmed(a, b, cap)


@lru_cache(maxsize=64)
def _distance_trimmed(a, b, cap):
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
                         max_span=2000000, max_calls=32, distance_cap=200000):
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
    cache = {}; cache_hits = 0
    def reconstruct(event_index, indices):
        nonlocal cache_hits
        key = (event_index, tuple(sorted(indices)))
        if key not in cache:
            cache[key] = _reconstruct(template, truth, event_index, calls, key[1], max_span, distance_cap)
        else:
            cache_hits += 1
        result, failure = cache[key]
        return (dict(result) if result else None), failure
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
        states = [0]; exclusions = []; evaluated = 0; invalid_groups = 0
        for _, event_index in candidates:
            nearby = [j for j, c in enumerate(calls) if c["chrom"] == call["chrom"] and
                      c.get("edit") is not None and by_index[j]["status"] in ("fp", "unresolved") and
                      _distance(_span(c), _span(truth[event_index], True)) <= candidate_flank]
            target = _net(_edit(truth[event_index]))
            allowance = _size(_edit(truth[event_index])) // 2
            net_lower, net_upper = _net_range(calls, nearby, index)
            net_error_lower = max(net_lower - target, target - net_upper, 0)
            if net_error_lower > allowance:
                exclusions.append(dict(event_id=truth[event_index]['id'], required_net_bp=target,
                    achievable_net_lower_bp=net_lower, achievable_net_upper_bp=net_upper,
                    dna_error_lower_bound_bp=net_error_lower, allowed_dna_error_bp=allowance))
                result, failure = reconstruct(event_index, (index,))
                if result:
                    result.update(marginal_improvement_proven=result['improvement_proven'],
                        marginal_improvement_exact_proven=result['baseline_improvement_exact_proven'],
                        without_focal_dna_edit_distance=result['baseline_dna_edit_distance'],
                        without_focal_dna_edit_distance_lower_bound=result['baseline_dna_edit_distance_lower_bound'],
                        without_focal_dna_edit_distance_censored=result['baseline_dna_edit_distance_censored'])
                    choices.append((3, result['dna_edit_distance_lower_bound'], 1, event_index,
                                    (index,), 'near_truth_unattributed', result))
                continue
            nearby.sort(key=lambda j: (_distance(_span(calls[j]), span), j))
            if len(nearby) > max_calls:
                bounded = True
            selected = nearby[:max_calls]
            if index not in selected:
                selected[-1:] = [index]
            selected.sort(key=lambda j: (_span(calls[j]), j))
            groups = _groups_in_net_range(calls, selected, index, target, allowance, states)
            try:
                candidate_groups = list(groups)
            except DiagnosticSearchLimit:
                bounded = True
                candidate_groups = []
            candidate_groups = sorted(set(candidate_groups + [(index,)]), key=lambda g: (len(g), g))
            for group in candidate_groups:
                if row["comparison_count"] >= _MAX_COMPARISONS_PER_CALL:
                    bounded = True
                    break
                row["comparison_count"] += 1
                result, failure = reconstruct(event_index, group)
                if failure:
                    if failure != 'overlapping_or_unordered_edits':
                        bounded = True
                    else:
                        invalid_groups += 1
                    continue
                evaluated += 1
                if ((result["comparison_work_bounded"] or result["dna_edit_distance_censored"]) and
                    result['dna_edit_distance_lower_bound'] <= allowance):
                    bounded = True
                distance = result["dna_edit_distance"]
                marginal_distance = result["baseline_dna_edit_distance"]
                marginal_lower = result["baseline_dna_edit_distance_lower_bound"]
                marginal_censored = result["baseline_dna_edit_distance_censored"]
                marginal_proven = result["improvement_proven"]
                marginal_exact = result["baseline_improvement_exact_proven"]
                if len(group) > 1:
                    without, failure = reconstruct(event_index, tuple(j for j in group if j != index))
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
                        if without["comparison_work_bounded"] and distance is not None and not marginal_proven:
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
        row.update(net_length_exclusions=exclusions, enumeration_states=states[0],
                   evaluated_groups=evaluated, invalid_groups_rejected=invalid_groups,
                   explanation_search_complete=not bounded)
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
        row['explanation_resolution'] = ('net_length_excludes_all_candidates' if len(exclusions) == len(candidates)
            else 'supported_explanation_found' if row['category'] in ('near_correct_true_event', 'distorted_true_event')
            else 'search_incomplete' if bounded else 'no_supported_local_explanation')
        rows.append(row)
    counts = Counter(r["category"] for r in rows)
    return {"rows": rows, "summary": {"n_calls": len(calls), "category_counts": dict(counts),
            "search_bounded_records": sum(r["search_bounded"] for r in rows),
            "algorithm_version": ALGORITHM_VERSION, "group_cache_entries": len(cache),
            "group_cache_hits": cache_hits, "max_enumeration_states_per_call": _MAX_ENUMERATION_STATES,
            "candidate_flank_bp": candidate_flank, "max_span_bp": max_span, "max_calls": max_calls,
            "distance_cap": distance_cap, "max_comparisons_per_call": _MAX_COMPARISONS_PER_CALL,
            "max_alignment_work": _MAX_ALIGNMENT_WORK,
            "interpretation": "Non-disjoint focal-record diagnostics; groups can reuse records and events. Off-local means no event destination within the specified flank; source support can be closer. Strict labels and precision remain unchanged."}}
