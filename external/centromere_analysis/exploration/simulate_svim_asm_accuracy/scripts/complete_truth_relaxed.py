"""Exploratory post-calibration closeness curves; strict scoring is unchanged."""

from itertools import combinations
import math
import edlib

from complete_truth_scoring import _apply, _distance, _edit, _net, _overlap, _size, _span


def _local_sequences(template, truth, event_index, calls, indices, max_span):
    event = truth[event_index]
    edits = [_edit(event)] + [calls[j]["edit"] for j in indices]
    left = max(0, min(e[0] for e in edits)-1)
    right = min(template.get_reference_length(event["chrom"]), max(e[1] for e in edits)+1)
    neighbors = {event_index}
    while True:
        old = (left, right, len(neighbors))
        for i, other in enumerate(truth):
            if other["chrom"] == event["chrom"] and _overlap(_span(other, True), (left, right)):
                neighbors.add(i)
                left = max(0, min(left, other["start"]-1))
                right = min(template.get_reference_length(event["chrom"]), max(right, other["end"]+1))
        if old == (left, right, len(neighbors)):
            break
    if right-left > max_span:
        return None, "reconstruction_span_bound"
    sequence = template.fetch(event["chrom"], left, right).upper()
    expected = _apply(sequence, left, [_edit(truth[i]) for i in neighbors])
    proposed = _apply(sequence, left, [_edit(truth[i]) for i in neighbors-{event_index}] +
                      [calls[j]["edit"] for j in indices])
    if expected is None or proposed is None:
        return None, "overlapping_or_unordered_edits"
    return (expected, proposed, left, right), None


def _assign_disjoint(options, max_states):
    order = sorted(options, key=lambda i: (len(options[i]), i))
    best, states, bounded = [], 0, False

    def visit(position, used, chosen):
        nonlocal best, states, bounded
        states += 1
        if states > max_states:
            bounded = True
            return
        if len(chosen) > len(best):
            best = list(chosen)
        if position == len(order) or len(chosen)+len(order)-position <= len(best):
            return
        event_index = order[position]
        for option in options[event_index]:
            indices = set(option["call_indices"])
            if not indices & used:
                visit(position+1, used | indices, chosen+[(event_index, option)])
                if bounded:
                    return
        visit(position+1, used, chosen)

    visit(0, set(), [])
    return best, bounded, states


def score_relaxed(template, truth, calls, strict_result, *,
                  breakpoint_tolerances=(100, 500, 1000), relative_size_tolerance=0.1,
                  relative_edit_thresholds=(0, 0.001, 0.01, 0.05), candidate_flank=1000,
                  max_subset_calls=8, max_comparisons=5000, max_assignment_states=5000,
                  max_reconstruction_span=2_000_000, minimum_distance_cap=1000):
    if (max_subset_calls < 1 or max_comparisons < 1 or max_assignment_states < 1 or
            candidate_flank < 0 or relative_size_tolerance < 0 or
            minimum_distance_cap < 0 or
            any(x < 0 for x in breakpoint_tolerances) or any(x < 0 for x in relative_edit_thresholds)):
        raise ValueError("Invalid relaxed-search bounds")
    reserved = {row.get("index", j) for j, row in enumerate(strict_result["calls"]) if row["status"] == "tp"}
    recovered_ids = {row["id"] for row in strict_result["events"] if row["recovered"]}
    strict_tp = len(recovered_ids)
    events, candidate_sets, comparisons = [], {}, 0
    for i, event in enumerate(truth):
        if event["id"] in recovered_ids:
            continue
        target = _edit(event)
        size = _size(target)
        distance_cap = max(minimum_distance_cap, math.ceil(max(relative_edit_thresholds, default=0)*size))
        eligible = [j for j, call in enumerate(calls) if j not in reserved and
                    call["chrom"] == event["chrom"] and
                    _distance(_span(call), (event["start"], event["end"])) <= candidate_flank]
        supported = [j for j in eligible if calls[j].get("edit") is not None]
        supported.sort(key=lambda j: (abs(calls[j]["edit"][0]-target[0]) +
                                      abs(calls[j]["edit"][1]-target[1]), j))
        bounded = len(supported) > max_subset_calls
        selected = supported[:max_subset_calls]
        local_options, failed_reconstructions = [], 0
        for n in range(1, len(selected)+1):
            for subset in combinations(selected, n):
                subset = tuple(sorted(subset))
                if comparisons >= max_comparisons:
                    bounded = True
                    break
                comparisons += 1
                reconstructed, failure = _local_sequences(template, truth, i, calls, subset,
                                                          max_reconstruction_span)
                if failure:
                    failed_reconstructions += 1
                    if failure == "reconstruction_span_bound":
                        bounded = True
                    continue
                expected, observed, left, right = reconstructed
                length_lower_bound = abs(len(expected)-len(observed))
                distance = (-1 if length_lower_bound > distance_cap else
                            edlib.align(expected, observed, mode="NW", task="distance",
                                        k=distance_cap)["editDistance"])
                censored = distance < 0
                lower_bound = max(distance_cap+1, length_lower_bound) if censored else distance
                exact_distance = None if censored else distance
                start_shift = min(calls[j]["edit"][0] for j in subset)-target[0]
                end_shift = max(calls[j]["edit"][1] for j in subset)-target[1]
                net_error = sum(_net(calls[j]["edit"]) for j in subset)-_net(target)
                local_options.append(dict(call_indices=list(subset),
                    call_ids=[calls[j].get("id", f"call_{j}") for j in subset],
                    called_types=[calls[j]["kind"] for j in subset],
                    dna_edit_distance=exact_distance, dna_edit_distance_censored=censored,
                    dna_edit_distance_lower_bound=lower_bound, edit_distance_cap=distance_cap,
                    relative_dna_edit_distance=None if censored else distance/size,
                    relative_dna_edit_distance_lower_bound=lower_bound/size,
                    relative_to_local_alt_length=None if censored else distance/max(len(expected), len(observed), 1),
                    relative_to_local_alt_length_lower_bound=lower_bound/max(len(expected), len(observed), 1),
                    start_shift_bp=start_shift, end_shift_bp=end_shift,
                    max_breakpoint_displacement_bp=max(abs(start_shift), abs(end_shift)),
                    net_length_error_bp=net_error, relative_net_length_error=abs(net_error)/size,
                    summed_call_edit_size_bp=sum(_size(calls[j]["edit"]) for j in subset),
                    truth_event_size_bp=size, truth_local_alt_length=len(expected),
                    called_local_alt_length=len(observed), reconstruction_start=left, reconstruction_end=right))
            if comparisons >= max_comparisons:
                break
        local_options.sort(key=lambda row: (row["relative_dna_edit_distance_lower_bound"],
                                            row["max_breakpoint_displacement_bp"],
                                            row["relative_net_length_error"], len(row["call_indices"]),
                                            row["call_indices"]))
        candidate_sets[i] = local_options
        row = dict(id=event["id"], chrom=event["chrom"], kind=event["kind"],
                   truth_event_size_bp=size, eligible_call_count=len(eligible),
                   supported_candidate_count=len(supported), unsupported_candidate_count=len(eligible)-len(supported),
                   candidate_call_indices=selected, valid_subsets=len(local_options),
                   failed_reconstructions=failed_reconstructions, search_bounded=bounded,
                   interpretation="exploratory closeness; strict result unchanged",
                   candidate_ranking="DNA-distance lower bound; censored candidates are not proven nearest",
                   status="candidate_quantified" if local_options else "no_reconstructable_candidate")
        if local_options:
            row.update(local_options[0])
        events.append(row)
    curves = []
    for bp in breakpoint_tolerances:
        for threshold in relative_edit_thresholds:
            options = {i: [candidate for candidate in group
                           if candidate["max_breakpoint_displacement_bp"] <= bp and
                           candidate["relative_net_length_error"] <= relative_size_tolerance and
                           not candidate["dna_edit_distance_censored"] and
                           candidate["relative_dna_edit_distance"] <= threshold]
                       for i, group in candidate_sets.items()}
            options = {i: group for i, group in options.items() if group}
            assigned, assignment_bounded, states = _assign_disjoint(options, max_assignment_states)
            assigned_calls = {j for _, option in assigned for j in option["call_indices"]}
            combined = strict_tp+len(assigned)
            curves.append(dict(breakpoint_tolerance_bp=bp, relative_net_size_tolerance=relative_size_tolerance,
                relative_dna_edit_threshold=threshold, n_truth=len(truth), n_calls=len(calls),
                strict_recovered_events=strict_tp, additional_approximate_events=len(assigned),
                combined_detected_events=combined,
                approximate_recovery_fraction=combined/len(truth) if truth else None,
                strict_reserved_call_records=len(reserved), additional_assigned_call_records=len(assigned_calls),
                unmatched_call_records=len(calls)-len(reserved)-len(assigned_calls),
                apparent_record_precision_under_relaxed_criterion=(len(reserved)+len(assigned_calls))/len(calls) if calls else None,
                search_bounded=any(row["search_bounded"] for row in events),
                assignment_bounded=assignment_bounded, assignment_states=states,
                matches=[dict(truth_id=truth[i]["id"], **option) for i, option in assigned]))
    return dict(events=events, curves=curves, comparisons=comparisons,
                conventions=dict(primary="unchanged strict complete-truth scoring",
                    analysis="exploratory after calibration; not sequence certification",
                    threshold_application="Only residual unmatched truth is relaxed; all strict matches remain reserved regardless of these breakpoint tolerances",
                    dna_distance="bounded local edlib NW; exact within cap, otherwise explicit censored lower bound",
                    distance_cap_rule="max(minimum_distance_cap, ceil(maximum DNA threshold * true event size))",
                    minimum_distance_cap=minimum_distance_cap,
                    dna_denominator="true event edit size, not full reconstructed chromosome/window length",
                    size_criterion="absolute net-length error divided by true event edit size",
                    candidate_rule="same chromosome and interval overlap within fixed flank; no genome-wide allele search",
                    candidate_flank=candidate_flank, max_subset_calls=max_subset_calls,
                    max_comparisons=max_comparisons, max_assignment_states=max_assignment_states,
                    negative_controls="No structural truth means all call records remain off-ledger; neutral pairs do not create TP events"))
