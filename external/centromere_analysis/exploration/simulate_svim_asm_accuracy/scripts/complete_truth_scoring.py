"""Score a complete structural-edit ledger on a known substitution template."""

from collections import Counter
from itertools import combinations


def _apply(sequence, origin, edits):
    parts, cursor, insertions = [], origin, set()
    for start, end, allele in sorted(edits, key=lambda e: (e[0], e[1])):
        if start < cursor or end < start or end > origin + len(sequence) or start < origin:
            return None
        if start == end:
            if start in insertions:
                return None
            insertions.add(start)
        parts.extend((sequence[cursor-origin:start-origin], allele))
        cursor = end
    parts.append(sequence[cursor-origin:])
    return "".join(parts)


def _edit(event):
    return (int(event["start"]), int(event["end"]), event["alt"].upper())


def _span(item, is_truth=False):
    edit = _edit(item) if is_truth else item.get("edit")
    if edit is not None:
        return int(edit[0]), int(edit[1])
    start = int(item.get("start", item.get("pos", 1) - 1))
    return start, max(start, int(item.get("end", start)))


def _intervals(event):
    intervals = [_span(event, True)]
    if event.get("support_start") is not None and event.get("support_end") is not None:
        intervals.append((event["support_start"], event["support_end"]))
    if event.get("source_start") is not None and event.get("source_end") is not None:
        intervals.append((event["source_start"], event["source_end"]))
    return intervals


def _distance(a, b):
    return max(a[0] - b[1], b[0] - a[1], 0)


def _overlap(a, b):
    if a[0] == a[1]:
        return b[0] <= a[0] < b[1]
    return a[0] < b[1] and a[1] > b[0]


def _net(edit):
    return len(edit[2]) - (edit[1] - edit[0])


def _size(edit):
    return max(edit[1] - edit[0], len(edit[2]))


def _effect(edit):
    return "INS" if edit[0] == edit[1] else "DEL" if not edit[2] else "REPLACEMENT"


def _expected_type(kind):
    return {"TANDEM": "DUP:TANDEM", "DONOR_INS": "INS", "DISPERSED": "DUP:INT",
            "DISPERSED_COPY": "DUP:INT"}.get(kind, kind)


def _compatible_type(truth_kind, called_kind):
    expected = _expected_type(truth_kind)
    return expected == called_kind or (expected in ("DUP:TANDEM", "DUP:INT") and called_kind == "INS")


def _sequence_pieces(edits, length):
    pieces, cursor, insertions = [], 0, set()
    for start, end, allele in sorted(edits, key=lambda e: (e[0], e[1])):
        if not cursor <= start <= end <= length:
            return None
        if start == end:
            if start in insertions:
                return None
            insertions.add(start)
        if start > cursor:
            pieces.append(("reference", cursor, start))
        if allele:
            pieces.append(("literal", allele, len(allele)))
        cursor = end
    if cursor < length:
        pieces.append(("reference", cursor, length))
    return pieces


def _piece_sequences_equal(template, chrom, first, second):
    i = j = ai = bj = 0
    while i < len(first) and j < len(second):
        a, b = first[i], second[j]
        alen = a[2]-a[1] if a[0] == "reference" else a[2]
        blen = b[2]-b[1] if b[0] == "reference" else b[2]
        count = min(alen-ai, blen-bj)
        identical_reference = a[0] == b[0] == "reference" and a[1]+ai == b[1]+bj
        if not identical_reference:
            count = min(count, 1_000_000)
            aa = (template.fetch(chrom, a[1]+ai, a[1]+ai+count).upper()
                  if a[0] == "reference" else a[1][ai:ai+count])
            bb = (template.fetch(chrom, b[1]+bj, b[1]+bj+count).upper()
                  if b[0] == "reference" else b[1][bj:bj+count])
            if aa != bb:
                return False
        ai += count
        bj += count
        if ai == alen:
            i += 1
            ai = 0
        if bj == blen:
            j += 1
            bj = 0
    return i == len(first) and j == len(second)


def _complete_reconstruction(template, truth, calls):
    if any(c.get("edit") is None for c in calls):
        return dict(status="unresolved_unsupported_calls", equivalent=None)
    for chrom in sorted({t["chrom"] for t in truth} | {c["chrom"] for c in calls}):
        length = template.get_reference_length(chrom)
        expected = _sequence_pieces([_edit(t) for t in truth if t["chrom"] == chrom], length)
        observed = _sequence_pieces([c["edit"] for c in calls if c["chrom"] == chrom], length)
        if expected is None:
            raise ValueError("Truth ledger has ambiguous edit ordering")
        if observed is None:
            return dict(status="unresolved_overlapping_or_unordered_calls", equivalent=None)
        if not _piece_sequences_equal(template, chrom, expected, observed):
            return dict(status="not_equivalent", equivalent=False)
    return dict(status="equivalent", equivalent=True)


def _equivalent(template, truth, selected_truth, calls, selected_calls):
    target = [truth[i] for i in selected_truth]
    proposed = [calls[j] for j in selected_calls]
    if any(c.get("edit") is None for c in proposed):
        return False
    chrom = target[0]["chrom"]
    if any(t["chrom"] != chrom for t in target) or any(c["chrom"] != chrom for c in proposed):
        return False
    edits = [_edit(t) for t in target] + [c["edit"] for c in proposed]
    left = max(0, min(e[0] for e in edits) - 1)
    right = min(template.get_reference_length(chrom), max(e[1] for e in edits) + 1)
    neighbors = set(selected_truth)
    while True:
        previous = (left, right, len(neighbors))
        for i, event in enumerate(truth):
            if event["chrom"] == chrom and _overlap(_span(event, True), (left, right)):
                neighbors.add(i)
                left = max(0, min(left, event["start"] - 1))
                right = min(template.get_reference_length(chrom), max(right, event["end"] + 1))
        if previous == (left, right, len(neighbors)):
            break
    sequence = template.fetch(chrom, left, right).upper()
    expected = _apply(sequence, left, [_edit(truth[i]) for i in neighbors])
    alternate = _apply(sequence, left, [_edit(truth[i]) for i in neighbors-set(selected_truth)] +
                       [c["edit"] for c in proposed])
    return expected is not None and alternate is not None and expected == alternate


def _maximum_matching(edges, quality=None):
    import numpy as np
    from scipy.optimize import linear_sum_assignment

    ts = sorted(edges)
    cs = sorted({c for values in edges.values() for c in values})
    if not ts or not cs:
        return {}
    ci = {c: i for i, c in enumerate(cs)}
    costs = quality or {(t, c): rank for t, values in edges.items() for rank, c in enumerate(values)}
    highest = max(costs.values(), default=0)
    missing = (highest + 1) * (len(ts) + 1)
    forbidden = (missing + 1) * (len(ts) + 1)
    matrix = np.full((len(ts), len(cs) + len(ts)), forbidden, dtype=np.int64)
    matrix[:, len(cs):] = missing
    for row, t in enumerate(ts):
        for c in edges[t]:
            matrix[row, ci[c]] = costs[(t, c)]
    rows, columns = linear_sum_assignment(matrix)
    return {ts[r]: cs[c] for r, c in zip(rows, columns)
            if c < len(cs) and matrix[r, c] < missing}


def _overlapping_call_indices(calls, indices):
    ambiguous = set()
    for a, b in combinations(sorted(indices), 2):
        x, y = calls[a].get("edit"), calls[b].get("edit")
        if x is None or y is None or calls[a]["chrom"] != calls[b]["chrom"]:
            continue
        xi, yi = x[0] == x[1], y[0] == y[1]
        conflict = (x[0] == y[0] if xi and yi else
                    y[0] < x[0] < y[1] if xi else
                    x[0] < y[0] < x[1] if yi else
                    x[0] < y[1] and y[0] < x[1])
        if conflict:
            ambiguous.update((a, b))
    return ambiguous


def score_complete_truth(template_fasta, truth_events, calls, negative_windows, *,
                         max_representation_shift=1_000_000, max_component_calls=12,
                         compound_flank=10_000, max_component_truth=4,
                         max_compound_comparisons=10_000, approximate_distance=100,
                         approximate_length_fraction=0.1, negative_buffer=10_000):
    truth = [dict(t, alt=t["alt"].upper()) for t in truth_events]
    call_data = [dict(c, edit=(int(c["edit"][0]), int(c["edit"][1]), c["edit"][2].upper())
                      if c.get("edit") is not None else None) for c in calls]
    if len({t["id"] for t in truth}) != len(truth):
        raise ValueError("Truth IDs must be unique")
    for event in truth:
        start, end, allele = _edit(event)
        if not 0 <= start <= end <= template_fasta.get_reference_length(event["chrom"]):
            raise ValueError(f"Invalid truth interval: {event['id']}")
        if template_fasta.fetch(event["chrom"], start, end).upper() == allele:
            raise ValueError(f"Truth event has no sequence effect: {event['id']}")
    for chrom in {t["chrom"] for t in truth}:
        previous_end, insertion_starts = -1, set()
        for event in sorted((t for t in truth if t["chrom"] == chrom), key=lambda t: (t["start"], t["end"])):
            if event["start"] < previous_end:
                raise ValueError("Truth edits overlap; encode a replacement as one edit")
            if event["start"] == event["end"]:
                if event["start"] in insertion_starts:
                    raise ValueError("Multiple truth insertions at one boundary require a combined allele")
                insertion_starts.add(event["start"])
            previous_end = event["end"]
    event_rows = [dict(id=t["id"], chrom=t["chrom"], kind=t["kind"], start=t["start"], end=t["end"],
                       size_bp=_size(_edit(t)), inserted_bp=len(t["alt"]), deleted_bp=t["end"]-t["start"],
                       status="not_recovered", recovered=False, exact=False, matched_call_ids=[],
                       breakpoint_shift_bp=None, end_shift_bp=None, size_error_bp=None,
                       net_length_error_bp=None, called_type=None, type_match=None,
                       literal_type_match=None, type_compatible=None,
                       localization_within_100bp=None, approximate_match=False,
                       closest_candidate_call_id=None, closest_candidate_start_shift_bp=None,
                       closest_candidate_end_shift_bp=None, closest_candidate_size_error_bp=None,
                       closest_candidate_net_length_error_bp=None) for t in truth]
    call_rows = [dict(id=c.get("id", f"call_{j}"), index=j, chrom=c["chrom"], kind=c["kind"],
                      start=_span(c)[0], end=_span(c)[1], status="fp", reason="unmatched",
                      matched_truth_ids=[]) for j, c in enumerate(call_data)]
    individual, approximate, near_truth, near_calls = {}, {}, {}, {}
    for i, event in enumerate(truth):
        te = _edit(event)
        edges, app = [], []
        for j, call in enumerate(call_data):
            if call["chrom"] != event["chrom"]:
                continue
            ce = call.get("edit")
            close = min(_distance(_span(call), interval) for interval in _intervals(event)) <= compound_flank
            equiv = False
            if ce is not None and abs(ce[0]-te[0]) <= max_representation_shift:
                if _net(ce) == _net(te):
                    equiv = ce == te or _equivalent(template_fasta, truth, {i}, call_data, {j})
                if (_effect(ce) == _effect(te) and abs(ce[0]-te[0]) <= approximate_distance and
                        abs(_size(ce)-_size(te)) <= max(10, _size(te)*approximate_length_fraction)):
                    app.append(j)
            if equiv:
                edges.append(j)
            if close or equiv:
                near_truth.setdefault(i, set()).add(j)
                near_calls.setdefault(j, set()).add(i)
        individual[i] = sorted(edges, key=lambda j: (call_data[j]["edit"] != te,
                                                      abs(call_data[j]["edit"][0]-te[0]), j))
        approximate[i] = sorted(app, key=lambda j: abs(call_data[j]["edit"][0]-te[0]))
        candidates = [j for j in near_truth.get(i, ()) if call_data[j].get("edit") is not None]
        if candidates:
            closest = min(candidates, key=lambda j: (abs(call_data[j]["edit"][0]-te[0]),
                                                     abs(_size(call_data[j]["edit"])-_size(te)), j))
            ce = call_data[closest]["edit"]
            event_rows[i].update(closest_candidate_call_id=call_rows[closest]["id"],
                                 closest_candidate_start_shift_bp=ce[0]-te[0],
                                 closest_candidate_end_shift_bp=ce[1]-te[1],
                                 closest_candidate_size_error_bp=_size(ce)-_size(te),
                                 closest_candidate_net_length_error_bp=_net(ce)-_net(te))
    approximate_matches = _maximum_matching(approximate)
    for i in approximate_matches:
        event_rows[i]["approximate_match"] = True
    displacement = {(i, j): abs(call_data[j]["edit"][0]-truth[i]["start"]) +
                    abs(call_data[j]["edit"][1]-truth[i]["end"])
                    for i, edges in individual.items() for j in edges}
    exact_penalty = (max(displacement.values(), default=0) + 1) * (len(truth) + 1)
    quality = {(i, j): distance + exact_penalty * (call_data[j]["edit"] != _edit(truth[i]))
               for (i, j), distance in displacement.items()}
    matched = _maximum_matching(individual, quality)
    used_truth, used_calls, groups = set(), set(), []

    def assign(ts, cs, status):
        ts, cs = set(ts), set(cs)
        used_truth.update(ts)
        used_calls.update(cs)
        ids = [call_rows[j]["id"] for j in sorted(cs)]
        tids = [truth[i]["id"] for i in sorted(ts)]
        groups.append(dict(status=status, truth_ids=tids, call_ids=ids,
                           n_truth=len(ts), n_calls=len(cs)))
        for i in ts:
            event_rows[i].update(status=status, recovered=True, matched_call_ids=ids)
        for j in cs:
            call_rows[j].update(status="tp", reason=status, matched_truth_ids=tids)
        if len(ts) == len(cs) == 1:
            i, j = next(iter(ts)), next(iter(cs))
            te, ce = _edit(truth[i]), call_data[j]["edit"]
            event_rows[i].update(exact=ce == te, breakpoint_shift_bp=ce[0]-te[0], end_shift_bp=ce[1]-te[1],
                                 size_error_bp=_size(ce)-_size(te), net_length_error_bp=_net(ce)-_net(te),
                                 called_type=call_data[j]["kind"],
                                 type_match=_expected_type(truth[i]["kind"]) == call_data[j]["kind"],
                                 literal_type_match=_expected_type(truth[i]["kind"]) == call_data[j]["kind"],
                                 type_compatible=_compatible_type(truth[i]["kind"], call_data[j]["kind"]),
                                 localization_within_100bp=max(abs(ce[0]-te[0]), abs(ce[1]-te[1])) <= 100)

    for i, j in matched.items():
        assign({i}, {j}, "exact" if _edit(truth[i]) == call_data[j]["edit"] else "sequence_equivalent")
    visited, components = set(), []
    for seed in near_truth:
        if seed in visited:
            continue
        ts, cs, pending = set(), set(), [seed]
        while pending:
            i = pending.pop()
            if i in ts:
                continue
            ts.add(i)
            for j in near_truth.get(i, ()):
                cs.add(j)
                pending.extend(near_calls.get(j, set())-ts)
        visited.update(ts)
        remaining_t, remaining_c = ts-used_truth, cs-used_calls
        unsupported = {j for j in remaining_c if call_data[j]["edit"] is None}
        supported = sorted(remaining_c-unsupported)
        bounded = len(supported) > max_component_calls or len(remaining_t) > max_component_truth
        comparisons = 0
        if remaining_t and not bounded:
            targets = [set(s) for n in range(len(remaining_t), 0, -1)
                       for s in combinations(sorted(remaining_t), n)]
            for target in targets:
                if target & used_truth:
                    continue
                available = [j for j in supported if j not in used_calls]
                net = sum(_net(_edit(truth[i])) for i in target)
                found = None
                for n in range(1, len(available)+1):
                    for subset in combinations(available, n):
                        if sum(_net(call_data[j]["edit"]) for j in subset) != net:
                            continue
                        comparisons += 1
                        if comparisons > max_compound_comparisons:
                            bounded = True
                            break
                        if _equivalent(template_fasta, truth, target, call_data, set(subset)):
                            found = set(subset)
                            break
                    if found is not None or bounded:
                        break
                if found is not None:
                    status = ("split_equivalent" if len(target) == 1 and len(found) > 1 else
                              "merged_equivalent" if len(target) > 1 and len(found) == 1 else "compound_equivalent")
                    assign(target, found, status)
                if bounded:
                    break
        overlapping = _overlapping_call_indices(call_data, cs)
        duplicates = {j for j in cs-used_calls if any(j in edges for edges in individual.values())}
        unresolved_overlap = (overlapping-used_calls)-duplicates
        unresolved_c = unsupported | unresolved_overlap | ((cs-used_calls) if bounded else set())
        for j in unresolved_c-used_calls:
            reason = ("unsupported_call" if j in unsupported else "overlapping_or_unordered_calls"
                      if j in unresolved_overlap else "compound_search_bound")
            call_rows[j].update(status="unresolved", reason=reason)
        unresolved_t = ts-used_truth if unresolved_c-used_calls else set()
        for i in unresolved_t:
            event_rows[i]["status"] = "unresolved"
        components.append(dict(truth_ids=[truth[i]["id"] for i in sorted(ts)],
                               call_ids=[call_rows[j]["id"] for j in sorted(cs)],
                               search_bounded=bounded, comparisons=comparisons,
                               unsupported_calls=len(unsupported), overlapping_calls=len(overlapping)))
    for j, row in enumerate(call_rows):
        if row["status"] == "fp":
            if any(j in edges for edges in individual.values()):
                row["reason"] = "duplicate_equivalent_call"
            elif j in near_calls:
                row["reason"] = "no_sequence_equivalent_truth_representation"
            else:
                row["reason"] = "no_truth_in_positive_context" if truth else "no_structural_truth_control"
    window_rows = []
    for k, window in enumerate(negative_windows):
        interval = (window["start"], window["end"])
        positive = any(t["chrom"] == window["chrom"] and
                       any(_distance(interval, v) < negative_buffer or _overlap(v, interval)
                           for v in _intervals(t)) for t in truth)
        fp, uncertain = [], []
        for j, call in enumerate(call_data):
            if call["chrom"] == window["chrom"] and _overlap(_span(call), interval):
                if call_rows[j]["status"] == "fp":
                    fp.append(call_rows[j]["id"])
                elif call_rows[j]["status"] == "unresolved":
                    uncertain.append(call_rows[j]["id"])
        window_rows.append(dict(id=window.get("id", f"window_{k}"), chrom=window["chrom"],
                                start=interval[0], end=interval[1], eligible_negative=not positive,
                                fp_call_ids=fp, unresolved_call_ids=uncertain,
                                false_positive_window=bool(fp), unresolved_window=bool(uncertain)))
    counts = Counter(r["status"] for r in call_rows)
    tp = sum(r["recovered"] for r in event_rows)
    uncertain_events = sum(r["status"] == "unresolved" for r in event_rows)
    negatives = [w for w in window_rows if w["eligible_negative"]]
    fp_windows = sum(w["false_positive_window"] for w in negatives)
    uncertain_only = sum(w["unresolved_window"] and not w["false_positive_window"] for w in negatives)
    ratio = lambda x, n: x/n if n else None
    metrics = dict(n_truth=len(truth), n_calls=len(call_data), tp_events=tp,
                   fn_events=len(truth)-tp-uncertain_events, unresolved_events=uncertain_events,
                   fn_events_strict=len(truth)-tp, tp_calls=counts["tp"], fp_calls=counts["fp"],
                   unresolved_calls=counts["unresolved"], exact_events=sum(r["exact"] for r in event_rows),
                   literal_type_correct_events=sum(r["literal_type_match"] is True for r in event_rows),
                   type_compatible_events=sum(r["type_compatible"] is True for r in event_rows),
                   localized_within_100bp_events=sum(r["localization_within_100bp"] is True for r in event_rows),
                   approximate_events=len(approximate_matches),
                   sensitivity_lower=ratio(tp, len(truth)), sensitivity_upper=ratio(tp+uncertain_events, len(truth)),
                   precision_lower=ratio(counts["tp"], len(call_data)),
                   precision_upper=ratio(counts["tp"]+counts["unresolved"], len(call_data)),
                   n_negative_windows=len(negatives), fp_windows=fp_windows,
                   unresolved_negative_windows=uncertain_only,
                   specificity_lower=ratio(len(negatives)-fp_windows-uncertain_only, len(negatives)),
                   specificity_upper=ratio(len(negatives)-fp_windows, len(negatives)))
    reconstruction = _complete_reconstruction(template_fasta, truth, call_data)
    metrics.update(full_sequence_reconstruction=reconstruction["equivalent"],
                   full_sequence_reconstruction_status=reconstruction["status"])
    return dict(events=event_rows, calls=call_rows, components=components, matches=groups,
                metrics=metrics, negative_windows=window_rows,
                conventions=dict(endpoint="background-SNP-restored structural sequence",
                                 bounds_interpretation="representation/search uncertainty, not statistical confidence intervals",
                                 recall_unit="truth event", precision_unit="call record",
                                 specificity_unit="fixed unaffected reference window",
                                 max_representation_shift=max_representation_shift,
                                 compound_flank=compound_flank, max_component_calls=max_component_calls,
                                 max_component_truth=max_component_truth,
                                 max_compound_comparisons=max_compound_comparisons,
                                 negative_buffer=negative_buffer))
