"""Sufficient-only exact contradiction checks for simple BND junction records."""

import hashlib
import re


_BREAKEND = re.compile(r"^([^\[\]]*)([\[\]])([^\[\]]+):(\d+)\2([^\[\]]*)$")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _reverse_complement(sequence):
    return sequence.translate(_COMPLEMENT)[::-1]


def _sha(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def _parse(call):
    alleles = call.get("literal_alt", [])
    if isinstance(alleles, str):
        alleles = [alleles]
    if len(alleles) != 1:
        raise ValueError("BND must have exactly one selected literal ALT")
    match = _BREAKEND.fullmatch(alleles[0])
    if match is None:
        raise ValueError("Unsupported or malformed breakend ALT")
    prefix, bracket, remote_chrom, remote_position, suffix = match.groups()
    if bool(prefix) == bool(suffix):
        raise ValueError("Breakend must have exactly one local anchor")
    anchor = prefix or suffix
    if len(anchor) != 1 or anchor.upper() not in "ACGTN":
        raise ValueError("Inserted or ambiguous BND allele sequence is unsupported")
    local = (str(call["chrom"]), int(call["pos"]), "left" if prefix else "right")
    remote = (remote_chrom, int(remote_position), "right" if bracket == "[" else "left")
    if local[1] < 1 or remote[1] < 1:
        raise ValueError("Breakend position must be positive")
    return dict(local=local, remote=remote, prefix=bool(prefix), bracket=bracket,
                adjacency=tuple(sorted((local, remote))))


def _flank(reference, endpoint, length):
    chrom, pos, side = endpoint
    if side == "left":
        start, end = pos - length, pos
    else:
        start, end = pos - 1, pos - 1 + length
    if start < 0 or end > reference.get_reference_length(chrom):
        raise ValueError("Full diagnostic flank extends beyond a contig boundary")
    sequence = reference.fetch(chrom, start, end).upper()
    if len(sequence) != length or set(sequence) - set("ACGT"):
        raise ValueError("Diagnostic flank contains ambiguous or missing bases")
    return sequence, (chrom, start, end)


def _touches(start, end, flank):
    _, left, right = flank
    return left <= start <= right if start == end else start < right and end > left


def _truth_edits(truth_events):
    by_chrom = {}
    ids = set()
    for event in truth_events:
        if event["id"] in ids:
            raise ValueError("Truth IDs must be unique; namespace chromosome IDs first")
        ids.add(event["id"])
        start, end = int(event["start"]), int(event["end"])
        if start < 0 or end < start:
            raise ValueError("Invalid truth edit coordinates")
        by_chrom.setdefault(event["chrom"], []).append(
            dict(id=event["id"], start=start, end=end, alt=event["alt"].upper()))
    for edits in by_chrom.values():
        edits.sort(key=lambda e: (e["start"], e["end"]))
        cursor, insertions = 0, set()
        for event in edits:
            start, end = event["start"], event["end"]
            if start < cursor or (start == end and start in insertions):
                raise ValueError("Overlapping or unordered truth edits")
            if start == end:
                insertions.add(start)
            cursor = end
    return by_chrom


def _window(reference, chrom, position, shift, flank_length, edits):
    length = reference.get_reference_length(chrom)
    if position > length:
        raise ValueError("Breakend exceeds reference contig length")
    left = max(0, position - 1 - shift - flank_length)
    right = min(length, position + shift + flank_length)
    while True:
        previous = left, right
        for event in edits:
            start, end = event["start"], event["end"]
            if (start <= right and end >= left):
                left = max(0, min(left, start - flank_length))
                right = min(length, max(right, end + flank_length))
        if previous == (left, right):
            break
    source = reference.fetch(chrom, left, right).upper()
    if len(source) != right - left:
        raise ValueError("Incomplete reference window")
    parts, cursor, used = [], left, []
    for event in edits:
        start, end = event["start"], event["end"]
        if left <= start <= end <= right:
            parts.extend((source[cursor-left:start-left], event["alt"]))
            cursor = end
            used.append(event["id"])
    parts.append(source[cursor-left:])
    return "".join(parts), dict(chrom=chrom, start=left, end=right,
                                truth_event_ids=used, background_sha256=_sha(source))


def _scan_complete_chromosome(reference, chrom, edits, junction):
    length = reference.get_reference_length(chrom)
    digest, tail, offset = hashlib.sha256(), "", 0
    forward_at = reverse_at = -1
    ambiguous = False
    reverse = _reverse_complement(junction)

    def chunks():
        cursor = 0
        for event in edits:
            start, end = event["start"], event["end"]
            if end > length:
                raise ValueError("Truth edit exceeds reference contig")
            for left in range(cursor, start, 1_000_000):
                right = min(start, left + 1_000_000)
                sequence = reference.fetch(chrom, left, right).upper()
                if len(sequence) != right - left:
                    raise ValueError("Incomplete reference during global confirmation")
                yield sequence
            for left in range(0, len(event["alt"]), 1_000_000):
                yield event["alt"][left:left+1_000_000]
            cursor = end
        for left in range(cursor, length, 1_000_000):
            right = min(length, left + 1_000_000)
            sequence = reference.fetch(chrom, left, right).upper()
            if len(sequence) != right - left:
                raise ValueError("Incomplete reference during global confirmation")
            yield sequence

    for sequence in chunks():
        digest.update(sequence.encode())
        ambiguous |= bool(set(sequence) - set("ACGT"))
        combined = tail + sequence
        if forward_at < 0:
            at = combined.find(junction)
            if at >= 0:
                forward_at = offset - len(tail) + at
        if reverse_at < 0:
            at = combined.find(reverse)
            if at >= 0:
                reverse_at = offset - len(tail) + at
        offset += len(sequence)
        tail = combined[-(len(junction)-1):]
    return dict(chrom=chrom, reference_length=length, edited_length=offset,
                edited_sha256=digest.hexdigest(), forward_match_offset=forward_at,
                reverse_match_offset=reverse_at, ambiguous_bases_present=ambiguous,
                truth_event_ids=[e["id"] for e in edits])


def validate_false_junctions(reference, truth_events, calls, *,
                             max_representation_shift=1_000_000, flank_length=50):
    if max_representation_shift < 0 or flank_length < 1:
        raise ValueError("Invalid representation shift or flank length")
    truth = _truth_edits(truth_events)
    if len({c["id"] for c in calls}) != len(calls):
        raise ValueError("Call IDs must be unique; namespace query-origin IDs first")
    parsed = {}
    for call in calls:
        if call.get("kind") == "BND":
            try:
                parsed[call["id"]] = _parse(call)
            except (ValueError, KeyError, TypeError):
                pass
    receipts, window_cache, complete_cache = {}, {}, {}
    for call in calls:
        if call.get("kind") != "BND":
            continue
        receipt = dict(call_id=call["id"], status="unassessed", definite_fp=False,
                       full_allele_recovered=False, within_domain_contradiction=False,
                       complete_endpoint_chromosome_absence=False,
                       max_representation_shift=max_representation_shift,
                       flank_length=flank_length,
                       proof_scope="local shift domains followed by complete edited endpoint chromosomes; both orientations")
        receipts[call["id"]] = receipt
        if call["id"] not in parsed:
            receipt["reason"] = "unsupported_breakend_allele"
            continue
        p = parsed[call["id"]]
        receipt["adjacency_id"] = _sha(repr(p["adjacency"]))
        receipt["endpoints"] = [list(e) for e in p["adjacency"]]
        try:
            local, local_span = _flank(reference, p["local"], flank_length)
            remote, remote_span = _flank(reference, p["remote"], flank_length)
        except (ValueError, KeyError) as error:
            receipt["reason"] = str(error)
            continue
        flanks = (local_span, remote_span)
        blockers = []
        for chrom, events in truth.items():
            for event in events:
                if any(chrom == f[0] and _touches(event["start"], event["end"], f) for f in flanks):
                    blockers.append("truth:" + event["id"])
        for other in calls:
            if other["id"] == call["id"]:
                continue
            other_bnd = parsed.get(other["id"])
            if other_bnd:
                if other_bnd["adjacency"] == p["adjacency"]:
                    continue
                spans = [(e[0], e[1]-1, e[1]) for e in (other_bnd["local"], other_bnd["remote"])]
            elif other.get("kind") == "BND":
                blockers.append("unparsed_bnd:" + other["id"])
                continue
            elif other.get("edit") is not None:
                spans = [(other["chrom"], int(other["edit"][0]), int(other["edit"][1]))]
            else:
                blockers.append("unsupported_call:" + other["id"])
                continue
            if any(chrom == f[0] and _touches(start, end, f)
                   for chrom, start, end in spans for f in flanks):
                blockers.append("call:" + other["id"])
        if blockers:
            receipt.update(reason="diagnostic_flank_may_be_modified", blockers=sorted(set(blockers)))
            continue
        if p["prefix"]:
            junction = local + (remote if p["bracket"] == "[" else _reverse_complement(remote))
        else:
            junction = (_reverse_complement(remote) if p["bracket"] == "[" else remote) + local
        reverse = _reverse_complement(junction)
        receipt.update(junction_sequence=junction, junction_sha256=_sha(junction),
                       diagnostic_flanks=[list(f) for f in flanks], windows=[])
        found, ambiguous = False, False
        for endpoint in p["adjacency"]:
            chrom, position, _ = endpoint
            key = (chrom, position)
            if key not in window_cache:
                window_cache[key] = _window(reference, chrom, position, max_representation_shift,
                                             flank_length, truth.get(chrom, []))
            sequence, metadata = window_cache[key]
            forward_at, reverse_at = sequence.find(junction), sequence.find(reverse)
            ambiguous |= bool(set(sequence) - set("ACGT"))
            found |= forward_at >= 0 or reverse_at >= 0
            receipt["windows"].append(dict(metadata, edited_length=len(sequence),
                edited_sha256=_sha(sequence), forward_match_offset=forward_at,
                reverse_match_offset=reverse_at))
        if found:
            receipt.update(status="junction_present_not_full_allele", reason="required_word_present")
        elif ambiguous:
            receipt["reason"] = "ambiguous_bases_in_search_domain"
        else:
            receipt.update(within_domain_contradiction=True, complete_chromosomes=[])
            for chrom in sorted({e[0] for e in p["adjacency"]}):
                key = chrom, junction
                if key not in complete_cache:
                    complete_cache[key] = _scan_complete_chromosome(reference, chrom, truth.get(chrom, []), junction)
                receipt["complete_chromosomes"].append(dict(complete_cache[key]))
            entire = receipt["complete_chromosomes"]
            if any(w["forward_match_offset"] >= 0 or w["reverse_match_offset"] >= 0 for w in entire):
                receipt.update(status="junction_present_outside_local_domain",
                               reason="global_confirmation_found_required_junction")
            elif any(w["ambiguous_bases_present"] for w in entire):
                receipt["reason"] = "ambiguous_bases_in_complete_endpoint_chromosome"
            else:
                receipt.update(status="definite_false_junction", definite_fp=True,
                               complete_endpoint_chromosome_absence=True,
                               reason="required_exact_junction_absent_from_both_complete_edited_endpoint_chromosomes")
    return receipts
