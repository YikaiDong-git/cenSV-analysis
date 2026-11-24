"""Conservative, independent projection of reference calls into native assemblies."""

import math
import re

import edlib


def reverse_complement(sequence):
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _unique_position(sequence, word):
    first = sequence.find(word)
    if first < 0 or sequence.find(word, first + 1) >= 0:
        return None
    return first


def build_correspondence(query_seq, ref_seq, query_offset=0, ref_offset=0, *,
                         kmer_size=31, max_edit_fraction=0.15, anchor_stride=31):
    query_seq, ref_seq = query_seq.upper(), ref_seq.upper()
    if kmer_size < 1 or anchor_stride < 1 or not 0 <= max_edit_fraction <= 1:
        raise ValueError("Invalid anchor or edit-distance parameters")
    result = dict(status="empty_sequence", query_offset=int(query_offset),
                  ref_offset=int(ref_offset), query_length=len(query_seq),
                  ref_length=len(ref_seq), kmer_size=kmer_size,
                  uniqueness_scope="supplied_windows", boundary_map={},
                  anchors=[], gaps=[])
    if not query_seq or not ref_seq:
        return result
    max_distance = math.ceil(len(query_seq) * max_edit_fraction)
    attempts = []
    for reverse, sequence in [(False, query_seq), (True, reverse_complement(query_seq))]:
        alignment = edlib.align(sequence, ref_seq, mode="HW", task="path", k=max_distance)
        if alignment["editDistance"] >= 0:
            attempts.append((alignment["editDistance"], reverse, sequence, alignment))
    if not attempts:
        result.update(status="alignment_failed", max_edit_distance=max_distance)
        return result
    best_distance = min(a[0] for a in attempts)
    winners = [a for a in attempts if a[0] == best_distance]
    result.update(edit_distance=best_distance, edit_fraction=best_distance / len(query_seq))
    if len(winners) != 1:
        result.update(status="ambiguous_orientation")
        return result
    _, reverse, oriented, alignment = winners[0]
    locations = alignment["locations"]
    if len(locations) != 1:
        result.update(status="ambiguous_alignment", n_optimal_locations=len(locations))
        return result
    start, inclusive_end = locations[0]
    if start is None or not alignment.get("cigar"):
        result.update(status="alignment_failed")
        return result
    result.update(reverse=reverse, cigar=alignment["cigar"],
                  aligned_ref_start=ref_offset + start,
                  aligned_ref_end=ref_offset + inclusive_end + 1)
    boundary_map = {}

    def add_boundary(rpos, qpos):
        native = query_offset + (len(query_seq) - qpos if reverse else qpos)
        absolute = ref_offset + rpos
        if absolute not in boundary_map:
            boundary_map[absolute] = native
        elif boundary_map[absolute] != native:
            boundary_map[absolute] = None

    qpos, rpos = 0, start
    add_boundary(rpos, qpos)
    for length, operation in re.findall(r"(\d+)([=XIDM])", alignment["cigar"]):
        length = int(length)
        if operation in "=XM":
            for step in range(1, length + 1):
                add_boundary(rpos + step, qpos + step)
            rpos += length
            qpos += length
        elif operation == "I":
            result["gaps"].append((ref_offset + rpos, ref_offset + rpos, "I"))
            qpos += length
            add_boundary(rpos, qpos)
        elif operation == "D":
            result["gaps"].append((ref_offset + rpos, ref_offset + rpos + length, "D"))
            for step in range(1, length + 1):
                add_boundary(rpos + step, qpos)
            rpos += length
    if qpos != len(query_seq) or rpos != inclusive_end + 1:
        raise ValueError("edlib CIGAR does not consume the declared alignment")
    result["boundary_map"] = boundary_map
    anchor_starts = set(range(0, max(0, len(oriented) - kmer_size + 1), anchor_stride))
    if len(oriented) >= kmer_size:
        anchor_starts.add(len(oriented) - kmer_size)
    for qstart in sorted(anchor_starts):
        word = oriented[qstart:qstart + kmer_size]
        if set(word) - set("ACGT") or _unique_position(oriented, word) != qstart:
            continue
        rstart = _unique_position(ref_seq, word)
        if rstart is None:
            continue
        expected = [query_offset + (len(query_seq) - qstart - j if reverse else qstart + j)
                    for j in range(kmer_size + 1)]
        if any(boundary_map.get(ref_offset + rstart + j) != q
               for j, q in enumerate(expected)):
            continue
        result["anchors"].append(dict(ref_start=ref_offset + rstart,
                                      ref_end=ref_offset + rstart + kmer_size,
                                      query_start=min(expected), query_end=max(expected)))
    result["anchors"].sort(key=lambda a: a["ref_start"])
    result.update(status="ok" if len(result["anchors"]) >= 2 else "insufficient_unique_anchors",
                  n_unique_anchors=len(result["anchors"]))
    return result


def map_query_boundary(correspondence, query_position):
    if correspondence.get("status") != "ok":
        return None
    matches = [r for r, q in correspondence["boundary_map"].items() if q == query_position]
    return matches[0] if len(matches) == 1 else None


def project_call_detailed(correspondence, edit, *, max_anchor_distance=25000):
    if correspondence.get("status") != "ok":
        return dict(status="correspondence_unresolved", edit=None)
    start, end, allele = edit
    if start > end or start < 0:
        raise ValueError("Invalid reference edit interval")
    anchors = correspondence["anchors"]
    left = [a for a in anchors if a["ref_end"] <= start]
    right = [a for a in anchors if a["ref_start"] >= end]
    if not left or not right:
        return dict(status="missing_flanking_anchors", edit=None)
    left, right = left[-1], right[0]
    if start - left["ref_end"] > max_anchor_distance or right["ref_start"] - end > max_anchor_distance:
        return dict(status="distant_flanking_anchors", edit=None)
    for gap_start, gap_end, operation in correspondence["gaps"]:
        overlaps = (left["ref_end"] <= gap_start <= right["ref_start"] if operation == "I"
                    else gap_start < right["ref_start"] and gap_end > left["ref_end"])
        if overlaps:
            return dict(status="background_indel_between_anchors", edit=None)
    mapping = correspondence["boundary_map"]
    qstart, qend = mapping.get(start), mapping.get(end)
    if qstart is None or qend is None:
        return dict(status="unmapped_or_ambiguous_boundary", edit=None)
    if abs(qend - qstart) != end - start:
        return dict(status="noncontiguous_projection", edit=None)
    if correspondence["reverse"]:
        qstart, qend = qend, qstart
        allele = reverse_complement(allele)
    return dict(status="projected", edit=(qstart, qend, allele.upper()),
                left_anchor=left, right_anchor=right)


def project_call(correspondence, edit, *, max_anchor_distance=25000):
    return project_call_detailed(correspondence, edit,
                                 max_anchor_distance=max_anchor_distance)["edit"]


def verify_query_edit(baseline_window, offset, truth_edit, projected_call):
    if projected_call is None:
        return False
    sequence = baseline_window.upper()

    def apply(edit):
        start, end, allele = edit
        start, end = start - offset, end - offset
        if not 0 <= start <= end <= len(sequence):
            return None
        return sequence[:start] + allele.upper() + sequence[end:]

    truth, called = apply(truth_edit), apply(projected_call)
    return truth is not None and called is not None and truth == called
