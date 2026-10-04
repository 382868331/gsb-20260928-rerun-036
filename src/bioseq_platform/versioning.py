from __future__ import annotations

from .features import Feature
from .io import stable_hash
from .validation import normalize_sequence

ROTATE = "rotate"
DELETE = "delete"
INSERT = "insert"


def create_chain(reference_id, sequence, chain_id=None):
    normalized = normalize_sequence(sequence)
    chain_id = chain_id or "chain-" + stable_hash([reference_id, normalized])[:12]
    return {
        "id": chain_id,
        "reference_id": reference_id,
        "versions": [{
            "index": 0,
            "sequence": normalized,
            "edit": None,
            "provenance": [
                {"kind": "origin", "position": position}
                for position in range(1, len(normalized) + 1)
            ],
        }],
    }


def current_sequence(chain):
    return chain["versions"][-1]["sequence"]


def _append_version(chain, sequence, edit, provenance):
    version = {
        "index": len(chain["versions"]),
        "sequence": sequence,
        "edit": edit,
        "provenance": provenance,
    }
    return {**chain, "versions": chain["versions"] + [version]}


def rotate(chain, offset):
    sequence = current_sequence(chain)
    length = len(sequence)
    shift = offset % length
    rotated = sequence[shift:] + sequence[:shift]
    provenance = [
        {"kind": "carried", "position": (index + shift) % length + 1}
        for index in range(length)
    ]
    edit = {"kind": ROTATE, "offset": shift, "length": length}
    return _append_version(chain, rotated, edit, provenance)


def delete(chain, start, end):
    sequence = current_sequence(chain)
    length = len(sequence)
    if not 1 <= start <= end <= length:
        raise ValueError(
            "deletion interval %d..%d outside version length %d" % (start, end, length)
        )
    shortened = sequence[: start - 1] + sequence[end:]
    if not shortened:
        raise ValueError("deletion would remove the entire reference")
    provenance = [
        {"kind": "carried", "position": position}
        for position in list(range(1, start)) + list(range(end + 1, length + 1))
    ]
    edit = {"kind": DELETE, "start": start, "end": end}
    return _append_version(chain, shortened, edit, provenance)


def insert(chain, after_position, bases):
    sequence = current_sequence(chain)
    length = len(sequence)
    if not 0 <= after_position <= length:
        raise ValueError(
            "insertion position %d outside version length %d" % (after_position, length)
        )
    inserted = normalize_sequence(bases)
    extended = sequence[:after_position] + inserted + sequence[after_position:]
    version_index = len(chain["versions"])
    provenance = (
        [
            {"kind": "carried", "position": position}
            for position in range(1, after_position + 1)
        ]
        + [
            {"kind": "inserted", "edit_version": version_index, "offset": offset}
            for offset in range(len(inserted))
        ]
        + [
            {"kind": "carried", "position": position}
            for position in range(after_position + 1, length + 1)
        ]
    )
    edit = {"kind": INSERT, "after_position": after_position, "bases": inserted}
    return _append_version(chain, extended, edit, provenance)


def undo_last(chain):
    if len(chain["versions"]) <= 1:
        raise ValueError("cannot undo the initial version")
    return {**chain, "versions": chain["versions"][:-1]}


def _normalize_provenance(chain):
    versions = []
    for version in chain["versions"]:
        provenance = []
        for entry in version["provenance"]:
            if entry["kind"] == "origin":
                provenance.append({"kind": "carried", "position": entry["position"]})
            else:
                provenance.append(entry)
        versions.append({**version, "provenance": provenance})
    return {**chain, "versions": versions}


def _step_forward_map(version):
    forward = {}
    for new_position, entry in enumerate(version["provenance"], start=1):
        if entry["kind"] == "carried":
            forward[entry["position"]] = new_position
    return forward


def composed_position_map(chain, target_version):
    if not 0 <= target_version < len(chain["versions"]):
        raise IndexError("version %d out of range" % target_version)
    chain = _normalize_provenance(chain)
    initial = len(chain["versions"][0]["sequence"])
    mapping = {position: position for position in range(1, initial + 1)}
    for version in chain["versions"][1 : target_version + 1]:
        forward = _step_forward_map(version)
        mapping = {
            source: forward[current]
            for source, current in mapping.items()
            if current in forward
        }
    return mapping


def _edit_transform(edit, old_length):
    kind = edit["kind"]
    if kind == ROTATE:
        shift = edit["offset"] % old_length
        return {
            position: (position - 1 - shift) % old_length + 1
            for position in range(1, old_length + 1)
        }, old_length
    if kind == DELETE:
        start, end = edit["start"], edit["end"]
        removed = end - start + 1
        forward = {}
        for position in range(1, old_length + 1):
            if position < start:
                forward[position] = position
            elif position > end:
                forward[position] = position - removed
        return forward, old_length - removed
    if kind == INSERT:
        after = edit["after_position"]
        count = len(edit["bases"])
        return {
            position: position if position <= after else position + count
            for position in range(1, old_length + 1)
        }, old_length + count
    raise ValueError("unknown edit kind: %s" % kind)


def direct_position_map(chain, target_version):
    if not 0 <= target_version < len(chain["versions"]):
        raise IndexError("version %d out of range" % target_version)
    length = len(chain["versions"][0]["sequence"])
    mapping = {position: position for position in range(1, length + 1)}
    for version in chain["versions"][1 : target_version + 1]:
        forward, length = _edit_transform(version["edit"], length)
        mapping = {
            source: forward[current]
            for source, current in mapping.items()
            if current in forward
        }
    return mapping


def verify_chain(chain):
    mismatches = []
    insertion_mismatches = []
    normalized = _normalize_provenance(chain)
    for index in range(1, len(normalized["versions"])):
        direct = direct_position_map(chain, index)
        composed = composed_position_map(chain, index)
        if direct != composed:
            for source in sorted(set(direct) | set(composed)):
                if direct.get(source) != composed.get(source):
                    mismatches.append({
                        "version": index,
                        "source_position": source,
                        "direct": direct.get(source),
                        "composed": composed.get(source),
                    })
    length = len(normalized["versions"][0]["sequence"])
    carried_inserted = []
    for version in normalized["versions"][1:]:
        index = version["index"]
        edit = version["edit"]
        forward_direct, new_length = _edit_transform(edit, length)
        forward_composed = _step_forward_map(version)
        carried_inserted = [
            (label, forward_direct.get(direct), forward_composed.get(composed))
            for label, direct, composed in carried_inserted
        ]
        if edit["kind"] == INSERT:
            after = edit["after_position"]
            for offset in range(len(edit["bases"])):
                expected = {"kind": "inserted", "edit_version": index, "offset": offset}
                actual = version["provenance"][after + offset]
                if actual != expected:
                    insertion_mismatches.append({
                        "version": index,
                        "inserted_offset": offset,
                        "expected": expected,
                        "actual": actual,
                    })
                position = after + 1 + offset
                carried_inserted.append(
                    ("inserted:%d:%d" % (index, offset), position, position)
                )
        length = new_length
        for label, direct, composed in carried_inserted:
            if direct != composed:
                insertion_mismatches.append({
                    "version": index,
                    "source_position": label,
                    "direct": direct,
                    "composed": composed,
                })
    return {
        "chain_id": chain["id"],
        "checked_versions": len(chain["versions"]) - 1,
        "consistent": not mismatches and not insertion_mismatches,
        "mismatches": mismatches,
        "insertion_mismatches": insertion_mismatches,
    }


def _element_provenance(element):
    if element["kind"] == "origin":
        return {"kind": "origin", "position": element["position"]}
    return {
        "kind": "inserted",
        "edit_version": element["edit_version"],
        "offset": element["offset"],
    }


def _segments(ordered_pairs, length):
    segments = []
    for element, mapped in ordered_pairs:
        provenance = _element_provenance(element)
        source = element["position"] if element["kind"] == "origin" else None
        if segments and mapped == segments[-1]["positions"][-1] % length + 1:
            segments[-1]["positions"].append(mapped)
            segments[-1]["source_positions"].append(source)
            segments[-1]["provenance"].append(provenance)
        else:
            segments.append({
                "positions": [mapped],
                "source_positions": [source],
                "provenance": [provenance],
            })
    for segment in segments:
        segment["start"] = segment["positions"][0]
        segment["end"] = segment["positions"][-1]
        segment["wraps_origin"] = segment["start"] > segment["end"]
    return segments


def _insertion_impact(positions, edit, length):
    after = edit["after_position"]
    following = after % length + 1 if after >= 1 else 1
    before_inside = after in positions
    after_inside = following in positions
    if before_inside and after_inside:
        return "internal"
    if before_inside or after_inside:
        return "boundary"
    return "outside"


def trace_feature(chain, feature, through_version=None):
    chain = _normalize_provenance(chain)
    versions = chain["versions"]
    if through_version is None:
        through_version = len(versions) - 1
    if not 0 <= through_version < len(versions):
        raise IndexError("version %d out of range" % through_version)

    ordered_sources = list(feature.positions(len(versions[0]["sequence"])))
    elements = [{"kind": "origin", "position": source} for source in ordered_sources]
    current = list(ordered_sources)
    snapshots = []

    for index, version in enumerate(versions):
        if index > through_version:
            break
        length = len(version["sequence"])
        edit = version["edit"]
        insertion_impact = None
        if index > 0:
            forward = _step_forward_map(version)
            previous = current
            current = [forward.get(position) if position is not None else None for position in previous]
            if edit["kind"] == INSERT:
                present = {position for position in previous if position is not None}
                location = _insertion_impact(present, edit, len(versions[index - 1]["sequence"]))
                insertion_impact = {"edit": edit, "location": location}
                if location == "internal":
                    after = edit["after_position"]
                    count = len(edit["bases"])
                    anchor = next(
                        i for i, position in enumerate(previous) if position == after
                    )
                    new_elements = [
                        {"kind": "inserted", "edit_version": index, "offset": offset}
                        for offset in range(count)
                    ]
                    elements[anchor + 1 : anchor + 1] = new_elements
                    current[anchor + 1 : anchor + 1] = [
                        after + 1 + offset for offset in range(count)
                    ]
        pairs = [
            (element, mapped)
            for element, mapped in zip(elements, current)
            if mapped is not None
        ]
        deleted_sources = [
            element["position"]
            for element, mapped in zip(elements, current)
            if mapped is None and element["kind"] == "origin"
        ]
        deleted_insertions = [
            _element_provenance(element)
            for element, mapped in zip(elements, current)
            if mapped is None and element["kind"] == "inserted"
        ]
        if not deleted_sources and not deleted_insertions:
            status = "projected"
        elif pairs:
            status = "partially_deleted"
        else:
            status = "deleted"
        snapshots.append({
            "version": index,
            "sequence_length": length,
            "status": status,
            "positions": tuple(mapped for _, mapped in pairs),
            "inserted_positions": tuple(
                mapped
                for element, mapped in pairs
                if element["kind"] == "inserted"
            ),
            "segments": _segments(pairs, length),
            "deleted_source_positions": tuple(deleted_sources),
            "deleted_insertions": tuple(deleted_insertions),
            "edit": edit,
            "insertion_impact": insertion_impact,
        })

    return {
        "chain_id": chain["id"],
        "feature_id": feature.id,
        "feature_name": feature.name,
        "feature_strand": feature.strand,
        "feature_start": feature.start,
        "feature_end": feature.end,
        "feature_wraps_origin": feature.wraps_origin,
        "source_positions": tuple(ordered_sources),
        "through_version": through_version,
        "versions": snapshots,
    }


def apply_operation(chain, operation, **kwargs):
    if operation == ROTATE:
        return rotate(chain, kwargs["offset"])
    if operation == DELETE:
        return delete(chain, kwargs["start"], kwargs["end"])
    if operation == INSERT:
        return insert(chain, kwargs["after_position"], kwargs["bases"])
    if operation == "undo":
        return undo_last(chain)
    raise ValueError("unknown chain operation: %s" % operation)
