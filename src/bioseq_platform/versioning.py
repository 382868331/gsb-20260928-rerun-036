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


def _origin_atom(position):
    return {"kind": "origin", "position": position}


def _inserted_atom(edit_version, offset):
    return {"kind": "inserted", "edit_version": edit_version, "offset": offset}


def _atom_key(atom):
    if atom["kind"] == "origin":
        return ("origin", atom["position"])
    return ("inserted", atom["edit_version"], atom["offset"])


def _step_forward_map(version):
    forward = {}
    for new_position, entry in enumerate(version["provenance"], start=1):
        if entry["kind"] == "carried":
            forward[entry["position"]] = new_position
    return forward


def composed_atom_map(chain, target_version):
    """Position (or None) of every atom ever present up to ``target_version``.

    Built by walking each version's provenance entries, so inserted atoms
    survive subsequent rotations, deletions and further insertions.
    """
    chain = _normalize_provenance(chain)
    versions = chain["versions"]
    if not 0 <= target_version < len(versions):
        raise IndexError("version %d out of range" % target_version)
    current = {}
    initial_length = len(versions[0]["sequence"])
    for position in range(1, initial_length + 1):
        current[_atom_key(_origin_atom(position))] = position
    for version in versions[1 : target_version + 1]:
        by_position = {position: key for key, position in current.items()}
        moved = {}
        for new_position, entry in enumerate(version["provenance"], start=1):
            if entry["kind"] == "carried":
                key = by_position.get(entry["position"])
                if key is not None:
                    moved[key] = new_position
            else:
                atom = _inserted_atom(entry["edit_version"], entry["offset"])
                moved[_atom_key(atom)] = new_position
        current = moved
    return current


def direct_atom_map(chain, target_version):
    """Same mapping as :func:`composed_atom_map`, derived from the edits."""
    if not 0 <= target_version < len(chain["versions"]):
        raise IndexError("version %d out of range" % target_version)
    versions = chain["versions"]
    current = {}
    length = len(versions[0]["sequence"])
    for position in range(1, length + 1):
        current[_atom_key(_origin_atom(position))] = position
    for version in versions[1 : target_version + 1]:
        edit = version["edit"]
        forward, new_length = _edit_transform(edit, length)
        moved = {}
        for key, position in current.items():
            if position in forward:
                moved[key] = forward[position]
        if edit["kind"] == INSERT:
            count = len(edit["bases"])
            index = version["index"]
            for offset in range(count):
                moved[_atom_key(_inserted_atom(index, offset))] = edit["after_position"] + 1 + offset
        current = moved
        length = new_length
    return current


def _label_atom(key):
    if key[0] == "origin":
        return {"kind": "origin", "position": key[1]}
    return {"kind": "inserted", "edit_version": key[1], "offset": key[2]}


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
    for index in range(1, len(chain["versions"])):
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
        direct_atoms = direct_atom_map(chain, index)
        composed_atoms = composed_atom_map(chain, index)
        if direct_atoms != composed_atoms:
            for key in sorted(set(direct_atoms) | set(composed_atoms)):
                if direct_atoms.get(key) != composed_atoms.get(key):
                    mismatches.append({
                        "version": index,
                        "atom": _label_atom(key),
                        "direct": direct_atoms.get(key),
                        "composed": composed_atoms.get(key),
                    })
    return {
        "chain_id": chain["id"],
        "checked_versions": len(chain["versions"]) - 1,
        "consistent": not mismatches,
        "mismatches": mismatches,
    }


def _segments(ordered_members, length):
    """Group surviving members into consecutive circular segments.

    ``ordered_members`` is a list of ``(atom, position)`` pairs in feature
    direction.  Each segment carries both the original positions
    (``source_positions``; ``None`` for inserted bases) and the full atom
    descriptors (``sources``), so inserted bases stay traceable after later
    rotations and deletions.
    """
    segments = []
    for atom, mapped in ordered_members:
        if segments and mapped == segments[-1]["positions"][-1] % length + 1:
            segment = segments[-1]
        else:
            segment = {
                "positions": [],
                "source_positions": [],
                "sources": [],
            }
            segments.append(segment)
        segment["positions"].append(mapped)
        segment["sources"].append(atom)
        segment["source_positions"].append(
            atom["position"] if atom["kind"] == "origin" else None
        )
    for segment in segments:
        segment["start"] = segment["positions"][0]
        segment["end"] = segment["positions"][-1]
        segment["wraps_origin"] = segment["start"] > segment["end"]
    return segments


def _insertion_location(feature_positions, edit, length):
    after = edit["after_position"]
    following = after % length + 1 if after >= 1 else 1
    before_inside = after in feature_positions
    after_inside = following in feature_positions
    if before_inside and after_inside:
        return "internal"
    if before_inside or after_inside:
        return "boundary"
    return "outside"


def _order_members(atom_positions, atom_rank):
    """Return surviving ``(atom, position)`` pairs in feature direction.

    Circular order follows increasing positions; the cut is anchored at the
    surviving member nearest the feature's first original position, which
    survives rotations because the cut is defined by members rather than the
    numeric origin of each version.
    """
    surviving = [(position, key) for key, position in atom_positions.items()]
    if not surviving:
        return []
    surviving.sort(key=lambda item: atom_rank[item[1]])
    return [(_label_atom(key), position) for position, key in surviving]


def trace_feature(chain, feature, through_version=None):
    chain = _normalize_provenance(chain)
    versions = chain["versions"]
    if through_version is None:
        through_version = len(versions) - 1
    if not 0 <= through_version < len(versions):
        raise IndexError("version %d out of range" % through_version)

    ordered_sources = list(feature.positions(len(versions[0]["sequence"])))
    feature_atoms = {
        _atom_key(_origin_atom(source)) for source in ordered_sources
    }
    atom_rank = {
        _atom_key(_origin_atom(source)): rank
        for rank, source in enumerate(ordered_sources)
    }
    deleted_atoms = []
    snapshots = []

    for index, version in enumerate(versions):
        length = len(version["sequence"])
        edit = version["edit"]

        if edit is not None:
            old_positions = composed_atom_map(chain, index - 1)
            position_to_key = {position: key for key, position in old_positions.items()}
            if edit["kind"] == DELETE:
                for position in range(edit["start"], edit["end"] + 1):
                    key = position_to_key.get(position)
                    if key in feature_atoms:
                        feature_atoms.discard(key)
                        if key not in deleted_atoms:
                            deleted_atoms.append(key)
        atom_positions = composed_atom_map(chain, index)
        feature_positions_map = {
            key: atom_positions[key] for key in feature_atoms if key in atom_positions
        }
        ordered_members = _order_members(feature_positions_map, atom_rank)

        insertion_impact = None
        if edit and edit["kind"] == INSERT:
            previous_map = composed_atom_map(chain, index - 1)
            feature_previous = frozenset(
                previous_map[key] for key in feature_atoms if key in previous_map
            )
            previous_length = len(versions[index - 1]["sequence"])
            location = _insertion_location(feature_previous, edit, previous_length)
            attributed_sources = []
            if location == "internal":
                after = edit["after_position"]
                ordered_previous = _order_members(
                    {key: previous_map[key] for key in feature_atoms if key in previous_map},
                    atom_rank,
                )
                predecessor_rank = None
                successor_rank = None
                for order_index, (_, member_position) in enumerate(ordered_previous):
                    if member_position == after:
                        predecessor_rank = atom_rank[
                            _atom_key(ordered_previous[order_index][0])
                        ]
                        if order_index + 1 < len(ordered_previous):
                            successor_rank = atom_rank[
                                _atom_key(ordered_previous[order_index + 1][0])
                            ]
                        break
                feature_ranks = [atom_rank[key] for key in feature_atoms if key in atom_rank]
                if successor_rank is None:
                    successor_rank = (max(feature_ranks) + 1.0) if feature_ranks else 0.0
                if predecessor_rank is None:
                    predecessor_rank = (min(feature_ranks) - 1.0) if feature_ranks else -1.0
                for offset in range(len(edit["bases"])):
                    atom = _inserted_atom(index, offset)
                    key = _atom_key(atom)
                    feature_atoms.add(key)
                    feature_positions_map[key] = after + 1 + offset
                    atom_rank[key] = predecessor_rank + (successor_rank - predecessor_rank) * (offset + 1) / (
                        len(edit["bases"]) + 1
                    )
                attributed_sources = [
                    _inserted_atom(index, offset) for offset in range(len(edit["bases"]))
                ]
                ordered_members = _order_members(feature_positions_map, atom_rank)
            insertion_impact = {
                "edit": edit,
                "location": location,
                "attributed": location == "internal",
                "inserted_sources": tuple(attributed_sources),
            }

        surviving_atoms = feature_atoms & set(atom_positions)
        if not surviving_atoms:
            status = "deleted"
        elif deleted_atoms:
            status = "partially_deleted"
        else:
            status = "projected"
        deleted_in_order = sorted(deleted_atoms, key=lambda key: atom_rank[key])
        deleted_sources = [_label_atom(key) for key in deleted_in_order]
        deleted_source_positions = tuple(
            atom["position"]
            for atom in deleted_sources
            if atom["kind"] == "origin"
        )
        snapshots.append({
            "version": index,
            "sequence_length": length,
            "status": status,
            "positions": tuple(position for _, position in ordered_members),
            "sources": tuple(atom for atom, _ in ordered_members),
            "segments": _segments(ordered_members, length),
            "deleted_sources": tuple(deleted_sources),
            "deleted_source_positions": deleted_source_positions,
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
        "versions": snapshots[: through_version + 1],
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
