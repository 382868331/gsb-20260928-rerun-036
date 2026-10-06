from __future__ import annotations

from itertools import combinations

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


def insert(chain, after_position, bases, candidate_id=None):
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
            {
                "kind": "inserted",
                "edit_version": version_index,
                **({"source_edit_version": candidate_id} if candidate_id is not None else {}),
                "offset": offset,
            }
            for offset in range(len(inserted))
        ]
        + [
            {"kind": "carried", "position": position}
            for position in range(after_position + 1, length + 1)
        ]
    )
    edit = {"kind": INSERT, "after_position": after_position, "bases": inserted}
    if candidate_id is not None:
        edit["candidate_id"] = candidate_id
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


def _inserted_atom(edit_version, offset, source_edit_version=None):
    atom = {"kind": "inserted", "edit_version": edit_version, "offset": offset}
    if source_edit_version is not None:
        atom["source_edit_version"] = source_edit_version
    return atom


def _atom_key(atom):
    if atom["kind"] == "origin":
        return ("origin", atom["position"])
    source = atom.get("source_edit_version", atom["edit_version"])
    return ("inserted", source, atom["offset"])


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
                atom = _inserted_atom(
                    entry["edit_version"],
                    entry["offset"],
                    entry.get("source_edit_version"),
                )
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
                moved[_atom_key(_inserted_atom(
                    index,
                    offset,
                    edit.get("candidate_id"),
                ))] = edit["after_position"] + 1 + offset
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


def _segments(ordered_members, length, strand="+"):
    """Group surviving members into consecutive circular segments.

    ``ordered_members`` is a list of ``(atom, position)`` pairs in feature
    direction.  Each segment carries both the original positions
    (``source_positions``; ``None`` for inserted bases) and the full atom
    descriptors (``sources``), so inserted bases stay traceable after later
    rotations and deletions.
    """
    segments = []
    step = 1 if strand == "+" else length - 1
    for atom, mapped in ordered_members:
        previous = segments[-1]["positions"][-1] if segments else None
        adjacent = previous is not None and mapped == (previous - 1 + step) % length + 1
        crosses_origin = previous is not None and (
            (strand == "+" and previous == length and mapped == 1)
            or (strand == "-" and previous == 1 and mapped == length)
        )
        if adjacent:
            segment = segments[-1]
            if crosses_origin:
                segment["wraps_origin"] = True
        else:
            segment = {
                "positions": [],
                "source_positions": [],
                "sources": [],
                "wraps_origin": False,
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
        if strand == "-":
            segment["wraps_origin"] = segment.get("wraps_origin", False)
        elif segment["start"] > segment["end"]:
            segment["wraps_origin"] = True
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

    ordered_sources = [
        atom["position"]
        for atom in _feature_directed_atoms(feature, len(versions[0]["sequence"]))
    ]
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
                previous_by_position = {position: key for key, position in previous_map.items()}
                after_key = previous_by_position.get(after)
                following_key = previous_by_position.get(after % previous_length + 1)
                predecessor_rank = atom_rank.get(after_key) if feature.strand == "+" else atom_rank.get(following_key)
                successor_rank = atom_rank.get(following_key) if feature.strand == "+" else atom_rank.get(after_key)
                feature_ranks = [atom_rank[key] for key in feature_atoms if key in atom_rank]
                if successor_rank is None:
                    successor_rank = (max(feature_ranks) + 1.0) if feature_ranks else 0.0
                if predecessor_rank is None:
                    predecessor_rank = (min(feature_ranks) - 1.0) if feature_ranks else -1.0
                ordered_offsets = (
                    range(len(edit["bases"]))
                    if feature.strand == "+"
                    else range(len(edit["bases"]) - 1, -1, -1)
                )
                for order_index, offset in enumerate(ordered_offsets):
                    atom = _inserted_atom(index, offset, edit.get("candidate_id"))
                    key = _atom_key(atom)
                    feature_atoms.add(key)
                    feature_positions_map[key] = after + 1 + offset
                    fraction = (order_index + 1) / (len(edit["bases"]) + 1)
                    atom_rank[key] = predecessor_rank + (successor_rank - predecessor_rank) * fraction
                attributed_sources = [
                    _inserted_atom(index, offset, edit.get("candidate_id"))
                    for offset in ordered_offsets
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
            "segments": _segments(ordered_members, length, feature.strand),
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


def _circular_positions(start, end, length):
    if not 1 <= start <= length or not 1 <= end <= length:
        raise ValueError("circular interval %d..%d outside length %d" % (start, end, length))
    if start <= end:
        return tuple(range(start, end + 1))
    return tuple(range(start, length + 1)) + tuple(range(1, end + 1))


def normalize_edit_candidates(candidates, reference_length):
    normalized = []
    seen = set()
    for index, candidate in enumerate(candidates):
        candidate = dict(candidate)
        candidate_id = candidate.get("id") or "candidate-%d" % index
        if candidate_id in seen:
            raise ValueError("duplicate edit candidate id: %s" % candidate_id)
        seen.add(candidate_id)
        candidate["id"] = candidate_id
        kind = candidate.get("kind")
        if kind == ROTATE:
            candidate["offset"] = int(candidate["offset"]) % reference_length
        elif kind == DELETE:
            deleted_atoms = []
            positions = candidate.get("positions")
            if positions is not None:
                deleted_atoms = [
                    {"kind": "origin", "position": int(position)}
                    for position in positions
                ]
            elif "start" not in candidate and candidate.get("deleted_sources"):
                deleted_atoms = []
            else:
                start = int(candidate["start"])
                end = int(candidate.get("end", start))
                deleted_atoms = [
                    {"kind": "origin", "position": position}
                    for position in _circular_positions(start, end, reference_length)
                ]
            for atom in candidate.get("deleted_sources", ()):
                atom = dict(atom)
                if atom.get("kind") not in ("origin", "inserted"):
                    raise ValueError("unsupported deletion source: %s" % atom.get("kind"))
                deleted_atoms.append(atom)
            candidate["deleted_atoms"] = deleted_atoms
        elif kind == INSERT:
            after_source = candidate.get("after_source")
            if after_source is None:
                position = int(candidate.get("after_position", 0))
                if position == 0:
                    after_source = {"kind": "origin", "position": reference_length}
                else:
                    after_source = {"kind": "origin", "position": position}
            candidate["after_source"] = dict(after_source)
            if candidate["after_source"].get("kind") not in ("origin", "inserted"):
                raise ValueError("unsupported insertion anchor source")
            candidate["bases"] = normalize_sequence(candidate.get("bases", ""))
            candidate["atom_sources"] = [
                {"kind": "inserted", "edit_version": candidate_id, "offset": offset}
                for offset in range(len(candidate["bases"]))
            ]
        else:
            raise ValueError("unknown edit candidate kind: %s" % kind)
        normalized.append(candidate)
    return normalized


def _feature_directed_atoms(feature, length):
    positions = feature.positions(length)
    if feature.strand == "-":
        positions = tuple(reversed(positions))
    return [_origin_atom(position) for position in positions]


def _allowed_atoms(allowed_break, feature, length):
    if allowed_break is None:
        allowed_break = {}
    if allowed_break is True or allowed_break == "all":
        return None
    if allowed_break is False or allowed_break == "none":
        return frozenset()
    allowed = set()
    feature_positions = set(feature.positions(length))
    positions = set()
    if "start" in allowed_break and "end" in allowed_break:
        positions.update(_circular_positions(allowed_break["start"], allowed_break["end"], length))
    positions.update(int(position) for position in allowed_break.get("positions", ()))
    allowed.update(
        _atom_key(_origin_atom(position))
        for position in positions
        if position in feature_positions
    )
    allowed.update(_atom_key(atom) for atom in allowed_break.get("sources", ()))
    return frozenset(allowed)


def _candidate_deleted_keys(candidate):
    return {_atom_key(atom) for atom in candidate.get("deleted_atoms", ())}


def _candidate_inserted_keys(candidate):
    return {_atom_key(atom) for atom in candidate.get("atom_sources", ())}


def evaluate_edit_candidate_set(
    candidates,
    selected_ids,
    features,
    reference_length,
    allowed_breaks=None,
    candidate_order=None,
):
    allowed_breaks = allowed_breaks or {}
    order = list(candidate_order or [candidate["id"] for candidate in candidates])

    def ordered_ids(candidate_ids):
        return tuple(sorted(candidate_ids, key=order.index))

    by_id = {candidate["id"]: candidate for candidate in candidates}
    missing = [candidate_id for candidate_id in selected_ids if candidate_id not in by_id]
    if missing:
        raise KeyError("unknown edit candidates: %s" % ", ".join(missing))
    selected = [by_id[candidate_id] for candidate_id in ordered_ids(selected_ids)]
    rotations = [candidate for candidate in selected if candidate["kind"] == ROTATE]
    conflicts = []
    if len({candidate["offset"] for candidate in rotations}) > 1:
        conflicts.append({
            "type": "rotation_mismatch",
            "candidate_ids": ordered_ids(candidate["id"] for candidate in rotations),
            "affected_sources": (),
        })

    inserts = [candidate for candidate in selected if candidate["kind"] == INSERT]
    inserted_by_source = {}
    junction_insertions = {}
    for candidate in inserts:
        for key in _candidate_inserted_keys(candidate):
            inserted_by_source.setdefault(key, []).append(candidate["id"])
        junction_insertions.setdefault(_atom_key(candidate["after_source"]), []).append(candidate["id"])
    for source, candidate_ids in inserted_by_source.items():
        if len(candidate_ids) > 1:
            conflicts.append({
                "type": "duplicate_inserted_source",
                "candidate_ids": ordered_ids(candidate_ids),
                "affected_sources": (_label_atom(source),),
            })
    inserted_keys = set(inserted_by_source)

    deleted_by_source = {}
    deletion_candidates_by_key = {}
    for candidate in selected:
        if candidate["kind"] != DELETE:
            continue
        for key in _candidate_deleted_keys(candidate):
            deleted_by_source.setdefault(key, []).append(candidate["id"])
            deletion_candidates_by_key.setdefault(key, []).append(candidate["id"])
            if key[0] == "inserted" and key not in inserted_keys:
                conflicts.append({
                    "type": "missing_inserted_source",
                "candidate_ids": (candidate["id"],),
                "affected_sources": (_label_atom(key),),
            })

    for candidate in inserts:
        anchor = _atom_key(candidate["after_source"])
        if anchor in deleted_by_source:
            conflicts.append({
                "type": "deleted_insertion_anchor",
                "feature_id": None,
                "candidate_ids": ordered_ids((candidate["id"],) + tuple(deleted_by_source[anchor])),
                "affected_sources": (_label_atom(anchor),),
            })

    for feature in features:
        feature_keys = {_atom_key(atom) for atom in _feature_directed_atoms(feature, reference_length)}
        allowed = _allowed_atoms(allowed_breaks.get(feature.id), feature, reference_length)
        for key in sorted(feature_keys & set(deleted_by_source)):
            if allowed is not None and key not in allowed:
                conflicts.append({
                    "type": "protected_feature_deletion",
                    "feature_id": feature.id,
                    "candidate_ids": ordered_ids(deleted_by_source[key]),
                    "affected_sources": (_label_atom(key),),
                })
    deleted_count = len(deleted_by_source)
    inserted_count = sum(len(candidate["bases"]) for candidate in inserts)
    surviving = reference_length + inserted_count - deleted_count
    if surviving == 0:
        conflicts.append({
            "type": "entire_sequence_deleted",
            "candidate_ids": ordered_ids(selected_ids),
            "affected_sources": (),
        })
    return {
        "feasible": not conflicts,
        "candidate_ids": ordered_ids(selected_ids),
        "rotation_offset": rotations[0]["offset"] if rotations else 0,
        "deleted_sources": tuple(_label_atom(key) for key in sorted(deleted_by_source)),
        "deletion_candidates": {
            str(key): ordered_ids(ids)
            for key, ids in sorted(deletion_candidates_by_key.items())
        },
        "inserted_sources": tuple(_label_atom(key) for key in sorted(inserted_by_source)),
        "shared_junction_insertions": {
            str(key): ordered_ids(ids)
            for key, ids in sorted(junction_insertions.items()) if len(ids) > 1
        },
        "change_cost": {"deleted": deleted_count, "inserted": inserted_count, "total": deleted_count + inserted_count},
        "conflicts": tuple(conflicts),
    }


def plan_edit_candidate_combination(chain, candidates, features, *, selected_feature_ids=None, allowed_breaks=None):
    reference_length = len(chain["versions"][0]["sequence"])
    if selected_feature_ids is not None:
        wanted = set(selected_feature_ids)
        features = [feature for feature in features if feature.id in wanted]
    normalized = normalize_edit_candidates(candidates, reference_length)
    possibilities = []
    conflict_index = {}
    all_ids = sorted(candidate["id"] for candidate in normalized)
    for size in range(len(normalized), -1, -1):
        for subset in combinations(normalized, size):
            ids = tuple(candidate["id"] for candidate in subset)
            result = evaluate_edit_candidate_set(
                normalized, ids, features, reference_length, allowed_breaks,
                all_ids,
            )
            for conflict in result["conflicts"]:
                signature = tuple(sorted((conflict.get("type"), str(conflict.get("affected_sources")))))
                conflict_index[signature] = conflict
            if result["feasible"]:
                possibilities.append(result)
    possibilities.sort(key=lambda result: (
        -len(result["candidate_ids"]),
        result["change_cost"]["total"],
        result["candidate_ids"],
    ))
    best = possibilities[0]
    return {
        "chain_id": chain["id"],
        "candidate_ids": tuple(all_ids),
        "selected_feature_ids": tuple(feature.id for feature in features),
        "candidates": tuple(normalized),
        "best": best,
        "feasible_combinations": tuple(possibilities),
        "conflicts": tuple(conflict_index.values()),
    }


def _delete_atoms(chain, atom_keys):
    atom_positions = composed_atom_map(chain, len(chain["versions"]) - 1)
    positions = sorted(atom_positions[key] for key in atom_keys if key in atom_positions)
    missing = set(atom_keys) - set(atom_positions)
    if missing:
        raise KeyError("cannot delete missing atom sources: %s" % sorted(missing))
    if not positions:
        return chain
    length = len(current_sequence(chain))
    intervals = []
    start = previous = positions[0]
    for position in positions[1:]:
        if position == previous + 1:
            previous = position
        else:
            intervals.append((start, previous))
            start = previous = position
    intervals.append((start, previous))
    merged = []
    for start, end in intervals:
        offset = sum(size for old_start, old_end, size in merged if old_start < start)
        chain = delete(chain, start - offset, end - offset)
        merged.append((start, end, end - start + 1))
    return chain


def apply_confirmed_edit_plan(chain, plan, new_chain_id):
    base = create_chain(chain["reference_id"], chain["versions"][0]["sequence"], new_chain_id)
    best = plan["best"]
    candidates_by_id = {candidate["id"]: candidate for candidate in plan.get("candidates", ())}
    selected = [candidates_by_id[candidate_id] for candidate_id in best["candidate_ids"]]
    deletion_keys = set()
    insertion_candidates = []
    for candidate in selected:
        if candidate["kind"] == DELETE:
            deletion_keys.update(_candidate_deleted_keys(candidate))
        elif candidate["kind"] == INSERT:
            insertion_candidates.append(candidate)

    origin_deletions = {key for key in deletion_keys if key[0] == "origin"}
    if origin_deletions:
        base = _delete_atoms(base, origin_deletions)

    candidate_order = list(best["candidate_ids"])
    pending = {candidate["id"]: candidate for candidate in insertion_candidates}
    while pending:
        progressed = False
        for candidate_id in sorted(pending, key=candidate_order.index, reverse=True):
            candidate = pending[candidate_id]
            anchor = _atom_key(candidate["after_source"])
            atom_positions = composed_atom_map(base, len(base["versions"]) - 1)
            if anchor not in atom_positions:
                continue
            base = insert(
                base,
                atom_positions[anchor],
                candidate["bases"],
                candidate_id=candidate["id"],
            )
            del pending[candidate_id]
            progressed = True
            break
        if not progressed:
            missing = list(pending)
            raise KeyError("insertion anchor sources are unavailable: %s" % ", ".join(missing))

    inserted_deletions = {key for key in deletion_keys if key[0] == "inserted"}
    if inserted_deletions:
        base = _delete_atoms(base, inserted_deletions)

    offset = best.get("rotation_offset", 0)
    if offset:
        base = rotate(base, offset)
    verification = verify_chain(base)
    if not verification["consistent"]:
        raise ValueError("confirmed edit plan produced an inconsistent version chain: %s" % verification["mismatches"])
    return base
