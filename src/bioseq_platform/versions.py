from __future__ import annotations

from .features import Feature
from .io import stable_hash
from .storage import ProjectStore
from .validation import normalize_sequence

OPERATION_TYPES = ("rotate", "delete", "insert")


def validate_operation(operation: dict, reference_length: int) -> dict:
    kind = operation.get("type")
    if kind == "rotate":
        offset = int(operation.get("offset", 0))
        if offset < 0:
            raise ValueError("rotate offset must be non-negative")
        return {"type": "rotate", "offset": offset % reference_length}
    if kind == "delete":
        start = int(operation["start"])
        end = int(operation["end"])
        for position in (start, end):
            if not 1 <= position <= reference_length:
                raise ValueError(
                    f"delete interval {start}..{end} outside reference length {reference_length}"
                )
        deleted = _arc_positions(start, end, reference_length)
        if len(deleted) >= reference_length:
            raise ValueError("delete operation must leave at least one base")
        return {"type": "delete", "start": start, "end": end}
    if kind == "insert":
        after = int(operation["after_position"])
        if not 1 <= after <= reference_length:
            raise ValueError(
                f"insert after_position {after} outside reference length {reference_length}"
            )
        bases = normalize_sequence(operation["bases"])
        return {"type": "insert", "after_position": after, "bases": bases}
    raise ValueError(f"unsupported operation type: {kind!r} (expected one of {OPERATION_TYPES})")


def _arc_positions(start: int, end: int, length: int) -> tuple[int, ...]:
    if start <= end:
        return tuple(range(start, end + 1))
    return tuple(range(start, length + 1)) + tuple(range(1, end + 1))


def apply_operation(sequence: str, operation: dict) -> tuple[str, dict[int, int | None], tuple[int, ...]]:
    """Apply one version operation to a circular reference.

    Returns the new sequence, a position map from old 1-based coordinates to
    new 1-based coordinates (None for deleted positions), and the positions
    of inserted bases in new coordinates.
    """
    length = len(sequence)
    kind = operation["type"]
    if kind == "rotate":
        offset = operation["offset"] % length
        position_map = {p: ((p - 1 - offset) % length) + 1 for p in range(1, length + 1)}
        return sequence[offset:] + sequence[:offset], position_map, ()
    if kind == "delete":
        start, end = operation["start"], operation["end"]
        deleted = set(_arc_positions(start, end, length))
        if start <= end:
            survivors = [p for p in range(1, length + 1) if p not in deleted]
        else:
            survivors = list(range(end + 1, start))
        position_map = {p: (index + 1 if p not in deleted else None) for index, p in enumerate(survivors)}
        for p in range(1, length + 1):
            position_map.setdefault(p, None)
        new_sequence = "".join(sequence[p - 1] for p in survivors)
        return new_sequence, position_map, ()
    if kind == "insert":
        after = operation["after_position"]
        bases = operation["bases"]
        count = len(bases)
        position_map = {p: p if p <= after else p + count for p in range(1, length + 1)}
        inserted = tuple(after + step for step in range(1, count + 1))
        return sequence[:after] + bases + sequence[after:], position_map, inserted
    raise ValueError(f"unsupported operation type: {kind!r}")


def initial_units(feature: Feature, reference_length: int) -> list[dict]:
    return [
        {"position": position, "source": "original", "origin": position}
        for position in feature.positions(reference_length)
    ]


def segments_from_units(units: list[dict]) -> list[dict]:
    """Group ordered units into linear segments, splitting at the origin wrap.

    Segment order follows the feature arc (start -> end across the origin when
    wrapping), so a cross-origin feature keeps its two legs in order instead of
    collapsing into a sorted position set.
    """
    runs: list[list[dict]] = []
    current: list[dict] = []
    for unit in units:
        if current and unit["position"] != current[-1]["position"] + 1:
            runs.append(current)
            current = []
        current.append(unit)
    if current:
        runs.append(current)
    return [
        {
            "start": run[0]["position"],
            "end": run[-1]["position"],
            "positions": [unit["position"] for unit in run],
            "sources": [unit["source"] for unit in run],
            "origins": [unit["origin"] for unit in run],
        }
        for run in runs
    ]


def _map_feature_units(
    feature: Feature,
    units: list[dict],
    operation: dict,
    position_map: dict[int, int | None],
    inserted_positions: tuple[int, ...],
    version_id: str,
) -> dict:
    surviving: list[dict] = []
    deleted_units: list[dict] = []
    for unit in units:
        mapped = position_map[unit["position"]]
        if mapped is None:
            deleted_units.append(dict(unit))
        else:
            surviving.append({**unit, "position": mapped})

    impacts: list[dict] = []
    if operation["type"] == "delete" and deleted_units:
        impacts.append({
            "type": "deletion",
            "start": operation["start"],
            "end": operation["end"],
            "location": "covers_feature" if not surviving else "partial",
            "removed_positions": [unit["position"] for unit in deleted_units],
            "removed_origins": [unit["origin"] for unit in deleted_units],
        })

    if operation["type"] == "insert":
        after = operation["after_position"]
        parent_length = len(position_map)
        following = after % parent_length + 1
        parent_positions = {unit["position"] for unit in units}
        before_inside = after in parent_positions
        after_inside = following in parent_positions
        if before_inside and after_inside:
            inserted_units = [
                {
                    "position": position,
                    "source": "inserted",
                    "origin": f"{version_id}:{index}",
                }
                for index, position in enumerate(inserted_positions, start=1)
            ]
            anchor = next(
                index for index, unit in enumerate(surviving) if unit["position"] == after
            )
            surviving[anchor + 1:anchor + 1] = inserted_units
            impacts.append({
                "type": "insertion",
                "after_position": after,
                "inserted_bases": operation["bases"],
                "location": "internal",
                "inserted_origins": [unit["origin"] for unit in inserted_units],
            })
        elif before_inside or after_inside:
            impacts.append({
                "type": "insertion",
                "after_position": after,
                "inserted_bases": operation["bases"],
                "location": "boundary",
            })

    if not surviving:
        status = "deleted"
    elif deleted_units:
        status = "partially_deleted"
    else:
        status = "mapped"

    segments = segments_from_units(surviving)
    wraps = bool(segments) and (
        len(segments) > 1 or segments[0]["start"] > segments[-1]["end"]
    )
    return {
        "feature_id": feature.id,
        "feature_name": feature.name,
        "feature_strand": feature.strand,
        "original": {
            "start": feature.start,
            "end": feature.end,
            "strand": feature.strand,
            "wraps_origin": feature.wraps_origin,
        },
        "status": status,
        "wraps_origin": wraps,
        "segments": segments,
        "units": surviving,
        "deleted_units": deleted_units,
        "impacts": impacts,
    }


def add_version(store: ProjectStore, root_reference_id: str, operation: dict) -> dict:
    chain = store.get_version_chain(root_reference_id)
    versions = chain["versions"]
    parent_reference_id = versions[-1]["reference_id"] if versions else root_reference_id
    parent_version_id = versions[-1]["version_id"] if versions else None
    parent_reference = store.get_reference(parent_reference_id)

    operation = validate_operation(operation, len(parent_reference.sequence))
    version_id = f"ver-{stable_hash([root_reference_id, parent_version_id, operation])[:10]}"
    new_sequence, position_map, inserted_positions = apply_operation(
        parent_reference.sequence, operation
    )
    new_reference = store.add_reference(
        new_sequence,
        sequence_type=parent_reference.sequence_type,
        name=f"{parent_reference.name or parent_reference.id}:{operation['type']}",
        metadata={
            "derived_from": parent_reference_id,
            "version_id": version_id,
            "operation": operation,
        },
    )

    features = store.get_features(root_reference_id)
    if versions:
        previous = {
            mapping["feature_id"]: mapping["units"]
            for mapping in versions[-1]["feature_mappings"]
        }
    else:
        previous = {
            feature.id: initial_units(feature, len(store.get_reference(root_reference_id).sequence))
            for feature in features
        }

    mappings = [
        _map_feature_units(
            feature,
            previous.get(feature.id, []),
            operation,
            position_map,
            inserted_positions,
            version_id,
        )
        for feature in features
    ]

    record = {
        "version_id": version_id,
        "reference_id": new_reference.id,
        "parent_version_id": parent_version_id,
        "operation": operation,
        "feature_mappings": mappings,
    }
    versions.append(record)
    store.save_version_chain(root_reference_id, chain)
    return record


def undo_last_version(store: ProjectStore, root_reference_id: str) -> dict:
    chain = store.get_version_chain(root_reference_id)
    if not chain["versions"]:
        raise ValueError(f"no versions to withdraw for reference {root_reference_id}")
    removed = chain["versions"].pop()
    store.save_version_chain(root_reference_id, chain)
    return removed


def export_mappings(store: ProjectStore, root_reference_id: str, version_index: int | None = None) -> dict:
    chain = store.get_version_chain(root_reference_id)
    versions = chain["versions"]
    if version_index is not None:
        if not 0 <= version_index < len(versions):
            raise IndexError(
                f"version index {version_index} out of range ({len(versions)} versions)"
            )
        versions = versions[: version_index + 1]
    return {
        "root_reference_id": root_reference_id,
        "exported_version_count": len(versions),
        "total_version_count": len(chain["versions"]),
        "versions": versions,
    }


def _compose(step_maps: list[dict[int, int | None]], start: int, length: int) -> dict[int, int]:
    composed: dict[int, int | None] = {p: p for p in range(1, length + 1)}
    for step_map in step_maps[start:]:
        composed = {
            p: step_map[q]
            for p, q in composed.items()
            if q is not None and step_map.get(q) is not None
        }
    return {p: q for p, q in composed.items() if q is not None}


def verify_chain(store: ProjectStore, root_reference_id: str) -> dict:
    """Cross-check coordinates with a direct composed map vs per-version tracking."""
    chain = store.get_version_chain(root_reference_id)
    versions = chain["versions"]
    root_reference = store.get_reference(root_reference_id)
    root_length = len(root_reference.sequence)

    sequences_match = True
    step_maps: list[dict[int, int | None]] = []
    sequence = root_reference.sequence
    for record in versions:
        sequence, position_map, _ = apply_operation(sequence, record["operation"])
        step_maps.append(position_map)
        stored = store.get_reference(record["reference_id"])
        if stored.sequence != sequence or len(stored.sequence) != len(sequence):
            sequences_match = False

    features = store.get_features(root_reference_id)
    feature_reports = []
    for feature in features:
        direct_map = _compose(step_maps, 0, root_length)
        direct_positions = [
            direct_map[p] for p in feature.positions(root_length) if p in direct_map
        ]
        final_mapping = (
            versions[-1]["feature_mappings"] if versions else None
        )
        if final_mapping is None:
            composed_units = initial_units(feature, root_length)
        else:
            composed_units = next(
                m["units"] for m in final_mapping if m["feature_id"] == feature.id
            )
        composed_originals = [u["position"] for u in composed_units if u["source"] == "original"]
        composed_inserted = [u for u in composed_units if u["source"] == "inserted"]

        inserted_match = True
        for unit in composed_inserted:
            version_id, index_text = str(unit["origin"]).rsplit(":", 1)
            version_index = next(
                i for i, record in enumerate(versions) if record["version_id"] == version_id
            )
            operation = versions[version_index]["operation"]
            intro_position = operation["after_position"] + int(index_text)
            expected = _compose_positions(step_maps[version_index + 1:], intro_position)
            if expected != unit["position"]:
                inserted_match = False

        feature_reports.append({
            "feature_id": feature.id,
            "direct_positions": direct_positions,
            "composed_positions": composed_originals,
            "direct_matches_composed": direct_positions == composed_originals,
            "inserted_positions_match": inserted_match,
        })

    ok = sequences_match and all(
        report["direct_matches_composed"] and report["inserted_positions_match"]
        for report in feature_reports
    )
    return {
        "root_reference_id": root_reference_id,
        "version_count": len(versions),
        "sequences_match": sequences_match,
        "features": feature_reports,
        "ok": ok,
    }


def _compose_positions(step_maps: list[dict[int, int | None]], position: int) -> int | None:
    current: int | None = position
    for step_map in step_maps:
        if current is None:
            return None
        current = step_map.get(current)
    return current
