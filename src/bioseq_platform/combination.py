from __future__ import annotations

import copy
from itertools import combinations

from .io import stable_hash
from .validation import normalize_sequence
from . import versioning


def _atom_key(atom):
    if atom["kind"] == "origin":
        return ("origin", atom["position"])
    return ("inserted", atom["edit_version"], atom["offset"])


def _atom_label(key):
    if key[0] == "origin":
        return {"kind": "origin", "position": key[1]}
    return {"kind": "inserted", "edit_version": key[1], "offset": key[2]}


def normalize_candidates(candidates):
    normalized = []
    for index, candidate in enumerate(candidates):
        kind = candidate.get("kind")
        candidate_id = candidate.get("id") or "cand-%d" % (index + 1)
        entry = {"id": candidate_id, "kind": kind}
        if kind == versioning.DELETE:
            entry["start"] = int(candidate["start"])
            entry["end"] = int(candidate["end"])
        elif kind == versioning.INSERT:
            entry["after_position"] = int(candidate["after_position"])
            entry["bases"] = normalize_sequence(candidate["bases"])
        elif kind == versioning.ROTATE:
            entry["offset"] = int(candidate["offset"])
        else:
            raise ValueError("unknown candidate kind: %s" % kind)
        normalized.append(entry)
    ids = [entry["id"] for entry in normalized]
    if len(set(ids)) != len(ids):
        raise ValueError("candidate ids must be unique")
    return normalized


def parse_allowed_ranges(ranges):
    parsed = []
    for item in ranges or []:
        start, end = int(item[0]), int(item[1])
        parsed.append((start, end))
    return parsed


def _in_ranges(position, ranges):
    for start, end in ranges:
        if start <= end:
            if start <= position <= end:
                return True
        elif position >= start or position <= end:
            return True
    return False


def _insertion_location(feature_positions, after, length):
    following = after % length + 1 if after >= 1 else 1
    before_inside = after in feature_positions
    after_inside = following in feature_positions
    if before_inside and after_inside:
        return "internal"
    if before_inside or after_inside:
        return "boundary"
    return "outside"


def _kept_feature_views(chain, features, keep_feature_ids, base_version):
    views = []
    for feature in features:
        if keep_feature_ids is not None and feature.id not in keep_feature_ids:
            continue
        trace = versioning.trace_feature(chain, feature, through_version=base_version)
        snapshot = trace["versions"][base_version]
        views.append({
            "feature": feature,
            "status": snapshot["status"],
            "member_keys": [_atom_key(atom) for atom in snapshot["sources"]],
            "positions": set(snapshot["positions"]),
            "segments": snapshot["segments"],
        })
    return views


def _candidate_effect(candidate, pos_to_key, length):
    kind = candidate["kind"]
    if kind == versioning.DELETE:
        start, end = candidate["start"], candidate["end"]
        if not 1 <= start <= end <= length:
            raise ValueError(
                "candidate %s deletion %d..%d outside base version length %d"
                % (candidate["id"], start, end, length)
            )
        positions = set(range(start, end + 1))
        return {
            "deleted_positions": positions,
            "deleted_keys": {pos_to_key[p] for p in positions},
        }
    if kind == versioning.INSERT:
        after = candidate["after_position"]
        if not 0 <= after <= length:
            raise ValueError(
                "candidate %s insertion position %d outside base version length %d"
                % (candidate["id"], after, length)
            )
        return {
            "after_position": after,
            "anchor_key": pos_to_key[after] if after >= 1 else None,
            "following_position": after % length + 1 if after >= 1 else 1,
            "bases": candidate["bases"],
        }
    if kind == versioning.ROTATE:
        shift = candidate["offset"] % length
        return {"shift": shift, "start_key": pos_to_key[shift + 1]}
    raise ValueError("unknown candidate kind: %s" % kind)


def _feature_violations(candidate, effect, views, atom_map, ranges, length):
    violations = []
    kind = candidate["kind"]
    for view in views:
        feature = view["feature"]
        if kind == versioning.DELETE:
            affected = [
                key for key in view["member_keys"] if key in effect["deleted_keys"]
            ]
            if not affected:
                continue
            disallowed = [
                key for key in affected if not _in_ranges(atom_map[key], ranges)
            ]
            if disallowed:
                violations.append({
                    "feature_id": feature.id,
                    "feature_strand": feature.strand,
                    "feature_wraps_origin": feature.wraps_origin,
                    "kind": kind,
                    "affected_sources": [_atom_label(key) for key in affected],
                    "disallowed_sources": [_atom_label(key) for key in disallowed],
                })
        elif kind == versioning.INSERT:
            after = effect["after_position"]
            location = _insertion_location(view["positions"], after, length)
            if location == "outside":
                continue
            involved = {
                position
                for position in (after, effect["following_position"])
                if position in view["positions"]
            }
            disallowed_positions = {
                position for position in involved if not _in_ranges(position, ranges)
            }
            if disallowed_positions:
                affected = [
                    key for key in view["member_keys"] if atom_map[key] in involved
                ]
                violations.append({
                    "feature_id": feature.id,
                    "feature_strand": feature.strand,
                    "feature_wraps_origin": feature.wraps_origin,
                    "kind": kind,
                    "location": location,
                    "affected_sources": [_atom_label(key) for key in affected],
                    "disallowed_positions": sorted(disallowed_positions),
                })
    return violations


def _pair_conflict(first, first_effect, second, second_effect):
    kinds = (first["kind"], second["kind"])
    if kinds == (versioning.DELETE, versioning.INSERT) or kinds == (versioning.INSERT, versioning.DELETE):
        delete_effect = first_effect if first["kind"] == versioning.DELETE else second_effect
        insert_effect = first_effect if first["kind"] == versioning.INSERT else second_effect
        anchor = insert_effect["anchor_key"]
        if anchor is not None and anchor in delete_effect["deleted_keys"]:
            return {
                "candidates": [first["id"], second["id"]],
                "reason": "insert_anchor_deleted",
                "affected_sources": [_atom_label(anchor)],
            }
    if kinds == (versioning.DELETE, versioning.ROTATE) or kinds == (versioning.ROTATE, versioning.DELETE):
        delete_effect = first_effect if first["kind"] == versioning.DELETE else second_effect
        rotate_effect = first_effect if first["kind"] == versioning.ROTATE else second_effect
        start_key = rotate_effect["start_key"]
        if start_key in delete_effect["deleted_keys"]:
            return {
                "candidates": [first["id"], second["id"]],
                "reason": "rotation_start_deleted",
                "affected_sources": [_atom_label(start_key)],
            }
    return None


def study_combinations(chain, features, keep_feature_ids, allowed_ranges, candidates, base_version=None):
    versions = chain["versions"]
    if base_version is None:
        base_version = len(versions) - 1
    if not 0 <= base_version < len(versions):
        raise IndexError("version %d out of range" % base_version)
    length = len(versions[base_version]["sequence"])
    candidates = normalize_candidates(candidates)
    ranges = parse_allowed_ranges(allowed_ranges)

    atom_map = versioning.composed_atom_map(chain, base_version)
    pos_to_key = {position: key for key, position in atom_map.items()}
    views = _kept_feature_views(chain, features, keep_feature_ids, base_version)

    effects = {}
    candidate_reports = []
    for candidate in candidates:
        effect = _candidate_effect(candidate, pos_to_key, length)
        effects[candidate["id"]] = effect
        violations = _feature_violations(candidate, effect, views, atom_map, ranges, length)
        report = {"id": candidate["id"], "kind": candidate["kind"], "violations": violations}
        if candidate["kind"] == versioning.DELETE:
            report["affected_sources"] = [
                _atom_label(pos_to_key[p]) for p in sorted(effect["deleted_positions"])
            ]
        candidate_reports.append(report)

    conflicts = []
    for first, second in combinations(candidates, 2):
        conflict = _pair_conflict(
            first, effects[first["id"]], second, effects[second["id"]]
        )
        if conflict:
            conflicts.append(conflict)

    violating_ids = {
        report["id"] for report in candidate_reports if report["violations"]
    }
    conflicted_pairs = {
        frozenset(conflict["candidates"]) for conflict in conflicts
    }

    def subset_valid(ids):
        if ids & violating_ids:
            return False
        for pair in combinations(ids, 2):
            if frozenset(pair) in conflicted_pairs:
                return False
        deleted = set()
        for candidate_id in ids:
            deleted |= effects[candidate_id].get("deleted_keys", set())
        if len(deleted) >= length:
            return False
        return True

    ranked = []
    for size in range(1, len(candidates) + 1):
        for subset in combinations([c["id"] for c in candidates], size):
            ids = frozenset(subset)
            if not subset_valid(ids):
                continue
            deleted_keys = set()
            inserted_bases = 0
            for candidate_id in ids:
                effect = effects[candidate_id]
                deleted_keys |= effect.get("deleted_keys", set())
                inserted_bases += len(effect.get("bases", ""))
            ranked.append({
                "candidate_ids": sorted(ids),
                "candidate_count": len(ids),
                "deleted_sources": len(deleted_keys),
                "inserted_bases": inserted_bases,
                "change_cost": len(deleted_keys) + inserted_bases,
            })
    ranked.sort(key=lambda item: (-item["candidate_count"], item["change_cost"], item["candidate_ids"]))

    return {
        "chain_id": chain["id"],
        "base_version": base_version,
        "base_length": length,
        "kept_features": [view["feature"].id for view in views],
        "allowed_ranges": [list(item) for item in ranges],
        "candidates": candidate_reports,
        "conflicts": conflicts,
        "combinations": ranked,
        "best": ranked[0] if ranked else None,
    }


def _merged_delete_ranges(positions):
    merged = []
    for position in sorted(positions):
        if merged and position == merged[-1][1] + 1:
            merged[-1][1] = position
        else:
            merged.append([position, position])
    return merged


def confirm_combination(chain, features, keep_feature_ids, allowed_ranges, candidates, selected_ids, new_chain_id=None, base_version=None):
    study = study_combinations(
        chain, features, keep_feature_ids, allowed_ranges, candidates, base_version
    )
    base_version = study["base_version"]
    length = study["base_length"]
    candidates = normalize_candidates(candidates)
    selected_ids = list(selected_ids)
    known = {candidate["id"] for candidate in candidates}
    unknown = [candidate_id for candidate_id in selected_ids if candidate_id not in known]
    if unknown:
        raise ValueError("unknown candidate ids: %s" % ", ".join(unknown))
    selected = [c for c in candidates if c["id"] in set(selected_ids)]

    chosen = frozenset(selected_ids)
    blocking_conflicts = [
        conflict for conflict in study["conflicts"]
        if set(conflict["candidates"]) <= chosen
    ]
    blocking_violations = [
        {"candidate": report["id"], "violations": report["violations"]}
        for report in study["candidates"]
        if report["id"] in chosen and report["violations"]
    ]
    if blocking_conflicts or blocking_violations:
        raise ValueError(
            "selected combination is not compatible: conflicts=%s violations=%s"
            % (blocking_conflicts, blocking_violations)
        )

    atom_map = versioning.composed_atom_map(chain, base_version)
    pos_to_key = {position: key for key, position in atom_map.items()}
    effects = {
        candidate["id"]: _candidate_effect(candidate, pos_to_key, length)
        for candidate in selected
    }

    new_id = new_chain_id or "chain-" + stable_hash([
        chain["id"], base_version, sorted(selected_ids),
    ])[:12]
    new_chain = {
        "id": new_id,
        "reference_id": chain["reference_id"],
        "versions": copy.deepcopy(chain["versions"][: base_version + 1]),
    }

    deleted_positions = set()
    for candidate in selected:
        deleted_positions |= effects[candidate["id"]].get("deleted_positions", set())
    for start, end in _merged_delete_ranges(deleted_positions):
        current = versioning.composed_atom_map(new_chain, len(new_chain["versions"]) - 1)
        keys = {pos_to_key[p] for p in range(start, end + 1)}
        current_positions = sorted(current[key] for key in keys)
        new_chain = versioning.delete(new_chain, current_positions[0], current_positions[-1])

    inserts = [c for c in selected if c["kind"] == versioning.INSERT]
    inserts.sort(key=lambda c: (c["after_position"], c["id"]), reverse=True)
    for candidate in inserts:
        effect = effects[candidate["id"]]
        anchor = effect["anchor_key"]
        if anchor is None:
            after = 0
        else:
            current = versioning.composed_atom_map(new_chain, len(new_chain["versions"]) - 1)
            after = current[anchor]
        new_chain = versioning.insert(new_chain, after, candidate["bases"])

    total_shift = sum(
        effects[c["id"]]["shift"] for c in selected if c["kind"] == versioning.ROTATE
    ) % length
    if total_shift:
        start_key = pos_to_key[total_shift + 1]
        current = versioning.composed_atom_map(new_chain, len(new_chain["versions"]) - 1)
        new_chain = versioning.rotate(new_chain, current[start_key] - 1)

    verification = versioning.verify_chain(new_chain)
    feature_traces = []
    keep = set(keep_feature_ids) if keep_feature_ids is not None else None
    for feature in features:
        if keep is not None and feature.id not in keep:
            continue
        trace = versioning.trace_feature(new_chain, feature)
        final = trace["versions"][-1]
        feature_traces.append({
            "feature_id": feature.id,
            "feature_strand": feature.strand,
            "status": final["status"],
            "positions": final["positions"],
            "sources": final["sources"],
            "segments": final["segments"],
        })

    return {
        "chain": new_chain,
        "new_chain_id": new_id,
        "base_chain_id": chain["id"],
        "base_version": base_version,
        "applied_candidates": sorted(chosen),
        "final_sequence": new_chain["versions"][-1]["sequence"],
        "verification": verification,
        "kept_features": feature_traces,
    }
