from __future__ import annotations

from dataclasses import dataclass, field

from .models import Alignment


@dataclass(frozen=True)
class Feature:
    id: str
    reference_id: str
    name: str
    start: int
    end: int
    strand: str = "+"
    metadata: dict = field(default_factory=dict)

    @property
    def wraps_origin(self) -> bool:
        return self.start > self.end

    def positions(self, reference_length: int) -> tuple[int, ...]:
        if not 1 <= self.start <= reference_length or not 1 <= self.end <= reference_length:
            raise ValueError(
                f"feature {self.id} interval {self.start}..{self.end} "
                f"outside reference length {reference_length}"
            )
        if self.strand not in ("+", "-"):
            raise ValueError(f"feature {self.id} strand must be '+' or '-'")
        if self.start <= self.end:
            return tuple(range(self.start, self.end + 1))
        return tuple(range(self.start, reference_length + 1)) + tuple(range(1, self.end + 1))

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "reference_id": self.reference_id,
            "name": self.name,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Feature":
        return cls(
            id=data["id"],
            reference_id=data["reference_id"],
            name=data.get("name", ""),
            start=int(data["start"]),
            end=int(data["end"]),
            strand=data.get("strand", "+"),
            metadata=data.get("metadata", {}) if isinstance(data.get("metadata"), dict) else {},
        )


def _deletion_positions(variant: dict, reference_length: int) -> tuple[int, ...]:
    count = len(variant.get("reference_bases", "")) or 1
    start = int(variant["start"])
    return tuple((start - 1 + step) % reference_length + 1 for step in range(count))


def _variant_impacts(feature_positions: frozenset[int], variants: tuple[dict, ...], reference_length: int) -> tuple[dict, ...]:
    impacts: list[dict] = []
    for variant in variants:
        kind = variant.get("type")
        if kind == "substitution":
            position = int(variant["position"])
            if position in feature_positions:
                impacts.append({**variant, "location": "internal"})
        elif kind == "deletion":
            deleted = set(_deletion_positions(variant, reference_length))
            overlap = deleted & feature_positions
            if not overlap:
                continue
            location = "internal" if deleted <= feature_positions else "partial"
            impacts.append({
                **variant,
                "location": location,
                "positions_in_feature": sorted(overlap),
                "covers_feature": feature_positions <= deleted,
            })
        elif kind == "insertion":
            after = int(variant.get("after_position", 0))
            following = after % reference_length + 1 if after >= 1 else 1
            before_inside = after in feature_positions
            after_inside = following in feature_positions
            if before_inside and after_inside:
                impacts.append({**variant, "location": "internal"})
            elif before_inside or after_inside:
                impacts.append({**variant, "location": "boundary"})
    return tuple(impacts)


def _circular_interval(positions: list[int], length: int) -> tuple[int, int, bool]:
    ordered = sorted(set(positions))
    if len(ordered) == 1:
        return ordered[0], ordered[0], False
    widest_gap = -1
    widest_index = 0
    for index, position in enumerate(ordered):
        following = ordered[(index + 1) % len(ordered)]
        gap = (following - position) % length
        if gap > widest_gap:
            widest_gap = gap
            widest_index = index
    start = ordered[(widest_index + 1) % len(ordered)]
    end = ordered[widest_index]
    return start, end, start > end


def project_feature(
    feature: Feature,
    alignment: Alignment,
    reference_length: int,
    query_length: int | None = None,
) -> dict:
    positions = feature.positions(reference_length)
    position_set = frozenset(positions)
    covered: list[int] = []
    mapped_query: list[int] = []
    for reference_coordinate, query_coordinate in zip(
        alignment.reference_coordinates, alignment.query_coordinates
    ):
        if reference_coordinate and reference_coordinate in position_set:
            covered.append(reference_coordinate)
            if query_coordinate:
                mapped_query.append(query_coordinate)

    impacts = _variant_impacts(position_set, alignment.variants, reference_length)
    fully_deleted = bool(covered) and not mapped_query
    if not covered:
        status = "not_covered"
    elif fully_deleted:
        status = "deleted"
    else:
        status = "projected"

    if status == "projected":
        if query_length is None:
            query_length = max(mapped_query)
        query_start, query_end, query_wraps = _circular_interval(mapped_query, query_length)
    else:
        query_start = None
        query_end = None
        query_wraps = False

    return {
        "feature_id": feature.id,
        "feature_name": feature.name,
        "feature_strand": feature.strand,
        "feature_wraps_origin": feature.wraps_origin,
        "status": status,
        "query_start": query_start,
        "query_end": query_end,
        "query_wraps_origin": query_wraps,
        "query_strand": alignment.strand,
        "rotation_offset": alignment.rotation_offset,
        "reference_covered": tuple(sorted(set(covered))),
        "impacts": impacts,
    }


def project_features(
    features: list[Feature],
    alignment: Alignment,
    reference_length: int,
    query_length: int | None = None,
) -> tuple[dict, ...]:
    return tuple(
        project_feature(feature, alignment, reference_length, query_length)
        for feature in features
    )
