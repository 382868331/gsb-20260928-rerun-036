from __future__ import annotations

from .models import Alignment, FunctionalRegion

# Insertions are anchored to the reference base immediately before the gap
# columns (the variant's "after_position"). A boundary insertion is
# therefore attributed, exactly once, to the region whose directed end is that
# anchor base; it is never re-counted as a full interior event in the
# neighbouring region that starts at the following base.


def region_members(region: FunctionalRegion, reference_length: int) -> list[int]:
    """Reference positions covered by the region, in annotation direction."""
    if region.start <= region.end:
        forward = list(range(region.start, region.end + 1))
    else:
        forward = list(range(region.start, reference_length + 1)) + list(range(1, region.end + 1))
    if region.strand == "reverse":
        forward.reverse()
    return forward


def _to_segments(directed_positions: list[int], length: int, direction: str) -> list[dict]:
    """Collapse consecutive directed coordinates into 1-based closed ranges.

    Forward ranges report [start, end] with start <= end, except a segment
    crossing the record end which is flagged with "wraps". Reverse ranges
    report [high, low] so the annotation/alignment reading order is retained.
    """
    if not directed_positions:
        return []
    groups: list[list[int]] = [[directed_positions[0]]]
    for position in directed_positions[1:]:
        previous = groups[-1][-1]
        if direction == "forward":
            consecutive = position == previous + 1 or (previous == length and position == 1)
        else:
            consecutive = position == previous - 1 or (previous == 1 and position == length)
        if consecutive:
            groups[-1].append(position)
        else:
            groups.append([position])
    return [
        {
            "start": group[0],
            "end": group[-1],
            "wraps": direction == "forward" and group[0] > group[-1],
        }
        for group in groups
    ]


def region_index(region: FunctionalRegion, reference_length: int, position: int) -> int:
    """1-based offset of a reference position inside the directed region."""
    return region_members(region, reference_length).index(position) + 1


def last_aligned_reference(ref_coords: list[int], column: int) -> int | None:
    index = column - 1
    while index >= 0 and not ref_coords[index]:
        index -= 1
    if index < 0:
        return None
    return ref_coords[index]


def is_region_end_anchor(region: FunctionalRegion, reference_length: int, anchor: int | None) -> bool:
    if anchor is None:
        return False
    return region_members(region, reference_length)[-1] == anchor


def project_region(region: FunctionalRegion, alignment: Alignment, reference_length: int) -> dict:
    """Project one directed reference region onto a single alignment.

    All returned coordinates stay on the original input sequences: reference
    coordinates are never rotated away and reverse-complement alignments use
    original (un-reversed) query coordinates with reversed ordering.
    """
    region.validate(reference_length)
    ordered_members = region_members(region, reference_length)
    members = set(ordered_members)
    ref_coords = list(alignment.reference_coordinates)
    query_coords = list(alignment.query_coordinates)
    query_length = max((coord for coord in query_coords if coord), default=0)
    column_count = len(ref_coords)

    deleted: list[int] = []
    deletion_runs: list[dict] = []
    substitutions: list[dict] = []
    insertions: list[dict] = []

    column = 0
    while column < column_count:
        ref_position = ref_coords[column]
        query_position = query_coords[column]
        in_region = bool(ref_position) and ref_position in members

        if in_region and not query_position:
            run_ref: list[int] = []
            while (
                column < column_count
                and ref_coords[column] in members
                and not query_coords[column]
            ):
                run_ref.append(ref_coords[column])
                column += 1
            deleted.extend(run_ref)
            deletion_runs.append({
                "type": "deletion",
                "reference_start": run_ref[0],
                "reference_end": run_ref[-1],
                "region_offset": min(
                    region_index(region, reference_length, position) for position in run_ref
                ),
                "region_length": len(members),
                "span": len(run_ref),
                "removes_region": set(run_ref) == members,
            })
            continue

        if in_region and query_position:
            reference_base = alignment.aligned_reference[column]
            query_base = alignment.aligned_query[column]
            if reference_base != query_base:
                substitutions.append({
                    "type": "substitution",
                    "reference_position": ref_position,
                    "query_position": query_position,
                    "region_offset": region_index(region, reference_length, ref_position),
                    "region_length": len(members),
                    "reference_base": reference_base,
                    "alternative_base": query_base,
                })
            column += 1
            continue

        if not ref_position:
            anchor = last_aligned_reference(ref_coords, column)
            if anchor in members:
                inserted_query: list[int] = []
                inserted_bases: list[str] = []
                while column < column_count and not ref_coords[column]:
                    if query_coords[column]:
                        inserted_query.append(query_coords[column])
                    inserted_bases.append(alignment.aligned_query[column])
                    column += 1
                boundary = is_region_end_anchor(region, reference_length, anchor)
                insertions.append({
                    "type": "insertion",
                    "after_reference_position": anchor,
                    "region_offset": region_index(region, reference_length, anchor),
                    "region_length": len(members),
                    "query_positions": inserted_query,
                    "inserted_bases": "".join(inserted_bases),
                    "location": "end_boundary" if boundary else "interior",
                })
            else:
                column += 1
            continue

        column += 1

    deleted_set = set(deleted)
    query_by_reference = {
        ref_coords[index]: query_coords[index]
        for index in range(column_count)
        if ref_coords[index] in members and query_coords[index]
    }
    surviving = [
        position
        for position in ordered_members
        if position not in deleted_set and position in query_by_reference
    ]

    if not surviving and deleted_set == members:
        status = "deleted"
        query_range = None
    elif not surviving:
        status = "unaligned"
        query_range = None
    else:
        status = "interrupted" if deleted else "covered"
        query_directed = [query_by_reference[position] for position in surviving]
        if alignment.strand == "reverse":
            query_directed.reverse()
        query_range = {
            "start": query_directed[0],
            "end": query_directed[-1],
            "strand": alignment.strand,
            "segments": _to_segments(
                query_directed,
                query_length,
                "reverse" if alignment.strand == "reverse" else "forward",
            ),
        }

    kept = [position for position in ordered_members if position not in deleted_set]
    deleted_directed = [position for position in ordered_members if position in deleted_set]
    return {
        "region_id": region.id,
        "region_name": region.name,
        "region_strand": region.strand,
        "region_start": region.start,
        "region_end": region.end,
        "region_wraps": region.end < region.start,
        "region_length": len(members),
        "status": status,
        "query_range": query_range,
        "projected_reference_segments": _to_segments(
            kept if region.strand == "forward" else list(reversed(kept)),
            reference_length,
            region.strand,
        ),
        "deleted_reference_segments": _to_segments(
            deleted_directed if region.strand == "forward" else list(reversed(deleted_directed)),
            reference_length,
            region.strand,
        ),
        "variants": {
            "substitutions": substitutions,
            "deletions": deletion_runs,
            "insertions": insertions,
        },
    }


def sort_key(region: FunctionalRegion, reference_length: int) -> tuple:
    return (
        region.start,
        0 if region.end >= region.start else 1,
        -len(region_members(region, reference_length)),
        region.id,
    )


def ordered_regions(regions, reference_length: int) -> list[FunctionalRegion]:
    return sorted(regions, key=lambda region: sort_key(region, reference_length))


def project_annotation(regions, alignment: Alignment, reference_length: int) -> list[dict]:
    return [
        project_region(region, alignment, reference_length)
        for region in ordered_regions(regions, reference_length)
    ]
