from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Literal

SequenceType = Literal["dna", "rna", "protein"]
AlignmentMode = Literal["local", "global"]
Topology = Literal["linear", "circular"]


@dataclass(frozen=True)
class Scoring:
    match: int = 2
    mismatch: int = -1
    gap_open: int = -2
    gap_extend: int = -2

    def validate(self) -> None:
        if self.match <= 0:
            raise ValueError("match score must be positive")
        if self.mismatch >= 0:
            raise ValueError("mismatch score must be negative")
        if self.gap_open > 0 or self.gap_extend > 0:
            raise ValueError("gap penalties must be zero or negative")


@dataclass(frozen=True)
class AlignmentParameters:
    mode: AlignmentMode = "local"
    topology: Topology = "linear"
    scoring: Scoring = field(default_factory=Scoring)
    max_alignments: int = 100

    def validate(self) -> None:
        if self.mode not in ("local", "global"):
            raise ValueError("mode must be 'local' or 'global'")
        if self.topology not in ("linear", "circular"):
            raise ValueError("topology must be 'linear' or 'circular'")
        self.scoring.validate()
        if self.max_alignments < 1:
            raise ValueError("max_alignments must be positive")

    def to_dict(self) -> dict:
        return {
            "mode": self.mode,
            "topology": self.topology,
            "scoring": asdict(self.scoring),
            "max_alignments": self.max_alignments,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "AlignmentParameters":
        return cls(
            mode=data["mode"],
            topology=data.get("topology", "linear"),
            scoring=Scoring(**data["scoring"]),
            max_alignments=data["max_alignments"],
        )


@dataclass(frozen=True)
class SequenceRecord:
    id: str
    sequence: str
    name: str = ""
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class Reference:
    id: str
    sequence: str
    name: str = ""
    sequence_type: SequenceType = "dna"
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class InvalidRecord:
    id: str
    raw_sequence: str
    reason: str
    name: str = ""
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class Alignment:
    score: int
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int
    aligned_reference: str
    aligned_query: str
    cigar: str
    strand: str = "forward"
    rotation_offset: int = 0
    wraps_reference: bool = False
    variants: tuple[dict, ...] = ()
    reference_coordinates: tuple[int, ...] = ()
    query_coordinates: tuple[int, ...] = ()
    trace_end_reference: int = 0
    trace_end_query: int = 0


@dataclass(frozen=True)
class SequenceResult:
    record_id: str
    status: str
    score: int
    alignment_count: int
    alignments: tuple[Alignment, ...] = ()
    alignment_limit_reached: bool = False
    selected_alignment_index: int | None = None


@dataclass(frozen=True)
class FunctionalRegion:
    """A directed interval on the original reference record.

    Coordinates are 1-based closed positions on the original (unrotated)
    reference. ``start``/``end`` follow the record start, so a region that
    crosses the record end is represented by ``end < start`` and stays one
    region. ``strand`` records the biological direction of the annotation
    (``forward``/``reverse``), independently of the alignment strand.
    """

    id: str
    start: int
    end: int
    strand: str = "forward"
    name: str = ""
    metadata: dict = field(default_factory=dict)

    def validate(self, reference_length: int) -> None:
        if not self.id:
            raise ValueError("functional region id must be non-empty")
        if self.strand not in ("forward", "reverse"):
            raise ValueError("region strand must be 'forward' or 'reverse'")
        for label, position in (("start", self.start), ("end", self.end)):
            if not 1 <= position <= reference_length:
                raise ValueError(
                    f"region {self.id} {label} {position} is outside the reference "
                    f"(length {reference_length})"
                )

    def member_positions(self, reference_length: int) -> tuple[int, ...]:
        """Reference positions covered by the region in annotation direction."""
        self.validate(reference_length)
        if self.start <= self.end:
            forward = range(self.start, self.end + 1)
        else:
            forward = list(range(self.start, reference_length + 1)) + list(range(1, self.end + 1))
        if self.strand == "reverse":
            forward = list(reversed(list(forward)))
        return tuple(forward)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "FunctionalRegion":
        return cls(
            id=str(data["id"]),
            start=int(data["start"]),
            end=int(data["end"]),
            strand=data.get("strand", "forward"),
            name=str(data.get("name", "")),
            metadata=data.get("metadata", {}) if isinstance(data.get("metadata"), dict) else {},
        )


@dataclass(frozen=True)
class AnnotationSet:
    regions: tuple[FunctionalRegion, ...]
    reference_id: str = ""
    annotation_id: str = ""
