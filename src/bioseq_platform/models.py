from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Literal

SequenceType = Literal["dna", "rna", "protein"]
AlignmentMode = Literal["local", "global"]


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
    scoring: Scoring = field(default_factory=Scoring)
    max_alignments: int = 100

    def validate(self) -> None:
        if self.mode not in ("local", "global"):
            raise ValueError("mode must be 'local' or 'global'")
        self.scoring.validate()
        if self.max_alignments < 1:
            raise ValueError("max_alignments must be positive")

    def to_dict(self) -> dict:
        return {
            "mode": self.mode,
            "scoring": asdict(self.scoring),
            "max_alignments": self.max_alignments,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "AlignmentParameters":
        return cls(
            mode=data["mode"],
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
