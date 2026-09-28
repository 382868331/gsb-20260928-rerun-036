from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class AlignmentParameters:
    match_score: int = 2
    mismatch_score: int = -1
    gap_penalty: int = 2

    def __post_init__(self) -> None:
        if self.match_score <= 0:
            raise ValueError("match_score must be positive")
        if self.mismatch_score >= 0:
            raise ValueError("mismatch_score must be negative")
        if self.gap_penalty <= 0:
            raise ValueError("gap_penalty must be positive")

    @property
    def gap_score(self) -> int:
        return -self.gap_penalty


@dataclass(frozen=True)
class AlignmentCandidate:
    score: int
    query_start: int
    query_end: int
    reference_start: int
    reference_end: int
    aligned_query: str
    aligned_reference: str
    cigar: str


Move = Literal["diag", "up", "left"]


def _cigar(aligned_query: str, aligned_reference: str) -> str:
    operations: list[str] = []
    current = ""
    count = 0
    for query_char, reference_char in zip(aligned_query, aligned_reference):
        if query_char == "-":
            operation = "D"
        elif reference_char == "-":
            operation = "I"
        elif query_char == reference_char:
            operation = "="
        else:
            operation = "X"
        if operation == current:
            count += 1
        else:
            if count:
                operations.append(f"{count}{current}")
            current = operation
            count = 1
    if count:
        operations.append(f"{count}{current}")
    return "".join(operations)


def smith_waterman(
    query: str,
    reference: str,
    parameters: AlignmentParameters,
    *,
    max_candidates: int = 100,
) -> list[AlignmentCandidate]:
    """Return optimal Smith-Waterman local alignment candidates.

    Tied tracebacks and tied terminal cells are both retained. Coordinates are
    zero-based, half-open intervals. Query gap means CIGAR deletion, reference
    gap means CIGAR insertion.
    """
    if max_candidates <= 0:
        raise ValueError("max_candidates must be positive")

    m = len(query)
    n = len(reference)
    scores = [[0] * (n + 1) for _ in range(m + 1)]
    traces: list[list[set[Move]]] = [
        [set() for _ in range(n + 1)] for _ in range(m + 1)
    ]
    best_score = 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = scores[i - 1][j - 1] + (
                parameters.match_score
                if query[i - 1] == reference[j - 1]
                else parameters.mismatch_score
            )
            up = scores[i - 1][j] + parameters.gap_score
            left = scores[i][j - 1] + parameters.gap_score
            cell_score = max(0, diagonal, up, left)
            scores[i][j] = cell_score
            if cell_score == 0:
                continue
            if diagonal == cell_score:
                traces[i][j].add("diag")
            if up == cell_score:
                traces[i][j].add("up")
            if left == cell_score:
                traces[i][j].add("left")
            best_score = max(best_score, cell_score)

    endpoint_cells = [
        (i, j)
        for i in range(m + 1)
        for j in range(n + 1)
        if scores[i][j] == best_score and traces[i][j]
    ]
    candidates: list[AlignmentCandidate] = []

    def walk(
        i: int,
        j: int,
        endpoint_i: int,
        endpoint_j: int,
        query_parts: list[str],
        reference_parts: list[str],
    ) -> None:
        if len(candidates) >= max_candidates:
            return
        if scores[i][j] == 0:
            aligned_query = "".join(reversed(query_parts))
            aligned_reference = "".join(reversed(reference_parts))
            candidates.append(
                AlignmentCandidate(
                    score=best_score,
                    query_start=i,
                    query_end=endpoint_i,
                    reference_start=j,
                    reference_end=endpoint_j,
                    aligned_query=aligned_query,
                    aligned_reference=aligned_reference,
                    cigar=_cigar(aligned_query, aligned_reference),
                )
            )
            return

        for move in sorted(traces[i][j]):
            if move == "diag":
                query_parts.append(query[i - 1])
                reference_parts.append(reference[j - 1])
                walk(i - 1, j - 1, endpoint_i, endpoint_j, query_parts, reference_parts)
                query_parts.pop()
                reference_parts.pop()
            elif move == "up":
                query_parts.append(query[i - 1])
                reference_parts.append("-")
                walk(i - 1, j, endpoint_i, endpoint_j, query_parts, reference_parts)
                query_parts.pop()
                reference_parts.pop()
            else:
                query_parts.append("-")
                reference_parts.append(reference[j - 1])
                walk(i, j - 1, endpoint_i, endpoint_j, query_parts, reference_parts)
                query_parts.pop()
                reference_parts.pop()
            if len(candidates) >= max_candidates:
                return

    for endpoint_i, endpoint_j in endpoint_cells:
        walk(endpoint_i, endpoint_j, endpoint_i, endpoint_j, [], [])

    candidates.sort(
        key=lambda item: (
            item.reference_start,
            item.query_start,
            item.aligned_query,
            item.aligned_reference,
        )
    )
    return candidates
