from __future__ import annotations

import hashlib

from .models import Alignment, AlignmentParameters, SequenceType

MATCH_STATE = 0
DELETION_STATE = 1
INSERTION_STATE = 2

DNA_COMPLEMENTS = str.maketrans("ACGTNacgtn", "TGCANtgcan")
RNA_COMPLEMENTS = str.maketrans("ACGUNacgun", "UGCANugcan")


def _compress_cigar(operations: list[str]) -> str:
    result: list[str] = []
    current = operations[0]
    count = 1
    for operation in operations[1:]:
        if operation == current:
            count += 1
        else:
            result.append(f"{count}{current}")
            current = operation
            count = 1
    result.append(f"{count}{current}")
    return "".join(result)


def _operation(reference_char: str, query_char: str) -> str:
    if reference_char == "-":
        return "I"
    if query_char == "-":
        return "D"
    return "=" if reference_char == query_char else "X"


def _variants(reference: str, query: str) -> tuple[dict, ...]:
    variants: list[dict] = []
    index = 0
    reference_position = 0
    while index < len(reference):
        if reference[index] != "-" and query[index] != "-":
            if reference[index] != query[index]:
                variants.append({
                    "type": "substitution",
                    "position": reference_position + 1,
                    "reference_base": reference[index],
                    "alternative_base": query[index],
                })
            reference_position += 1
            index += 1
        elif query[index] == "-":
            start = reference_position + 1
            bases: list[str] = []
            while index < len(reference) and query[index] == "-":
                bases.append(reference[index])
                reference_position += 1
                index += 1
            variants.append({
                "type": "deletion",
                "start": start,
                "end": reference_position,
                "reference_bases": "".join(bases),
            })
        else:
            inserted: list[str] = []
            after_position = reference_position
            while index < len(reference) and reference[index] == "-":
                inserted.append(query[index])
                index += 1
            variants.append({
                "type": "insertion",
                "after_position": after_position,
                "inserted_bases": "".join(inserted),
            })
    return tuple(variants)


def _reverse_complement(sequence: str, sequence_type: SequenceType) -> str:
    table = RNA_COMPLEMENTS if sequence_type == "rna" else DNA_COMPLEMENTS
    return sequence.translate(table)[::-1]


PathStep = tuple[str, str, int | None, int | None]


def _candidate_id(orientation: str, steps: tuple[PathStep, ...]) -> str:
    identity = (
        orientation,
        tuple((reference_position, query_position) for _, _, reference_position, query_position in steps),
    )
    digest = hashlib.sha256(repr(identity).encode("utf-8")).hexdigest()[:8]
    reference_positions = [position for _, _, position, _ in steps if position is not None]
    query_positions = [position for _, _, _, position in steps if position is not None]
    return (
        f"{orientation}:ref{min(reference_positions)}-{max(reference_positions)}"
        f":query{min(query_positions)}-{max(query_positions)}:{digest}"
    )


def _make_alignment(orientation: str, steps: tuple[PathStep, ...], score: int) -> Alignment:
    reference_string = "".join(reference_char for reference_char, _, _, _ in steps)
    query_string = "".join(query_char for _, query_char, _, _ in steps)
    reference_positions = [position for _, _, position, _ in steps if position is not None]
    query_positions = [position for _, _, _, position in steps if position is not None]
    return Alignment(
        candidate_id=_candidate_id(orientation, steps),
        score=score,
        reference_start=min(reference_positions),
        reference_end=max(reference_positions),
        query_start=min(query_positions),
        query_end=max(query_positions),
        orientation=orientation,
        reference_strand="+",
        query_strand="-" if orientation == "reverse" else "+",
        aligned_reference=reference_string,
        aligned_query=query_string,
        cigar=_compress_cigar([_operation(reference_char, query_char) for reference_char, query_char, _, _ in steps]),
        variants=_variants(reference_string, query_string),
        coordinates=tuple({
            "reference_position": reference_position,
            "query_position": query_position,
        } for _, _, reference_position, query_position in steps),
    )


def _align_one(
    reference: str,
    search_query: str,
    parameters: AlignmentParameters,
    orientation: str,
) -> tuple[list[Alignment], bool]:
    scoring = parameters.scoring
    m, n = len(reference), len(search_query)
    negative = float("-inf")
    scores = [[[negative] * (n + 1) for _ in range(m + 1)] for _ in range(3)]
    predecessors = [[[[] for _ in range(n + 1)] for _ in range(m + 1)] for _ in range(3)]

    scores[MATCH_STATE][0][0] = 0
    if parameters.mode == "local":
        for i in range(1, m + 1):
            scores[MATCH_STATE][i][0] = 0
        for j in range(1, n + 1):
            scores[MATCH_STATE][0][j] = 0
    if parameters.mode == "global":
        for i in range(1, m + 1):
            scores[DELETION_STATE][i][0] = scoring.gap_open + (i - 1) * scoring.gap_extend
            previous_state = DELETION_STATE if i > 1 else MATCH_STATE
            predecessors[DELETION_STATE][i][0].append((previous_state, i - 1, 0))
        for j in range(1, n + 1):
            scores[INSERTION_STATE][0][j] = scoring.gap_open + (j - 1) * scoring.gap_extend
            previous_state = INSERTION_STATE if j > 1 else MATCH_STATE
            predecessors[INSERTION_STATE][0][j].append((previous_state, 0, j - 1))

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            pair_score = scoring.match if reference[i - 1] == search_query[j - 1] else scoring.mismatch
            diagonal = [(scores[state][i - 1][j - 1], state) for state in range(3)]
            diagonal_best = max(value for value, _ in diagonal)
            match_score = diagonal_best + pair_score
            if parameters.mode == "local" and match_score < 0:
                match_score = 0
            scores[MATCH_STATE][i][j] = match_score
            if parameters.mode == "global" or match_score > 0:
                predecessors[MATCH_STATE][i][j] = [
                    (state, i - 1, j - 1)
                    for value, state in diagonal if value == diagonal_best
                ]

            deletion_sources = [
                (scores[state][i - 1][j] + scoring.gap_open, state)
                for state in (MATCH_STATE, INSERTION_STATE)
            ]
            deletion_sources.append((scores[DELETION_STATE][i - 1][j] + scoring.gap_extend, DELETION_STATE))
            deletion_best = max(value for value, _ in deletion_sources)
            if parameters.mode == "global" or deletion_best > 0:
                scores[DELETION_STATE][i][j] = deletion_best
                predecessors[DELETION_STATE][i][j] = [
                    (state, i - 1, j)
                    for value, state in deletion_sources if value == deletion_best
                ]

            insertion_sources = [
                (scores[state][i][j - 1] + scoring.gap_open, state)
                for state in (MATCH_STATE, DELETION_STATE)
            ]
            insertion_sources.append((scores[INSERTION_STATE][i][j - 1] + scoring.gap_extend, INSERTION_STATE))
            insertion_best = max(value for value, _ in insertion_sources)
            if parameters.mode == "global" or insertion_best > 0:
                scores[INSERTION_STATE][i][j] = insertion_best
                predecessors[INSERTION_STATE][i][j] = [
                    (state, i, j - 1)
                    for value, state in insertion_sources if value == insertion_best
                ]

    if parameters.mode == "global":
        ends = [(state, m, n) for state in range(3) if scores[state][m][n] != negative]
    else:
        best_score = max(scores[state][i][j] for state in range(3) for i in range(m + 1) for j in range(n + 1))
        if best_score <= 0:
            return [], False
        ends = [
            (state, i, j)
            for state in range(3)
            for i in range(m + 1)
            for j in range(n + 1)
            if scores[state][i][j] == best_score
        ]

    best_score = max(scores[state][i][j] for state, i, j in ends)
    ends = [(state, i, j) for state, i, j in ends if scores[state][i][j] == best_score]

    unique_paths: dict[tuple, tuple[PathStep, ...]] = {}
    limit_plus_one = parameters.max_alignments + 1

    def traceback(state: int, i: int, j: int, steps: list[PathStep]) -> None:
        if len(unique_paths) >= limit_plus_one:
            return
        edges = predecessors[state][i][j]
        if not edges:
            path = tuple(reversed(steps))
            identity = (
                orientation,
                tuple((reference_position, query_position) for _, _, reference_position, query_position in path),
            )
            unique_paths[identity] = path
            return
        for previous_state, previous_i, previous_j in edges:
            if previous_i < i:
                reference_char = reference[i - 1]
                reference_position: int | None = i
            else:
                reference_char = "-"
                reference_position = None
            if previous_j < j:
                query_char = search_query[j - 1]
                query_position: int | None = n - j + 1 if orientation == "reverse" else j
            else:
                query_char = "-"
                query_position = None
            steps.append((reference_char, query_char, reference_position, query_position))
            traceback(previous_state, previous_i, previous_j, steps)
            steps.pop()

    for state, i, j in ends:
        traceback(state, i, j, [])

    def path_order(path: tuple[PathStep, ...]) -> tuple[int, ...]:
        reference_positions = [position for _, _, position, _ in path if position is not None]
        query_positions = [position for _, _, _, position in path if position is not None]
        return (
            min(reference_positions),
            min(query_positions),
            max(reference_positions),
            max(query_positions),
        )

    limited = len(unique_paths) > parameters.max_alignments
    alignments = [
        _make_alignment(orientation, path, best_score)
        for path in sorted(unique_paths.values(), key=path_order)[:parameters.max_alignments]
    ]
    return alignments, limited


def align_sequences(
    reference: str,
    query: str,
    parameters: AlignmentParameters,
    sequence_type: SequenceType = "dna",
) -> tuple[list[Alignment], bool]:
    parameters.validate()
    search_queries = [("forward", query)]
    if sequence_type in {"dna", "rna"}:
        search_queries.append(("reverse", _reverse_complement(query, sequence_type)))

    candidates: list[Alignment] = []
    truncated = False
    for orientation, search_query in search_queries:
        orientation_candidates, orientation_truncated = _align_one(reference, search_query, parameters, orientation)
        candidates.extend(orientation_candidates)
        truncated = truncated or orientation_truncated

    if not candidates:
        return [], False
    best_score = max(candidate.score for candidate in candidates)
    best_candidates = [candidate for candidate in candidates if candidate.score == best_score]
    best_candidates.sort(key=lambda candidate: (
        candidate.orientation,
        candidate.reference_start,
        candidate.query_start,
        candidate.reference_end,
        candidate.query_end,
        candidate.candidate_id,
    ))
    truncated = truncated or len(best_candidates) > parameters.max_alignments
    return best_candidates[:parameters.max_alignments], truncated
