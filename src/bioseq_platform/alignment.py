from __future__ import annotations

from .models import Alignment, AlignmentParameters

MATCH_STATE = 0
DELETION_STATE = 1
INSERTION_STATE = 2


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


def _make_alignment(reference_key: str, query_key: str, score: int) -> Alignment:
    reference_columns = reference_key.split("\u241f")
    query_columns = query_key.split("\u241f")
    steps = [
        (reference_char, query_char, _operation(reference_char, query_char))
        for reference_char, query_char in zip(reference_columns, query_columns)
    ]
    reference_positions = [index + 1 for index, char in enumerate(reference_columns) if char != "-"]
    query_positions = [index + 1 for index, char in enumerate(query_columns) if char != "-"]
    reference_string = "".join(reference_columns)
    query_string = "".join(query_columns)
    return Alignment(
        score=score,
        reference_start=reference_positions[0],
        reference_end=reference_positions[-1],
        query_start=query_positions[0],
        query_end=query_positions[-1],
        aligned_reference=reference_string,
        aligned_query=query_string,
        cigar=_compress_cigar([step[2] for step in steps]),
        variants=_variants(reference_string, query_string),
    )


def align_sequences(reference: str, query: str, parameters: AlignmentParameters) -> tuple[list[Alignment], bool]:
    parameters.validate()
    scoring = parameters.scoring
    m, n = len(reference), len(query)
    negative = float("-inf")
    scores = [[[negative] * (n + 1) for _ in range(m + 1)] for _ in range(3)]
    predecessors = [[[[] for _ in range(n + 1)] for _ in range(m + 1)] for _ in range(3)]

    scores[MATCH_STATE][0][0] = 0
    if parameters.mode == "local":
        for i in range(1, m + 1):
            predecessors[MATCH_STATE][i][0].append((MATCH_STATE, i - 1, 0, "", ""))
        for j in range(1, n + 1):
            predecessors[MATCH_STATE][0][j].append((MATCH_STATE, 0, j - 1, "", ""))
    if parameters.mode == "global":
        for i in range(1, m + 1):
            scores[DELETION_STATE][i][0] = scoring.gap_open + (i - 1) * scoring.gap_extend
            previous_state = DELETION_STATE if i > 1 else MATCH_STATE
            predecessors[DELETION_STATE][i][0].append((previous_state, i - 1, 0, reference[i - 1], "-"))
        for j in range(1, n + 1):
            scores[INSERTION_STATE][0][j] = scoring.gap_open + (j - 1) * scoring.gap_extend
            previous_state = INSERTION_STATE if j > 1 else MATCH_STATE
            predecessors[INSERTION_STATE][0][j].append((previous_state, 0, j - 1, "-", query[j - 1]))

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            pair_score = scoring.match if reference[i - 1] == query[j - 1] else scoring.mismatch
            diagonal = [(scores[state][i - 1][j - 1], state) for state in range(3)]
            diagonal_best = max(value for value, _ in diagonal)
            match_score = diagonal_best + pair_score
            if parameters.mode == "local" and match_score < 0:
                match_score = 0
            scores[MATCH_STATE][i][j] = match_score
            if parameters.mode == "global" or match_score > 0:
                predecessors[MATCH_STATE][i][j] = [
                    (state, i - 1, j - 1, reference[i - 1], query[j - 1])
                    for value, state in diagonal if value == diagonal_best
                ]

            deletion_sources = [(scores[state][i - 1][j] + scoring.gap_open, state) for state in range(3)]
            deletion_sources.append((scores[DELETION_STATE][i - 1][j] + scoring.gap_extend, DELETION_STATE))
            deletion_best = max(value for value, _ in deletion_sources)
            if parameters.mode == "global" or deletion_best > 0:
                scores[DELETION_STATE][i][j] = deletion_best
                predecessors[DELETION_STATE][i][j] = [
                    (state, i - 1, j, reference[i - 1], "-")
                    for value, state in deletion_sources if value == deletion_best
                ]

            insertion_sources = [(scores[state][i][j - 1] + scoring.gap_open, state) for state in range(3)]
            insertion_sources.append((scores[INSERTION_STATE][i][j - 1] + scoring.gap_extend, INSERTION_STATE))
            insertion_best = max(value for value, _ in insertion_sources)
            if parameters.mode == "global" or insertion_best > 0:
                scores[INSERTION_STATE][i][j] = insertion_best
                predecessors[INSERTION_STATE][i][j] = [
                    (state, i, j - 1, "-", query[j - 1])
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
            for state in range(3) for i in range(m + 1) for j in range(n + 1)
            if scores[state][i][j] == best_score
        ]
    best_score = max(scores[state][i][j] for state, i, j in ends)
    ends = [(state, i, j) for state, i, j in ends if scores[state][i][j] == best_score]

    unique: set[tuple[str, str]] = set()
    limit_plus_one = parameters.max_alignments + 1

    def traceback(state: int, i: int, j: int, reference_steps: list[str], query_steps: list[str]) -> None:
        if len(unique) >= limit_plus_one:
            return
        edges = predecessors[state][i][j]
        if not edges:
            unique.add(("\u241f".join(reversed(reference_steps)), "\u241f".join(reversed(query_steps))))
            return
        for previous_state, previous_i, previous_j, reference_char, query_char in edges:
            reference_steps.append(reference_char)
            query_steps.append(query_char)
            traceback(previous_state, previous_i, previous_j, reference_steps, query_steps)
            reference_steps.pop()
            query_steps.pop()

    for state, i, j in ends:
        traceback(state, i, j, [], [])

    limited = len(unique) > parameters.max_alignments
    alignments = [
        _make_alignment(reference_key, query_key, best_score)
        for reference_key, query_key in sorted(unique)[:parameters.max_alignments]
    ]
    return alignments, limited
