from __future__ import annotations

from .models import Alignment, AlignmentParameters, SequenceType

MATCH_STATE = 0
DELETION_STATE = 1
INSERTION_STATE = 2

COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def supported_strands(sequence_type: SequenceType) -> tuple[str, ...]:
    if sequence_type in ("dna", "rna"):
        return ("forward", "reverse")
    return ("forward",)


def _path_sort_key(path: tuple[tuple[int | None, int | None], ...]) -> tuple:
    return tuple((ref_index if ref_index is not None else -1, query_index if query_index is not None else -1) for ref_index, query_index in path)


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


def _build_alignment(
    reference: str,
    query: str,
    score: int,
    strand: str,
    end_reference: int,
    end_query: int,
    path: tuple[tuple[int | None, int | None], ...],
) -> Alignment:
    ordered = list(reversed(path))
    reference_columns = [reference[ref_index - 1] if ref_index is not None else "-" for ref_index, _ in ordered]
    query_columns = [query[query_index - 1] if query_index is not None else "-" for _, query_index in ordered]
    reference_coords = [ref_index if ref_index is not None else 0 for ref_index, _ in ordered]
    query_coords = [query_index if query_index is not None else 0 for _, query_index in ordered]
    reference_aligned_positions = [ref_index for ref_index, _ in ordered if ref_index is not None]
    query_aligned_positions = [query_index for _, query_index in ordered if query_index is not None]

    reference_string = "".join(reference_columns)
    query_string = "".join(query_columns)
    operations = [_operation(ref_char, query_char) for ref_char, query_char in zip(reference_columns, query_columns)]

    if strand == "reverse":
        original_query_coords = [len(query) - coord + 1 if coord else 0 for coord in query_coords]
        query_start = len(query) - query_aligned_positions[0] + 1
        query_end = len(query) - query_aligned_positions[-1] + 1
        trace_end_query = len(query) - end_query + 1
    else:
        original_query_coords = query_coords
        query_start, query_end = query_aligned_positions[0], query_aligned_positions[-1]
        trace_end_query = end_query

    return Alignment(
        score=score,
        strand=strand,
        reference_start=reference_aligned_positions[0],
        reference_end=reference_aligned_positions[-1],
        query_start=query_start,
        query_end=query_end,
        aligned_reference=reference_string,
        aligned_query=query_string,
        cigar=_compress_cigar(operations),
        variants=_variants(reference_string, query_string),
        reference_coordinates=tuple(reference_coords),
        query_coordinates=tuple(original_query_coords),
        trace_end_reference=end_reference,
        trace_end_query=trace_end_query,
    )


def _align_one_strand(
    reference: str,
    query: str,
    parameters: AlignmentParameters,
) -> tuple[int, set[tuple[int, int, tuple[tuple[int | None, int | None], ...]]]]:
    parameters.validate()
    scoring = parameters.scoring
    m, n = len(reference), len(query)
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
            pair_score = scoring.match if reference[i - 1] == query[j - 1] else scoring.mismatch
            diagonal = [(scores[state][i - 1][j - 1], state) for state in range(3)]
            valid_diagonal = [(value, state) for value, state in diagonal if value != negative]
            diagonal_best = max((value for value, _ in valid_diagonal), default=negative)
            match_score = diagonal_best + pair_score
            restart = False
            if parameters.mode == "local" and (diagonal_best == negative or match_score < 0):
                match_score = 0
                restart = True
            scores[MATCH_STATE][i][j] = match_score
            if (parameters.mode == "global" or match_score > 0) and not restart:
                predecessors[MATCH_STATE][i][j] = [
                    (state, i - 1, j - 1)
                    for value, state in valid_diagonal if value == diagonal_best
                ]

            deletion_sources = [(scores[state][i - 1][j] + scoring.gap_open, state) for state in range(3)]
            deletion_sources.append((scores[DELETION_STATE][i - 1][j] + scoring.gap_extend, DELETION_STATE))
            deletion_best = max(value for value, _ in deletion_sources)
            if parameters.mode == "global" or deletion_best > 0:
                scores[DELETION_STATE][i][j] = deletion_best
                predecessors[DELETION_STATE][i][j] = [
                    (state, i - 1, j)
                    for value, state in deletion_sources if value == deletion_best
                ]

            insertion_sources = [(scores[state][i][j - 1] + scoring.gap_open, state) for state in range(3)]
            insertion_sources.append((scores[INSERTION_STATE][i][j - 1] + scoring.gap_extend, INSERTION_STATE))
            insertion_best = max(value for value, _ in insertion_sources)
            if parameters.mode == "global" or insertion_best > 0:
                scores[INSERTION_STATE][i][j] = insertion_best
                predecessors[INSERTION_STATE][i][j] = [
                    (state, i, j - 1)
                    for value, state in insertion_sources if value == insertion_best
                ]

    if parameters.mode == "global":
        best_score = max(scores[state][m][n] for state in range(3))
        ends = [(state, m, n) for state in range(3) if scores[state][m][n] == best_score]
    else:
        best_score = max(scores[state][i][j] for state in range(3) for i in range(m + 1) for j in range(n + 1))
        if best_score <= 0:
            return 0, set()
        ends = [
            (state, i, j)
            for state in range(3) for i in range(m + 1) for j in range(n + 1)
            if scores[state][i][j] == best_score
        ]

    candidates: set[tuple[int, int, tuple[tuple[int | None, int | None], ...]]] = set()
    limit_plus_one = parameters.max_alignments + 1

    def traceback(state: int, i: int, j: int, end_reference: int, end_query: int, path: list[tuple[int | None, int | None]]) -> None:
        if len(candidates) >= limit_plus_one:
            return
        edges = predecessors[state][i][j]
        if not edges:
            candidates.add((end_reference, end_query, tuple(path)))
            return
        for previous_state, previous_i, previous_j in edges:
            reference_step = i if previous_i != i else None
            query_step = j if previous_j != j else None
            path.append((reference_step, query_step))
            traceback(previous_state, previous_i, previous_j, end_reference, end_query, path)
            path.pop()

    for state, i, j in ends:
        traceback(state, i, j, i, j, [])

    return int(best_score), candidates


def align_sequences(
    reference: str,
    query: str,
    parameters: AlignmentParameters,
    sequence_type: SequenceType = "dna",
) -> tuple[list[Alignment], bool]:
    parameters.validate()
    searches = []
    for strand in supported_strands(sequence_type):
        if strand == "reverse":
            searches.append((strand, reverse_complement(query)))
        else:
            searches.append((strand, query))

    strand_candidates: list[tuple[str, str, int, set]] = []
    overall_best: int | float = 0 if parameters.mode == "local" else float("-inf")
    for strand, oriented_query in searches:
        best_score, candidates = _align_one_strand(reference, oriented_query, parameters)
        strand_candidates.append((strand, oriented_query, best_score, candidates))
        if parameters.mode == "global" or best_score > 0:
            overall_best = max(overall_best, best_score)

    if parameters.mode == "local" and overall_best <= 0:
        return [], False

    ordered_candidates: list[tuple[str, str, int, int, tuple]] = []
    for strand, oriented_query, best_score, candidates in strand_candidates:
        if best_score != overall_best:
            continue
        for end_reference, end_query, path in candidates:
            ordered_candidates.append((strand, oriented_query, end_reference, end_query, path))
    ordered_candidates.sort(key=lambda item: (item[0], item[2], item[3], _path_sort_key(item[4])))

    limited = len(ordered_candidates) > parameters.max_alignments
    selected = ordered_candidates[:parameters.max_alignments]
    alignments = [
        _build_alignment(
            reference,
            oriented_query,
            overall_best,
            strand,
            end_reference,
            end_query,
            path,
        )
        for strand, oriented_query, end_reference, end_query, path in selected
    ]
    return alignments, limited
