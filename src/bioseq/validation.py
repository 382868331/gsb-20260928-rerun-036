from __future__ import annotations

ALPHABETS = {
    "dna": frozenset("ACGT"),
    "rna": frozenset("ACGU"),
    "protein": frozenset("ACDEFGHIKLMNPQRSTVWY"),
}


def normalize_sequence(sequence: str) -> str:
    return "".join(sequence.strip().upper().split())


def validate_sequence(sequence: str, sequence_type: str) -> list[str]:
    try:
        alphabet = ALPHABETS[sequence_type]
    except KeyError as exc:
        raise ValueError(f"unsupported sequence_type: {sequence_type}") from exc
    normalized = normalize_sequence(sequence)
    if not normalized:
        return ["sequence is empty"]
    invalid = sorted({char for char in normalized if char not in alphabet})
    if invalid:
        return ["invalid characters: " + " ".join(invalid)]
    return []
