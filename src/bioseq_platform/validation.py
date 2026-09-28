from __future__ import annotations

from .models import InvalidRecord, SequenceRecord, SequenceType

DNA_ALPHABET = frozenset("ACGTN")
RNA_ALPHABET = frozenset("ACGUN")
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWY*")


def alphabet_for(sequence_type: SequenceType) -> frozenset[str]:
    if sequence_type == "dna":
        return DNA_ALPHABET
    if sequence_type == "rna":
        return RNA_ALPHABET
    if sequence_type == "protein":
        return PROTEIN_ALPHABET
    raise ValueError(f"unsupported sequence type: {sequence_type}")


def normalize_sequence(raw: object) -> str:
    if not isinstance(raw, str):
        raise ValueError("sequence must be a string")
    sequence = "".join(raw.strip().split()).upper().replace(".", "")
    if not sequence:
        raise ValueError("sequence is empty")
    if "-" in sequence:
        raise ValueError("input sequences cannot contain alignment gaps")
    return sequence


def validate_record(
    record: dict,
    sequence_type: SequenceType,
    default_id: str,
) -> tuple[SequenceRecord | None, InvalidRecord | None]:
    record_id = str(record.get("id") or default_id)
    name = str(record.get("name", ""))
    metadata = record.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {"raw_metadata": metadata}
    raw_sequence = record.get("sequence", "")
    try:
        sequence = normalize_sequence(raw_sequence)
        invalid_bases = sorted(set(sequence) - alphabet_for(sequence_type))
        if invalid_bases:
            raise ValueError(
                f"invalid {sequence_type} character(s): {''.join(invalid_bases)}"
            )
    except ValueError as exc:
        return None, InvalidRecord(
            id=record_id,
            raw_sequence=str(raw_sequence),
            reason=str(exc),
            name=name,
            metadata=metadata,
        )
    return SequenceRecord(
        id=record_id,
        sequence=sequence,
        name=name,
        metadata=metadata,
    ), None


def partition_records(
    records: list[dict],
    sequence_type: SequenceType,
) -> tuple[list[SequenceRecord], list[InvalidRecord]]:
    valid: list[SequenceRecord] = []
    invalid: list[InvalidRecord] = []
    seen: set[str] = set()
    for index, record in enumerate(records, start=1):
        default_id = f"record-{index:04d}"
        parsed, rejected = validate_record(record, sequence_type, default_id)
        target_id = parsed.id if parsed else rejected.id
        if target_id in seen:
            invalid.append(InvalidRecord(
                id=target_id,
                raw_sequence=record.get("sequence", ""),
                reason="duplicate record id",
                name=record.get("name", ""),
                metadata=record.get("metadata", {}) if isinstance(record.get("metadata"), dict) else {},
            ))
            continue
        seen.add(target_id)
        if parsed:
            valid.append(parsed)
        elif rejected:
            invalid.append(rejected)
    return valid, invalid
