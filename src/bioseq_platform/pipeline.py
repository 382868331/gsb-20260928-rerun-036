from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from .alignment import align_sequences
from .io import read_sequence_file, stable_hash
from .models import AlignmentParameters, InvalidRecord, Reference, SequenceRecord, SequenceResult
from .storage import ProjectStore
from .validation import normalize_sequence, partition_records


def run_batch(store: ProjectStore, *, reference: Reference, valid_records: list[SequenceRecord], invalid_records: list[InvalidRecord], parameters: AlignmentParameters, run_id: str | None = None, reuse: bool = True) -> tuple[str, dict, list[SequenceResult]]:
    parameters.validate()
    reference_payload = asdict(reference)
    dataset_payload = {"valid": [asdict(item) for item in valid_records], "invalid": [asdict(item) for item in invalid_records]}
    parameters_payload = parameters.to_dict()
    reference_id = store.add_reference(
        reference.sequence,
        sequence_type=reference.sequence_type,
        name=reference.name,
        metadata=reference.metadata,
        reference_id=reference.id,
    ).id
    dataset_id = store.add_dataset(valid_records, invalid_records)
    parameters_id = store.add_parameters(parameters)
    batch_hash = stable_hash({
        "reference_id": reference_id,
        "reference": reference_payload,
        "dataset_id": dataset_id,
        "dataset": dataset_payload,
        "parameters_id": parameters_id,
        "parameters": parameters_payload,
    })[:16]
    batch_id = run_id or f"run-{batch_hash}"
    if reuse and store.run_exists(batch_id):
        return batch_id, store.get_run_manifest(batch_id), store.get_run_results(batch_id)

    results: list[SequenceResult] = []
    for record in valid_records:
        alignments, limited = align_sequences(reference.sequence, record.sequence, parameters)
        results.append(SequenceResult(
            record_id=record.id,
            status="aligned",
            score=alignments[0].score if alignments else 0,
            alignment_count=len(alignments),
            alignments=tuple(alignments),
            alignment_limit_reached=limited,
        ))
    for record in invalid_records:
        results.append(SequenceResult(record_id=record.id, status=f"invalid: {record.reason}", score=0, alignment_count=0))

    manifest = {
        "run_id": batch_id,
        "reference_id": reference_id,
        "dataset_id": dataset_id,
        "parameters_id": parameters_id,
        "reference": reference_payload,
        "dataset": dataset_payload,
        "parameters": parameters_payload,
        "online_databases": False,
        "valid_record_count": len(valid_records),
        "invalid_record_count": len(invalid_records),
    }
    store.create_run(batch_id, manifest, results)
    return batch_id, manifest, results


def run_from_files(store: ProjectStore, *, reference_path: str, sequences_path: str, sequence_type: str, parameters: AlignmentParameters) -> tuple[str, dict, list[SequenceResult]]:
    reference_records = read_sequence_file(Path(reference_path))
    if len(reference_records) != 1:
        raise ValueError("reference file must contain exactly one FASTA/JSON record")
    raw_reference = reference_records[0]
    reference = Reference(
        id=str(raw_reference.get("id") or f"ref-{stable_hash(normalize_sequence(raw_reference['sequence']))[:12]}"),
        sequence=normalize_sequence(raw_reference["sequence"]),
        name=str(raw_reference.get("name", "")),
        sequence_type=sequence_type,
        metadata=raw_reference.get("metadata", {}) if isinstance(raw_reference.get("metadata"), dict) else {},
    )
    valid, invalid = partition_records(read_sequence_file(Path(sequences_path)), sequence_type)
    return run_batch(store, reference=reference, valid_records=valid, invalid_records=invalid, parameters=parameters)
