from __future__ import annotations

from dataclasses import asdict
from dataclasses import replace
from pathlib import Path

from .alignment import align_sequences
from .features import project_features
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
        alignments, limited = align_sequences(
            reference.sequence,
            record.sequence,
            parameters,
            sequence_type=reference.sequence_type,
        )
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


def select_alignment(store: ProjectStore, run_id: str, record_id: str, alignment_index: int) -> list[SequenceResult]:
    results = store.get_run_results(run_id)
    for index, result in enumerate(results):
        if result.record_id != record_id:
            continue
        if not 0 <= alignment_index < len(result.alignments):
            raise IndexError(
                f"alignment index {alignment_index} out of range for record {record_id} "
                f"({len(result.alignments)} candidates)"
            )
        results[index] = replace(result, selected_alignment_index=alignment_index)
        store.save_run_results(run_id, results)
        return results
    raise KeyError(f"record {record_id} not found in run {run_id}")


def _run_context(store: ProjectStore, run_id: str, record_id: str):
    manifest = store.get_run_manifest(run_id)
    reference_id = manifest["reference_id"]
    reference_length = len(manifest["reference"]["sequence"])
    query_length = 0
    for item in manifest["dataset"]["valid"]:
        if item["id"] == record_id:
            query_length = len(item["sequence"])
            break
    features = store.get_features(reference_id)
    results = store.get_run_results(run_id)
    for result in results:
        if result.record_id == record_id:
            return features, reference_length, query_length or None, result
    raise KeyError(f"record {record_id} not found in run {run_id}")


def project_run_record(store: ProjectStore, run_id: str, record_id: str, alignment_index: int | None = None) -> dict:
    features, reference_length, query_length, result = _run_context(store, run_id, record_id)
    indices = (
        [alignment_index]
        if alignment_index is not None
        else list(range(len(result.alignments)))
    )
    projections = []
    for index in indices:
        if not 0 <= index < len(result.alignments):
            raise IndexError(
                f"alignment index {index} out of range for record {record_id} "
                f"({len(result.alignments)} candidates)"
            )
        alignment = result.alignments[index]
        projections.append({
            "alignment_index": index,
            "selected": result.selected_alignment_index == index,
            "strand": alignment.strand,
            "rotation_offset": alignment.rotation_offset,
            "score": alignment.score,
            "features": project_features(features, alignment, reference_length, query_length),
        })
    return {
        "run_id": run_id,
        "record_id": record_id,
        "selected_alignment_index": result.selected_alignment_index,
        "feature_count": len(features),
        "projections": projections,
    }


def selection_projection_diff(store: ProjectStore, run_id: str, record_id: str, new_index: int) -> dict:
    features, reference_length, query_length, result = _run_context(store, run_id, record_id)
    old_index = result.selected_alignment_index
    if not 0 <= new_index < len(result.alignments):
        raise IndexError(
            f"alignment index {new_index} out of range for record {record_id} "
            f"({len(result.alignments)} candidates)"
        )
    changed: list[str] = []
    if old_index is not None and features:
        old_projection = project_features(features, result.alignments[old_index], reference_length, query_length)
        new_projection = project_features(features, result.alignments[new_index], reference_length, query_length)
        changed = [
            old["feature_id"]
            for old, new in zip(old_projection, new_projection)
            if old != new
        ]
    elif features:
        changed = [feature.id for feature in features]
    return {
        "run_id": run_id,
        "record_id": record_id,
        "previous_alignment_index": old_index,
        "new_alignment_index": new_index,
        "features_to_recompute": changed,
    }


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
