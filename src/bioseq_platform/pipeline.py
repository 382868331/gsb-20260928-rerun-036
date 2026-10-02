from __future__ import annotations

from dataclasses import asdict
from dataclasses import replace
from pathlib import Path

from .alignment import align_sequences
from .annotations import project_annotation
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


def annotation_from_file(
    store: ProjectStore,
    *,
    annotation_path: str,
    reference_id: str,
    annotation_id: str | None = None,
) -> str:
    from .io import read_annotation_file

    store.initialize()
    regions = read_annotation_file(Path(annotation_path))
    return store.add_annotation(regions, reference_id=reference_id, annotation_id=annotation_id)


def _projection_signature(region_projection: dict) -> dict:
    """Comparable content of one projected region, excluding bookkeeping."""
    return {
        "region_id": region_projection["region_id"],
        "status": region_projection["status"],
        "query_range": region_projection["query_range"],
        "variants": region_projection["variants"],
    }


def build_projections(
    store: ProjectStore,
    *,
    run_id: str,
    record_id: str,
    annotation_id: str,
) -> dict:
    """Project an annotation under every tied candidate of one record.

    Projections for all interpretations are retained (no automatic impact
    minimisation). The previous payload, if any, is embedded so old and new
    projections remain comparable after re-running or re-selecting.
    """
    stored_reference_id, regions = store.get_annotation(annotation_id)
    manifest = store.get_run_manifest(run_id)
    if stored_reference_id != manifest["reference_id"]:
        raise ValueError(
            f"annotation {annotation_id} belongs to reference {stored_reference_id}, "
            f"but run {run_id} uses {manifest['reference_id']}"
        )
    reference_length = len(manifest["reference"]["sequence"])
    results = store.get_run_results(run_id)
    result = next((item for item in results if item.record_id == record_id), None)
    if result is None:
        raise KeyError(f"record {record_id} not found in run {run_id}")

    candidates = []
    for index, alignment in enumerate(result.alignments):
        candidates.append({
            "alignment_index": index,
            "strand": alignment.strand,
            "rotation_offset": alignment.rotation_offset,
            "wraps_reference": alignment.wraps_reference,
            "score": alignment.score,
            "regions": project_annotation(regions, alignment, reference_length),
        })

    previous = store.load_projection(run_id, record_id, annotation_id)
    payload = {
        "run_id": run_id,
        "record_id": record_id,
        "annotation_id": annotation_id,
        "reference_id": stored_reference_id,
        "selected_alignment_index": result.selected_alignment_index,
        "candidates": candidates,
    }
    if previous is not None:
        payload["previous"] = {
            "selected_alignment_index": previous.get("selected_alignment_index"),
            "candidates": previous.get("candidates", []),
        }
    store.save_projection(run_id, record_id, annotation_id, payload)
    return payload


def projection_preview(
    store: ProjectStore,
    *,
    run_id: str,
    record_id: str,
    annotation_id: str,
    target_alignment_index: int,
) -> dict:
    """Preview which projected regions change when switching interpretation.

    The stored selection and projection are not modified; callers use
    ``select_alignment`` followed by ``build_projections`` to adopt the new
    interpretation while keeping the old payload for comparison.
    """
    payload = store.load_projection(run_id, record_id, annotation_id)
    if payload is None:
        raise KeyError(
            f"no projection for annotation {annotation_id} on record {record_id} in run {run_id}"
        )
    current_index = payload["selected_alignment_index"]
    candidates = {item["alignment_index"]: item for item in payload["candidates"]}
    if target_alignment_index not in candidates:
        raise IndexError(
            f"alignment index {target_alignment_index} out of range "
            f"({len(candidates)} projected candidates)"
        )
    if current_index is None:
        raise ValueError("no interpretation selected yet; use select first")
    if current_index not in candidates:
        raise KeyError(f"previously selected candidate {current_index} is no longer available")

    current_regions = {item["region_id"]: item for item in candidates[current_index]["regions"]}
    target_regions = {item["region_id"]: item for item in candidates[target_alignment_index]["regions"]}
    changed: list[dict] = []
    for region_id in sorted(set(current_regions) | set(target_regions)):
        before = current_regions.get(region_id)
        after = target_regions.get(region_id)
        before_signature = _projection_signature(before) if before else None
        after_signature = _projection_signature(after) if after else None
        if before_signature != after_signature:
            changed.append({
                "region_id": region_id,
                "before": before_signature,
                "after": after_signature,
                "needs_recompute": True,
            })
    return {
        "run_id": run_id,
        "record_id": record_id,
        "annotation_id": annotation_id,
        "current_alignment_index": current_index,
        "target_alignment_index": target_alignment_index,
        "changed_region_count": len(changed),
        "changed_regions": changed,
        "recompute_required": bool(changed),
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
