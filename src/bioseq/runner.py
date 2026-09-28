from __future__ import annotations

from dataclasses import asdict
from typing import Any

from .alignment import AlignmentParameters, smith_waterman
from .storage import PlatformStore, read_json, stable_hash, write_json
from .validation import normalize_sequence, validate_sequence


def _batch_id(
    dataset_id: str,
    reference_id: str,
    parameters_id: str,
    max_candidates: int,
) -> str:
    return "batch_" + stable_hash(
        {
            "dataset_id": dataset_id,
            "reference_id": reference_id,
            "parameters_id": parameters_id,
            "max_candidates": max_candidates,
            "algorithm": "smith-waterman-linear-gap",
            "schema_version": 1,
        }
    )[:16]


def run_batch(
    store: PlatformStore,
    dataset_id: str,
    reference_id: str,
    parameters: AlignmentParameters,
    *,
    max_candidates: int = 100,
) -> dict[str, Any]:
    dataset = store.get_dataset(dataset_id)
    reference_record = store.get_reference(reference_id)
    parameters_id = store.put_parameters(parameters)
    sequence_type = dataset["sequence_type"]
    if sequence_type != reference_record["sequence_type"]:
        raise ValueError("dataset and reference sequence_type differ")

    reference = normalize_sequence(reference_record["sequence"])
    reference_errors = validate_sequence(reference, sequence_type)
    if reference_errors:
        raise ValueError("invalid reference: " + "; ".join(reference_errors))

    seen_ids: set[str] = set()
    results: list[dict[str, Any]] = []
    quarantine: list[dict[str, Any]] = []

    for item in dataset["sequences"]:
        query_id = str(item["id"])
        raw_sequence = str(item["sequence"])
        query = normalize_sequence(raw_sequence)
        errors: list[str] = []
        if not query_id:
            errors.append("sequence id is empty")
        elif query_id in seen_ids:
            errors.append("duplicate sequence id")
        else:
            seen_ids.add(query_id)
        errors.extend(validate_sequence(raw_sequence, sequence_type))
        if errors:
            quarantine.append(
                {
                    "sequence_id": query_id,
                    "raw_sequence": raw_sequence,
                    "errors": errors,
                }
            )
            continue

        candidates = smith_waterman(
            query,
            reference,
            parameters,
            max_candidates=max_candidates,
        )
        results.append(
            {
                "sequence_id": query_id,
                "sequence": query,
                "candidate_count": len(candidates),
                "candidate_limit": max_candidates,
                "candidates_truncated": len(candidates) == max_candidates,
                "candidates": [asdict(candidate) for candidate in candidates],
            }
        )

    batch_id = _batch_id(dataset_id, reference_id, parameters_id, max_candidates)
    batch_directory = store.batch_dir / batch_id
    batch_directory.mkdir(parents=True, exist_ok=True)
    batch = {
        "schema_version": 1,
        "id": batch_id,
        "algorithm": "smith-waterman-linear-gap",
        "dataset_id": dataset_id,
        "reference_id": reference_id,
        "parameters_id": parameters_id,
        "max_candidates": max_candidates,
        "valid_sequence_count": len(results),
        "quarantined_sequence_count": len(quarantine),
    }
    write_json(store.batch_path(batch_id), batch)
    write_json(store.results_path(batch_id), {"batch_id": batch_id, "results": results})
    write_json(
        store.quarantine_path(batch_id),
        {"batch_id": batch_id, "quarantine": quarantine},
    )
    return {
        "batch": batch,
        "results": results,
        "quarantine": quarantine,
        "path": batch_directory,
    }


def load_batch(store: PlatformStore, batch_id: str) -> dict[str, Any]:
    return {
        "batch": read_json(store.batch_path(batch_id)),
        "results": read_json(store.results_path(batch_id)),
        "quarantine": read_json(store.quarantine_path(batch_id)),
    }
