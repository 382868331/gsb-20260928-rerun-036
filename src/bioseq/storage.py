from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from .alignment import AlignmentParameters


def canonical_json(data: Any) -> str:
    return json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def stable_hash(data: Any) -> str:
    return hashlib.sha256(canonical_json(data).encode("utf-8")).hexdigest()


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


class PlatformStore:
    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)
        self.reference_dir = self.root / "references"
        self.dataset_dir = self.root / "datasets"
        self.parameter_dir = self.root / "parameters"
        self.batch_dir = self.root / "batches"
        for directory in (
            self.root,
            self.reference_dir,
            self.dataset_dir,
            self.parameter_dir,
            self.batch_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)

    def reference_path(self, reference_id: str) -> Path:
        return self.reference_dir / f"{reference_id}.json"

    def dataset_path(self, dataset_id: str) -> Path:
        return self.dataset_dir / f"{dataset_id}.json"

    def parameters_path(self, parameters_id: str) -> Path:
        return self.parameter_dir / f"{parameters_id}.json"

    def batch_path(self, batch_id: str) -> Path:
        return self.batch_dir / batch_id / "batch.json"

    def results_path(self, batch_id: str) -> Path:
        return self.batch_dir / batch_id / "results.json"

    def quarantine_path(self, batch_id: str) -> Path:
        return self.batch_dir / batch_id / "quarantine.json"

    def put_reference(self, reference: dict[str, Any]) -> str:
        reference_id = reference.get("id") or "ref_" + stable_hash(
            {
                "sequence_type": reference["sequence_type"],
                "sequence": reference["sequence"],
            }
        )[:16]
        record = {
            "id": reference_id,
            "name": reference.get("name", reference_id),
            "sequence_type": reference["sequence_type"],
            "sequence": normalize_sequence(reference["sequence"]),
            "sequence_sha256": stable_hash(
                normalize_sequence(reference["sequence"])
            ),
        }
        write_json(self.reference_path(reference_id), record)
        return reference_id

    def put_dataset(self, dataset: dict[str, Any]) -> str:
        sequences = [
            {
                "id": item["id"],
                "sequence": normalize_sequence(item["sequence"]),
            }
            for item in dataset["sequences"]
        ]
        dataset_id = dataset.get("id") or "ds_" + stable_hash(
            {
                "sequence_type": dataset["sequence_type"],
                "sequences": sequences,
            }
        )[:16]
        record = {
            "id": dataset_id,
            "name": dataset.get("name", dataset_id),
            "sequence_type": dataset["sequence_type"],
            "sequences": sequences,
        }
        write_json(self.dataset_path(dataset_id), record)
        return dataset_id

    def put_parameters(self, parameters: AlignmentParameters) -> str:
        payload = {
            "match_score": parameters.match_score,
            "mismatch_score": parameters.mismatch_score,
            "gap_penalty": parameters.gap_penalty,
        }
        parameters_id = "par_" + stable_hash(payload)[:16]
        write_json(
            self.parameters_path(parameters_id),
            {"id": parameters_id, **payload},
        )
        return parameters_id

    def get_reference(self, reference_id: str) -> dict[str, Any]:
        return read_json(self.reference_path(reference_id))

    def get_dataset(self, dataset_id: str) -> dict[str, Any]:
        return read_json(self.dataset_path(dataset_id))

    def get_parameters(self, parameters_id: str) -> AlignmentParameters:
        payload = read_json(self.parameters_path(parameters_id))
        return AlignmentParameters(
            match_score=payload["match_score"],
            mismatch_score=payload["mismatch_score"],
            gap_penalty=payload["gap_penalty"],
        )


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_dataset(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    if source.suffix.lower() == ".json":
        data = load_json(source)
        return {
            "name": data.get("name", source.stem),
            "sequence_type": data["sequence_type"],
            "sequences": list(data["sequences"]),
        }
    if source.suffix.lower() != ".csv":
        raise ValueError("dataset must be a .json or .csv file")
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("dataset CSV is empty")
    sequence_type = rows[0].get("sequence_type", "dna")
    return {
        "name": source.stem,
        "sequence_type": sequence_type,
        "sequences": [
            {"id": row["id"], "sequence": row["sequence"]} for row in rows
        ],
    }


def load_reference(path: str | Path) -> dict[str, Any]:
    data = load_json(Path(path))
    return {
        "id": data.get("id"),
        "name": data.get("name", Path(path).stem),
        "sequence_type": data["sequence_type"],
        "sequence": data["sequence"],
    }
from .validation import normalize_sequence
