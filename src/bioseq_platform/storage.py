from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from .io import (
    atomic_write_json,
    dataclass_list,
    read_json,
    read_reference,
    result_from_dict,
    stable_hash,
    write_records,
    write_reference,
)
from .models import AlignmentParameters, InvalidRecord, Reference, SequenceRecord, SequenceResult
from .validation import normalize_sequence


class ProjectStore:
    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.references_dir = self.root / "references"
        self.datasets_dir = self.root / "datasets"
        self.parameters_dir = self.root / "parameters"
        self.runs_dir = self.root / "runs"

    def initialize(self) -> None:
        for directory in [self.root, self.references_dir, self.datasets_dir, self.parameters_dir, self.runs_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        manifest = self.root / "project.json"
        if not manifest.exists():
            atomic_write_json(manifest, {
                "platform": "bioseq-platform",
                "version": 1,
                "online_databases": False,
            })

    @staticmethod
    def _content_hash(payload: dict) -> str:
        return stable_hash(payload)[:16]

    def add_reference(self, sequence: str, *, sequence_type: str = "dna", name: str = "", metadata: dict | None = None, reference_id: str | None = None) -> Reference:
        normalized = normalize_sequence(sequence)
        reference = Reference(
            id=reference_id or f"ref-{stable_hash([sequence_type, normalized, name, metadata or {}])[:12]}",
            sequence=normalized,
            name=name,
            sequence_type=sequence_type,
            metadata=metadata or {},
        )
        write_reference(self.references_dir / f"{reference.id}.json", reference)
        return reference

    def get_reference(self, reference_id: str) -> Reference:
        return read_reference(self.references_dir / f"{reference_id}.json")

    def add_dataset(self, valid: list[SequenceRecord], invalid: list[InvalidRecord], *, dataset_id: str | None = None) -> str:
        payload = {"valid": [asdict(item) for item in valid], "invalid": [asdict(item) for item in invalid]}
        dataset_hash = dataset_id or f"ds-{self._content_hash(payload)}"
        write_records(self.datasets_dir / f"{dataset_hash}.json", valid, invalid)
        return dataset_hash

    def get_dataset(self, dataset_id: str) -> tuple[list[SequenceRecord], list[InvalidRecord]]:
        data = read_json(self.datasets_dir / f"{dataset_id}.json")
        return [SequenceRecord(**item) for item in data["valid"]], [InvalidRecord(**item) for item in data["invalid"]]

    def add_parameters(self, parameters: AlignmentParameters, *, parameters_id: str | None = None) -> str:
        parameters.validate()
        payload = parameters.to_dict()
        parameter_hash = parameters_id or f"param-{self._content_hash(payload)}"
        atomic_write_json(self.parameters_dir / f"{parameter_hash}.json", payload)
        return parameter_hash

    def get_parameters(self, parameters_id: str) -> AlignmentParameters:
        return AlignmentParameters.from_dict(read_json(self.parameters_dir / f"{parameters_id}.json"))

    def run_path(self, run_id: str) -> Path:
        return self.runs_dir / run_id

    def run_exists(self, run_id: str) -> bool:
        return (self.run_path(run_id) / "manifest.json").exists()

    def create_run(self, run_id: str, manifest: dict, results: list[SequenceResult]) -> None:
        run_dir = self.run_path(run_id)
        atomic_write_json(run_dir / "manifest.json", manifest)
        atomic_write_json(run_dir / "results.json", dataclass_list(results))

    def get_run_manifest(self, run_id: str) -> dict:
        return read_json(self.run_path(run_id) / "manifest.json")

    def get_run_results(self, run_id: str) -> list[SequenceResult]:
        return [result_from_dict(item) for item in read_json(self.run_path(run_id) / "results.json")]

    def save_run_results(self, run_id: str, results: list[SequenceResult]) -> None:
        if not self.run_exists(run_id):
            raise KeyError(f"unknown run: {run_id}")
        atomic_write_json(self.run_path(run_id) / "results.json", dataclass_list(results))
