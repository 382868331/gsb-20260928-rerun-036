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
from .features import Feature
from .validation import normalize_sequence
from . import versioning


class ProjectStore:
    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.references_dir = self.root / "references"
        self.datasets_dir = self.root / "datasets"
        self.parameters_dir = self.root / "parameters"
        self.runs_dir = self.root / "runs"
        self.features_dir = self.root / "features"
        self.chains_dir = self.root / "chains"

    def initialize(self) -> None:
        for directory in [self.root, self.references_dir, self.datasets_dir, self.parameters_dir, self.runs_dir, self.features_dir, self.chains_dir]:
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

    def _features_path(self, reference_id: str) -> Path:
        return self.features_dir / f"{reference_id}.json"

    def get_features(self, reference_id: str) -> list[Feature]:
        path = self._features_path(reference_id)
        if not path.exists():
            return []
        return [Feature.from_dict(item) for item in read_json(path).get("features", [])]

    def save_features(self, reference_id: str, features: list[Feature]) -> None:
        atomic_write_json(self._features_path(reference_id), {
            "reference_id": reference_id,
            "features": [feature.to_dict() for feature in features],
        })

    def add_feature(self, reference_id: str, *, name: str, start: int, end: int, strand: str = "+", metadata: dict | None = None) -> Feature:
        reference = self.get_reference(reference_id)
        feature = Feature(
            id=f"feat-{stable_hash([reference_id, name, start, end, strand, metadata or {}])[:10]}",
            reference_id=reference_id,
            name=name,
            start=start,
            end=end,
            strand=strand,
            metadata=metadata or {},
        )
        feature.positions(len(reference.sequence))
        features = self.get_features(reference_id)
        if any(existing.id == feature.id for existing in features):
            return feature
        features.append(feature)
        self.save_features(reference_id, features)
        return feature

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

    def _chain_path(self, chain_id: str) -> Path:
        return self.chains_dir / f"{chain_id}.json"

    def chain_exists(self, chain_id: str) -> bool:
        return self._chain_path(chain_id).exists()

    def save_chain(self, chain: dict) -> None:
        atomic_write_json(self._chain_path(chain["id"]), chain)

    def get_chain(self, chain_id: str) -> dict:
        path = self._chain_path(chain_id)
        if not path.exists():
            raise KeyError(f"unknown chain: {chain_id}")
        return read_json(path)

    def create_chain(self, reference_id: str, chain_id: str | None = None) -> dict:
        reference = self.get_reference(reference_id)
        chain = versioning.create_chain(reference.id, reference.sequence, chain_id)
        if self.chain_exists(chain["id"]):
            raise ValueError(f"chain already exists: {chain['id']}")
        self.save_chain(chain)
        return chain

    def update_chain(self, chain_id: str, operation: str, **kwargs) -> dict:
        chain = self.get_chain(chain_id)
        chain = versioning.apply_operation(chain, operation, **kwargs)
        self.save_chain(chain)
        return chain

    def confirm_edit_plan(self, chain_id: str, plan: dict, new_chain_id: str) -> dict:
        chain = self.get_chain(chain_id)
        confirmed = versioning.apply_confirmed_edit_plan(chain, plan, new_chain_id)
        if self.chain_exists(new_chain_id):
            raise ValueError(f"chain already exists: {new_chain_id}")
        self.save_chain(confirmed)
        return confirmed
