from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import asdict
from pathlib import Path
from typing import Any

from .models import Alignment, AlignmentParameters, InvalidRecord, Reference, SequenceRecord, SequenceResult


def stable_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def stable_hash(value: Any) -> str:
    return hashlib.sha256(stable_json(value).encode("utf-8")).hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_fasta(path: Path) -> list[dict]:
    records: list[dict] = []
    identifier = ""
    description = ""
    chunks: list[str] = []

    def flush() -> None:
        if identifier or chunks:
            records.append({"id": identifier or f"record-{len(records) + 1:04d}", "name": description, "sequence": "".join(chunks)})

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].split(maxsplit=1)
            identifier = header[0] if header else ""
            description = header[1] if len(header) > 1 else ""
            chunks = []
        else:
            chunks.append(line)
    flush()
    return records


def read_sequence_file(path: Path) -> list[dict]:
    suffix = path.suffix.lower()
    if suffix == ".json":
        data = read_json(path)
        if isinstance(data, dict):
            data = data.get("records", [])
        if not isinstance(data, list):
            raise ValueError("JSON input must be a list or {'records': [...]}")
        return data
    if suffix == ".jsonl":
        records = []
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    records.append(json.loads(line))
        return records
    if suffix == ".csv":
        with path.open(newline="", encoding="utf-8") as handle:
            return [{key: value for key, value in row.items() if value != ""} for row in csv.DictReader(handle)]
    if suffix in {".fa", ".fasta", ".fna", ".faa"}:
        return parse_fasta(path)
    raise ValueError(f"unsupported input format: {path.suffix}")


def dataclass_list(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return [dataclass_list(item) for item in value]
    if hasattr(value, "__dataclass_fields__"):
        return {key: dataclass_list(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {key: dataclass_list(item) for key, item in value.items()}
    return value


def alignment_from_dict(data: dict) -> Alignment:
    return Alignment(
        score=data["score"],
        reference_start=data["reference_start"],
        reference_end=data["reference_end"],
        query_start=data["query_start"],
        query_end=data["query_end"],
        aligned_reference=data["aligned_reference"],
        aligned_query=data["aligned_query"],
        cigar=data["cigar"],
        variants=tuple(data.get("variants", [])),
    )


def result_from_dict(data: dict) -> SequenceResult:
    return SequenceResult(
        record_id=data["record_id"],
        status=data["status"],
        score=data["score"],
        alignment_count=data["alignment_count"],
        alignments=tuple(alignment_from_dict(item) for item in data.get("alignments", [])),
        alignment_limit_reached=data.get("alignment_limit_reached", False),
    )


def write_reference(path: Path, reference: Reference) -> None:
    atomic_write_json(path, asdict(reference))


def read_reference(path: Path) -> Reference:
    return Reference(**read_json(path))


def write_records(path: Path, valid: list[SequenceRecord], invalid: list[InvalidRecord]) -> None:
    atomic_write_json(path, {"valid": [asdict(item) for item in valid], "invalid": [asdict(item) for item in invalid]})
