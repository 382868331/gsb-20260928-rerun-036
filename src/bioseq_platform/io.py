from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import asdict
from pathlib import Path
from typing import Any

from .models import (
    Alignment,
    AlignmentParameters,
    FunctionalRegion,
    InvalidRecord,
    Reference,
    SequenceRecord,
    SequenceResult,
)


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


def read_annotation_file(path: Path) -> list[FunctionalRegion]:
    """Read directed functional regions from JSON or a simple TSV.

    TSV columns are ``id start end [strand] [name]`` with 1-based closed
    coordinates; ``end < start`` marks a region crossing the record start.
    Lines starting with ``#`` are ignored. JSON accepts either a list of
    regions or an object with a ``regions`` list.
    """
    suffix = path.suffix.lower()
    if suffix == ".json":
        data = read_json(path)
        if isinstance(data, dict):
            data = data.get("regions", [])
        if not isinstance(data, list):
            raise ValueError("annotation JSON must be a list or {'regions': [...]}")
        regions = [FunctionalRegion.from_dict(item) for item in data]
    elif suffix in {".tsv", ".bed", ".txt"}:
        regions = []
        for raw_line in path.read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            if len(fields) < 3:
                raise ValueError(f"annotation line needs at least id/start/end: {raw_line!r}")
            identifier, start, end = fields[0], fields[1], fields[2]
            strand = fields[3] if len(fields) > 3 else "forward"
            if strand in {"+", "f", "plus"}:
                strand = "forward"
            elif strand in {"-", "r", "minus"}:
                strand = "reverse"
            name = fields[4] if len(fields) > 4 else ""
            regions.append(FunctionalRegion(
                id=identifier,
                start=int(start),
                end=int(end),
                strand=strand,
                name=name,
            ))
    else:
        raise ValueError(f"unsupported annotation format: {path.suffix}")
    if not regions:
        raise ValueError("annotation set contains no regions")
    identifiers = [region.id for region in regions]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("functional region ids must be unique")
    return regions


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
        cigar=data["cigar"],
        strand=data.get("strand", "forward"),
        reference_start=data["reference_start"],
        reference_end=data["reference_end"],
        query_start=data["query_start"],
        query_end=data["query_end"],
        aligned_reference=data["aligned_reference"],
        aligned_query=data["aligned_query"],
        rotation_offset=data.get("rotation_offset", 0),
        wraps_reference=data.get("wraps_reference", False),
        variants=tuple(data.get("variants", [])),
        reference_coordinates=tuple(data.get("reference_coordinates", ())),
        query_coordinates=tuple(data.get("query_coordinates", ())),
        trace_end_reference=data.get("trace_end_reference", data["reference_end"]),
        trace_end_query=data.get("trace_end_query", data["query_end"]),
    )


def result_from_dict(data: dict) -> SequenceResult:
    return SequenceResult(
        record_id=data["record_id"],
        status=data["status"],
        score=data["score"],
        alignment_count=data["alignment_count"],
        alignments=tuple(alignment_from_dict(item) for item in data.get("alignments", [])),
        alignment_limit_reached=data.get("alignment_limit_reached", False),
        selected_alignment_index=data.get("selected_alignment_index"),
    )


def write_reference(path: Path, reference: Reference) -> None:
    atomic_write_json(path, asdict(reference))


def read_reference(path: Path) -> Reference:
    return Reference(**read_json(path))


def write_records(path: Path, valid: list[SequenceRecord], invalid: list[InvalidRecord]) -> None:
    atomic_write_json(path, {"valid": [asdict(item) for item in valid], "invalid": [asdict(item) for item in invalid]})
