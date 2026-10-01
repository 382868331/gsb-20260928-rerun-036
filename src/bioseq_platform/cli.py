from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import dataclass_list
from .models import AlignmentParameters, Scoring
from .pipeline import run_from_files, select_alignment
from .storage import ProjectStore


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bioseq", description="Offline biological sequence alignment")
    subparsers = parser.add_subparsers(dest="command", required=True)
    init_parser = subparsers.add_parser("init", help="initialize a project directory")
    init_parser.add_argument("project", help="project directory")

    align_parser = subparsers.add_parser("align", help="run an offline alignment batch")
    align_parser.add_argument("project", help="project directory")
    align_parser.add_argument("--reference", required=True, help="one-record FASTA/JSON reference file")
    align_parser.add_argument("--sequences", required=True, help="FASTA/JSON/JSONL/CSV sequences")
    align_parser.add_argument("--type", choices=["dna", "rna", "protein"], default="dna")
    align_parser.add_argument("--mode", choices=["local", "global"], default="local")
    align_parser.add_argument("--topology", choices=["linear", "circular"], default="linear")
    align_parser.add_argument("--match", type=int, default=2)
    align_parser.add_argument("--mismatch", type=int, default=-1)
    align_parser.add_argument("--gap-open", type=int, default=-2)
    align_parser.add_argument("--gap-extend", type=int, default=-2)
    align_parser.add_argument("--max-alignments", type=int, default=100)

    show_parser = subparsers.add_parser("show-run", help="print stored run results as JSON")
    show_parser.add_argument("project", help="project directory")
    show_parser.add_argument("run_id", help="run identifier, e.g. run-...")
    select_parser = subparsers.add_parser("select", help="mark one tied candidate as the chosen alignment")
    select_parser.add_argument("project", help="project directory")
    select_parser.add_argument("run_id", help="run identifier")
    select_parser.add_argument("record_id", help="record identifier within the run")
    select_parser.add_argument("index", type=int, help="0-based candidate index to select")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    store = ProjectStore(args.project)
    if args.command == "init":
        store.initialize()
        print(json.dumps({"project": str(Path(args.project).resolve()), "status": "initialized"}, ensure_ascii=False))
        return 0
    if args.command == "align":
        store.initialize()
        parameters = AlignmentParameters(
            mode=args.mode,
            topology=args.topology,
            scoring=Scoring(match=args.match, mismatch=args.mismatch, gap_open=args.gap_open, gap_extend=args.gap_extend),
            max_alignments=args.max_alignments,
        )
        run_id, manifest, results = run_from_files(
            store,
            reference_path=args.reference,
            sequences_path=args.sequences,
            sequence_type=args.type,
            parameters=parameters,
        )
        print(json.dumps({
            "run_id": run_id,
            "reference_id": manifest["reference_id"],
            "dataset_id": manifest["dataset_id"],
            "parameters_id": manifest["parameters_id"],
            "results": dataclass_list(results),
        }, ensure_ascii=False, indent=2))
        return 0
    if args.command == "show-run":
        print(json.dumps({
            "manifest": store.get_run_manifest(args.run_id),
            "results": dataclass_list(store.get_run_results(args.run_id)),
        }, ensure_ascii=False, indent=2))
        return 0
    if args.command == "select":
        results = select_alignment(store, args.run_id, args.record_id, args.index)
        print(json.dumps({
            "run_id": args.run_id,
            "record_id": args.record_id,
            "selected_alignment_index": args.index,
            "results": dataclass_list(results),
        }, ensure_ascii=False, indent=2))
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
