from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import dataclass_list
from .models import AlignmentParameters, Scoring
from .pipeline import annotation_from_file, build_projections, projection_preview, run_from_files, select_alignment
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

    annotate_parser = subparsers.add_parser("annotate", help="save directed functional regions for a reference")
    annotate_parser.add_argument("project", help="project directory")
    annotate_parser.add_argument("--annotation", required=True, help="JSON/TSV functional region file")
    annotate_parser.add_argument("--reference-id", required=True, help="reference the regions belong to")

    project_parser = subparsers.add_parser("project", help="project an annotation across all alignment interpretations")
    project_parser.add_argument("project", help="project directory")
    project_parser.add_argument("run_id", help="run identifier")
    project_parser.add_argument("record_id", help="record identifier within the run")
    project_parser.add_argument("annotation_id", help="saved annotation identifier")

    preview_parser = subparsers.add_parser("preview-reselect", help="preview which projections change under another candidate")
    preview_parser.add_argument("project", help="project directory")
    preview_parser.add_argument("run_id", help="run identifier")
    preview_parser.add_argument("record_id", help="record identifier within the run")
    preview_parser.add_argument("annotation_id", help="saved annotation identifier")
    preview_parser.add_argument("index", type=int, help="candidate index that may be selected next")

    show_projection_parser = subparsers.add_parser("show-projection", help="print stored annotation projections")
    show_projection_parser.add_argument("project", help="project directory")
    show_projection_parser.add_argument("run_id", help="run identifier")
    show_projection_parser.add_argument("record_id", help="record identifier within the run")
    show_projection_parser.add_argument("annotation_id", help="saved annotation identifier")
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
    if args.command == "annotate":
        annotation_id = annotation_from_file(
            store,
            annotation_path=args.annotation,
            reference_id=args.reference_id,
        )
        print(json.dumps({
            "annotation_id": annotation_id,
            "reference_id": args.reference_id,
        }, ensure_ascii=False, indent=2))
        return 0
    if args.command == "project":
        payload = build_projections(
            store,
            run_id=args.run_id,
            record_id=args.record_id,
            annotation_id=args.annotation_id,
        )
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    if args.command == "preview-reselect":
        preview = projection_preview(
            store,
            run_id=args.run_id,
            record_id=args.record_id,
            annotation_id=args.annotation_id,
            target_alignment_index=args.index,
        )
        print(json.dumps(preview, ensure_ascii=False, indent=2))
        return 0
    if args.command == "show-projection":
        payload = store.load_projection(args.run_id, args.record_id, args.annotation_id)
        if payload is None:
            raise SystemExit(
                f"no projection for annotation {args.annotation_id} on {args.record_id}"
            )
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
