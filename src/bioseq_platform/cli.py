from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import dataclass_list
from .models import AlignmentParameters, Scoring
from .pipeline import (
    project_run_record,
    run_from_files,
    select_alignment,
    selection_projection_diff,
)
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
    select_parser.add_argument("--preview", action="store_true", help="only preview which feature projections would be recomputed")
    feature_add_parser = subparsers.add_parser("feature-add", help="annotate a directed interval on a reference")
    feature_add_parser.add_argument("project", help="project directory")
    feature_add_parser.add_argument("reference_id", help="reference identifier")
    feature_add_parser.add_argument("--name", required=True, help="feature name")
    feature_add_parser.add_argument("--start", type=int, required=True, help="1-based start (inclusive)")
    feature_add_parser.add_argument("--end", type=int, required=True, help="1-based end (inclusive); start > end wraps the origin")
    feature_add_parser.add_argument("--strand", choices=["+", "-"], default="+")
    feature_list_parser = subparsers.add_parser("feature-list", help="list features annotated on a reference")
    feature_list_parser.add_argument("project", help="project directory")
    feature_list_parser.add_argument("reference_id", help="reference identifier")
    project_parser = subparsers.add_parser("project", help="project annotated features onto query for each candidate")
    project_parser.add_argument("project", help="project directory")
    project_parser.add_argument("run_id", help="run identifier")
    project_parser.add_argument("record_id", help="record identifier within the run")
    project_parser.add_argument("--index", type=int, default=None, help="restrict to one 0-based candidate index")
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
        if args.preview:
            print(json.dumps(selection_projection_diff(store, args.run_id, args.record_id, args.index), ensure_ascii=False, indent=2))
            return 0
        preview = selection_projection_diff(store, args.run_id, args.record_id, args.index)
        results = select_alignment(store, args.run_id, args.record_id, args.index)
        print(json.dumps({
            "run_id": args.run_id,
            "record_id": args.record_id,
            "selected_alignment_index": args.index,
            "features_recomputed": preview["features_to_recompute"],
            "results": dataclass_list(results),
        }, ensure_ascii=False, indent=2))
        return 0
    if args.command == "feature-add":
        store.initialize()
        feature = store.add_feature(
            args.reference_id,
            name=args.name,
            start=args.start,
            end=args.end,
            strand=args.strand,
        )
        print(json.dumps(feature.to_dict(), ensure_ascii=False, indent=2))
        return 0
    if args.command == "feature-list":
        print(json.dumps([feature.to_dict() for feature in store.get_features(args.reference_id)], ensure_ascii=False, indent=2))
        return 0
    if args.command == "project":
        print(json.dumps(project_run_record(store, args.run_id, args.record_id, args.index), ensure_ascii=False, indent=2))
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
