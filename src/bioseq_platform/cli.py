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
from . import combination, versioning


def _parse_ranges(text):
    ranges = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        start, end = item.split(":")
        ranges.append((int(start), int(end)))
    return ranges


def _combo_inputs(args):
    with open(args.candidates, "r", encoding="utf-8") as handle:
        candidates = json.load(handle)
    keep = None
    if args.keep:
        keep = [item.strip() for item in args.keep.split(",") if item.strip()]
    ranges = _parse_ranges(args.allow) if args.allow else []
    return candidates, keep, ranges


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

    chain_create = subparsers.add_parser("chain-create", help="start a reference version chain")
    chain_create.add_argument("project", help="project directory")
    chain_create.add_argument("reference_id", help="reference identifier for version 0")
    chain_create.add_argument("--chain-id", default=None, help="explicit chain identifier")
    chain_edit = subparsers.add_parser("chain-edit", help="apply a rotate/delete/insert edit as a new version")
    chain_edit.add_argument("project", help="project directory")
    chain_edit.add_argument("chain_id", help="chain identifier")
    chain_edit.add_argument("operation", choices=["rotate", "delete", "insert"])
    chain_edit.add_argument("--offset", type=int, default=None, help="rotation offset (rotate)")
    chain_edit.add_argument("--start", type=int, default=None, help="1-based deletion start (delete)")
    chain_edit.add_argument("--end", type=int, default=None, help="1-based deletion end (delete)")
    chain_edit.add_argument("--after", type=int, default=None, help="insert after this position; 0 = before first (insert)")
    chain_edit.add_argument("--bases", default=None, help="inserted bases (insert)")
    chain_undo = subparsers.add_parser("chain-undo", help="withdraw the last edit (drop the last version)")
    chain_undo.add_argument("project", help="project directory")
    chain_undo.add_argument("chain_id", help="chain identifier")
    chain_trace = subparsers.add_parser("chain-trace", help="export feature mapping at every intermediate version")
    chain_trace.add_argument("project", help="project directory")
    chain_trace.add_argument("chain_id", help="chain identifier")
    chain_trace.add_argument("feature_id", help="feature annotated on the original reference")
    chain_trace.add_argument("--version", type=int, default=None, help="stop tracing at this version index")
    chain_verify = subparsers.add_parser("chain-verify", help="check direct mapping against per-version composition")
    chain_verify.add_argument("project", help="project directory")
    chain_verify.add_argument("chain_id", help="chain identifier")
    combo_study = subparsers.add_parser("combo-study", help="study compatible combinations of candidate edits")
    combo_study.add_argument("project", help="project directory")
    combo_study.add_argument("chain_id", help="base chain identifier")
    combo_study.add_argument("--candidates", required=True, help="JSON file with candidate edits")
    combo_study.add_argument("--keep", default=None, help="comma-separated feature ids to keep (default: all)")
    combo_study.add_argument("--allow", default=None, help="comma-separated start:end ranges where disruption is allowed")
    combo_study.add_argument("--base-version", type=int, default=None, help="base version index (default: latest)")
    combo_confirm = subparsers.add_parser("combo-confirm", help="apply a chosen combination as an independent chain")
    combo_confirm.add_argument("project", help="project directory")
    combo_confirm.add_argument("chain_id", help="base chain identifier")
    combo_confirm.add_argument("--candidates", required=True, help="JSON file with candidate edits")
    combo_confirm.add_argument("--select", required=True, help="comma-separated candidate ids to apply")
    combo_confirm.add_argument("--keep", default=None, help="comma-separated feature ids to keep (default: all)")
    combo_confirm.add_argument("--allow", default=None, help="comma-separated start:end ranges where disruption is allowed")
    combo_confirm.add_argument("--base-version", type=int, default=None, help="base version index (default: latest)")
    combo_confirm.add_argument("--new-chain-id", default=None, help="identifier for the resulting chain")
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

    if args.command == "chain-create":
        store.initialize()
        chain = store.create_chain(args.reference_id, args.chain_id)
        print(json.dumps({"chain_id": chain["id"], "versions": len(chain["versions"])}, ensure_ascii=False, indent=2))
        return 0
    if args.command == "chain-edit":
        kwargs = {}
        if args.operation == "rotate":
            if args.offset is None:
                raise SystemExit("chain-edit rotate requires --offset")
            kwargs["offset"] = args.offset
        elif args.operation == "delete":
            if args.start is None or args.end is None:
                raise SystemExit("chain-edit delete requires --start and --end")
            kwargs["start"], kwargs["end"] = args.start, args.end
        elif args.operation == "insert":
            if args.after is None or args.bases is None:
                raise SystemExit("chain-edit insert requires --after and --bases")
            kwargs["after_position"], kwargs["bases"] = args.after, args.bases
        chain = store.update_chain(args.chain_id, args.operation, **kwargs)
        print(json.dumps({
            "chain_id": chain["id"],
            "versions": len(chain["versions"]),
            "edit": chain["versions"][-1]["edit"],
            "sequence": chain["versions"][-1]["sequence"],
        }, ensure_ascii=False, indent=2))
        return 0
    if args.command == "chain-undo":
        chain = store.update_chain(args.chain_id, "undo")
        print(json.dumps({
            "chain_id": chain["id"],
            "versions": len(chain["versions"]),
            "sequence": chain["versions"][-1]["sequence"],
        }, ensure_ascii=False, indent=2))
        return 0
    if args.command == "chain-trace":
        chain = store.get_chain(args.chain_id)
        features = {feature.id: feature for feature in store.get_features(chain["reference_id"])}
        if args.feature_id not in features:
            raise SystemExit(f"unknown feature for reference {chain['reference_id']}: {args.feature_id}")
        trace = versioning.trace_feature(chain, features[args.feature_id], args.version)
        print(json.dumps(trace, ensure_ascii=False, indent=2))
        return 0
    if args.command == "chain-verify":
        print(json.dumps(versioning.verify_chain(store.get_chain(args.chain_id)), ensure_ascii=False, indent=2))
        return 0
    if args.command == "combo-study":
        chain = store.get_chain(args.chain_id)
        features = store.get_features(chain["reference_id"])
        candidates, keep, ranges = _combo_inputs(args)
        study = combination.study_combinations(
            chain, features, keep, ranges, candidates, args.base_version
        )
        print(json.dumps(study, ensure_ascii=False, indent=2))
        return 0
    if args.command == "combo-confirm":
        chain = store.get_chain(args.chain_id)
        features = store.get_features(chain["reference_id"])
        candidates, keep, ranges = _combo_inputs(args)
        selected = [item.strip() for item in args.select.split(",") if item.strip()]
        result = combination.confirm_combination(
            chain, features, keep, ranges, candidates, selected,
            new_chain_id=args.new_chain_id, base_version=args.base_version,
        )
        if store.chain_exists(result["new_chain_id"]):
            raise SystemExit(f"chain already exists: {result['new_chain_id']}")
        store.save_chain(result.pop("chain"))
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
