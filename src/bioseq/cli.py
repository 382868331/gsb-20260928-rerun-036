from __future__ import annotations

import argparse
import json
from pathlib import Path

from .alignment import AlignmentParameters
from .runner import load_batch, run_batch
from .storage import PlatformStore, load_dataset, load_reference


def _json_arg(value: str | int | dict | list) -> str:
    return json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bioseq")
    parser.add_argument("--store", default=".bioseq", help="platform data directory")
    subparsers = parser.add_subparsers(dest="command", required=True)

    align = subparsers.add_parser("align", help="store inputs and run a reproducible batch")
    align.add_argument("--reference", required=True)
    align.add_argument("--dataset", required=True)
    align.add_argument("--match", type=int, default=2)
    align.add_argument("--mismatch", type=int, default=-1)
    align.add_argument("--gap-penalty", type=int, default=2)
    align.add_argument("--max-candidates", type=int, default=100)

    show = subparsers.add_parser("show", help="show a stored batch")
    show.add_argument("batch_id")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    store = PlatformStore(Path(args.store))

    if args.command == "align":
        reference = store.put_reference(load_reference(Path(args.reference)))
        dataset = store.put_dataset(load_dataset(Path(args.dataset)))
        parameters = AlignmentParameters(
            match_score=args.match,
            mismatch_score=args.mismatch,
            gap_penalty=args.gap_penalty,
        )
        outcome = run_batch(
            store,
            dataset,
            reference,
            parameters,
            max_candidates=args.max_candidates,
        )
        print(_json_arg(outcome["batch"]))
        print(f"results: {outcome['path'] / 'results.json'}")
        print(f"quarantine: {outcome['path'] / 'quarantine.json'}")
        return 0

    if args.command == "show":
        print(_json_arg(load_batch(store, args.batch_id)))
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
