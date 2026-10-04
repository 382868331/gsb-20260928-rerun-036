from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioseq_platform.storage import ProjectStore
from bioseq_platform.versions import (
    add_version,
    export_mappings,
    undo_last_version,
    verify_chain,
)

SEQUENCE = "AAAACCCCGGGG"  # length 12


class VersionChainTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.store = ProjectStore(Path(self.tmp.name) / "proj")
        self.store.initialize()
        self.reference = self.store.add_reference(SEQUENCE, name="root")
        self.feature = self.store.add_feature(
            self.reference.id, name="wrap-feature", start=10, end=3, strand="-"
        )

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def _mapping(self, record: dict, feature_id: str | None = None) -> dict:
        feature_id = feature_id or self.feature.id
        return next(m for m in record["feature_mappings"] if m["feature_id"] == feature_id)

    def test_rotation_preserves_cross_origin_segments_and_order(self) -> None:
        record = add_version(self.store, self.reference.id, {"type": "rotate", "offset": 2})
        mapping = self._mapping(record)
        self.assertEqual(mapping["status"], "mapped")
        self.assertTrue(mapping["wraps_origin"])
        self.assertEqual(
            [(s["start"], s["end"]) for s in mapping["segments"]],
            [(8, 12), (1, 1)],
        )
        self.assertEqual(
            [p for s in mapping["segments"] for p in s["positions"]],
            [8, 9, 10, 11, 12, 1],
        )
        self.assertEqual(
            mapping["original"],
            {"start": 10, "end": 3, "strand": "-", "wraps_origin": True},
        )

    def test_rotate_delete_boundary_insert_chain(self) -> None:
        add_version(self.store, self.reference.id, {"type": "rotate", "offset": 4})
        deleted = add_version(
            self.store, self.reference.id, {"type": "delete", "start": 8, "end": 9}
        )
        mapping = self._mapping(deleted)
        self.assertEqual(mapping["status"], "partially_deleted")
        self.assertEqual(
            [p for s in mapping["segments"] for p in s["positions"]], [6, 7, 8, 9]
        )
        self.assertEqual(
            [u["origin"] for u in mapping["deleted_units"]], [12, 1]
        )
        self.assertEqual(mapping["impacts"][0]["location"], "partial")

        inserted = add_version(
            self.store,
            self.reference.id,
            {"type": "insert", "after_position": 9, "bases": "TT"},
        )
        mapping = self._mapping(inserted)
        self.assertEqual(mapping["impacts"][0]["location"], "boundary")
        self.assertEqual(
            [u["source"] for u in mapping["units"]], ["original"] * 4
        )

    def test_internal_insertion_attribution_and_segment_sources(self) -> None:
        add_version(self.store, self.reference.id, {"type": "rotate", "offset": 4})
        record = add_version(
            self.store,
            self.reference.id,
            {"type": "insert", "after_position": 7, "bases": "AA"},
        )
        mapping = self._mapping(record)
        self.assertEqual(mapping["impacts"][0]["location"], "internal")
        self.assertEqual(
            [u["source"] for u in mapping["units"]],
            ["original", "original", "inserted", "inserted", "original", "original", "original", "original"],
        )
        segment = mapping["segments"][0]
        self.assertEqual(segment["positions"], [6, 7, 8, 9, 10, 11, 12, 13])
        self.assertEqual(
            segment["sources"],
            ["original", "original", "inserted", "inserted", "original", "original", "original", "original"],
        )
        self.assertEqual(
            [u["origin"] for u in mapping["units"] if u["source"] == "inserted"],
            [f"{record['version_id']}:1", f"{record['version_id']}:2"],
        )

    def test_full_deletion_marked_deleted(self) -> None:
        feature = self.store.add_feature(self.reference.id, name="core", start=2, end=4)
        record = add_version(
            self.store, self.reference.id, {"type": "delete", "start": 1, "end": 5}
        )
        mapping = self._mapping(record, feature.id)
        self.assertEqual(mapping["status"], "deleted")
        self.assertEqual(mapping["segments"], [])
        self.assertEqual(mapping["impacts"][0]["location"], "covers_feature")

    def test_export_intermediate_and_undo_last(self) -> None:
        add_version(self.store, self.reference.id, {"type": "rotate", "offset": 4})
        add_version(self.store, self.reference.id, {"type": "delete", "start": 8, "end": 9})
        add_version(
            self.store, self.reference.id, {"type": "insert", "after_position": 9, "bases": "TT"}
        )
        exported = export_mappings(self.store, self.reference.id, version_index=1)
        self.assertEqual(exported["exported_version_count"], 2)
        self.assertEqual(exported["total_version_count"], 3)
        self.assertEqual(exported["versions"][-1]["operation"]["type"], "delete")

        removed = undo_last_version(self.store, self.reference.id)
        self.assertEqual(removed["operation"]["type"], "insert")
        chain = export_mappings(self.store, self.reference.id)
        self.assertEqual(chain["total_version_count"], 2)
        self.assertEqual(chain["versions"][-1]["operation"]["type"], "delete")

    def test_save_reopen_preserves_segment_sources_and_attribution(self) -> None:
        add_version(self.store, self.reference.id, {"type": "rotate", "offset": 4})
        add_version(
            self.store,
            self.reference.id,
            {"type": "insert", "after_position": 7, "bases": "AA"},
        )
        reopened = ProjectStore(self.store.root)
        chain = reopened.get_version_chain(self.reference.id)
        mapping = chain["versions"][-1]["feature_mappings"][0]
        self.assertEqual(
            mapping["segments"][0]["sources"],
            ["original", "original", "inserted", "inserted", "original", "original", "original", "original"],
        )
        self.assertEqual(mapping["impacts"][0]["location"], "internal")
        self.assertEqual(mapping["original"], {"start": 10, "end": 3, "strand": "-", "wraps_origin": True})

    def test_verify_direct_mapping_matches_per_version_composition(self) -> None:
        add_version(self.store, self.reference.id, {"type": "rotate", "offset": 4})
        add_version(self.store, self.reference.id, {"type": "delete", "start": 8, "end": 9})
        add_version(
            self.store, self.reference.id, {"type": "insert", "after_position": 7, "bases": "AA"}
        )
        add_version(self.store, self.reference.id, {"type": "rotate", "offset": 3})
        report = verify_chain(self.store, self.reference.id)
        self.assertTrue(report["sequences_match"])
        self.assertTrue(report["ok"])
        feature_report = report["features"][0]
        self.assertTrue(feature_report["direct_matches_composed"])
        self.assertTrue(feature_report["inserted_positions_match"])

    def test_verify_detects_tampered_sequence(self) -> None:
        record = add_version(self.store, self.reference.id, {"type": "rotate", "offset": 4})
        path = self.store.references_dir / f"{record['reference_id']}.json"
        import json

        data = json.loads(path.read_text(encoding="utf-8"))
        data["sequence"] = "TTTT" + data["sequence"][4:]
        path.write_text(json.dumps(data), encoding="utf-8")
        report = verify_chain(self.store, self.reference.id)
        self.assertFalse(report["sequences_match"])
        self.assertFalse(report["ok"])


if __name__ == "__main__":
    unittest.main()
