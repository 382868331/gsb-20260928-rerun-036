from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioseq_platform import versioning
from bioseq_platform.features import Feature
from bioseq_platform.storage import ProjectStore

REFERENCE = "ACGTACGTACGT"  # length 12


def make_feature(start=10, end=3, strand="-"):
    return Feature(
        id="feat-wrap", reference_id="ref-x", name="wrap",
        start=start, end=end, strand=strand,
    )


def build_chain():
    chain = versioning.create_chain("ref-x", REFERENCE, chain_id="chain-t")
    chain = versioning.rotate(chain, 4)                      # v1
    chain = versioning.delete(chain, 8, 9)                   # v2: partial feature deletion
    chain = versioning.insert(chain, 9, "GG")                # v3: boundary insertion
    return chain


class ChainEditTest(unittest.TestCase):
    def test_rotate_then_edit_sequences(self):
        chain = build_chain()
        self.assertEqual(chain["versions"][1]["sequence"], "ACGTACGT" + "ACGT")
        self.assertEqual(chain["versions"][1]["sequence"], REFERENCE[4:] + REFERENCE[:4])
        self.assertEqual(len(chain["versions"][2]["sequence"]), 10)
        self.assertEqual(len(chain["versions"][3]["sequence"]), 12)

    def test_inserted_provenance_attribution(self):
        chain = build_chain()
        inserted = [e for e in chain["versions"][3]["provenance"] if e["kind"] == "inserted"]
        self.assertEqual(inserted, [
            {"kind": "inserted", "edit_version": 3, "offset": 0},
            {"kind": "inserted", "edit_version": 3, "offset": 1},
        ])

    def test_undo_last_edit(self):
        chain = build_chain()
        before = versioning.trace_feature(chain, make_feature())
        chain = versioning.undo_last(chain)
        self.assertEqual(len(chain["versions"]), 3)
        self.assertEqual(chain["versions"][-1]["edit"]["kind"], "delete")
        after = versioning.trace_feature(chain, make_feature())
        truncated = build_chain()
        truncated["versions"] = truncated["versions"][:3]
        self.assertEqual(after, versioning.trace_feature(truncated, make_feature()))
        self.assertEqual(after["versions"], before["versions"][:3])
        with self.assertRaises(ValueError):
            versioning.undo_last(versioning.create_chain("ref-x", REFERENCE))


class TraceTest(unittest.TestCase):
    def test_cross_origin_definition_and_order_preserved(self):
        trace = versioning.trace_feature(build_chain(), make_feature())
        self.assertEqual(trace["feature_start"], 10)
        self.assertEqual(trace["feature_end"], 3)
        self.assertEqual(trace["feature_strand"], "-")
        self.assertTrue(trace["feature_wraps_origin"])
        v0 = trace["versions"][0]
        self.assertEqual(v0["positions"], (10, 11, 12, 1, 2, 3))
        self.assertNotEqual(v0["positions"], tuple(sorted(v0["positions"])))
        self.assertEqual([s["positions"] for s in v0["segments"]], [[10, 11, 12, 1, 2, 3]])
        self.assertEqual([s["source_positions"] for s in v0["segments"]], [[10, 11, 12, 1, 2, 3]])
        self.assertTrue(v0["segments"][0]["wraps_origin"])

    def test_intermediate_versions(self):
        trace = versioning.trace_feature(build_chain(), make_feature())
        self.assertEqual([v["version"] for v in trace["versions"]], [0, 1, 2, 3])
        v1 = trace["versions"][1]
        self.assertEqual(v1["status"], "projected")
        self.assertEqual(v1["positions"], (6, 7, 8, 9, 10, 11))
        self.assertEqual(len(v1["segments"]), 1)
        v2 = trace["versions"][2]
        self.assertEqual(v2["status"], "partially_deleted")
        self.assertEqual(v2["positions"], (6, 7, 8, 9))
        self.assertEqual(v2["deleted_source_positions"], (12, 1))
        self.assertEqual(v2["segments"][0]["source_positions"], [10, 11, 2, 3])
        v3 = trace["versions"][3]
        self.assertEqual(v3["status"], "partially_deleted")
        self.assertEqual(v3["positions"], (6, 7, 8, 9))
        self.assertEqual(v3["insertion_impact"]["location"], "boundary")
        self.assertEqual(v3["insertion_impact"]["edit"]["after_position"], 9)

    def test_internal_insertion_impact(self):
        chain = versioning.create_chain("ref-x", REFERENCE)
        chain = versioning.insert(chain, 11, "GG")
        trace = versioning.trace_feature(chain, make_feature())
        self.assertEqual(trace["versions"][1]["insertion_impact"]["location"], "internal")

    def test_export_intermediate_mapping(self):
        chain = build_chain()
        partial = versioning.trace_feature(chain, make_feature(), through_version=2)
        self.assertEqual([v["version"] for v in partial["versions"]], [0, 1, 2])
        full = versioning.trace_feature(chain, make_feature())
        self.assertEqual(partial["versions"], full["versions"][:3])

    def test_full_deletion_marked(self):
        chain = versioning.create_chain("ref-x", REFERENCE)
        chain = versioning.delete(chain, 10, 12)
        chain = versioning.delete(chain, 1, 3)
        trace = versioning.trace_feature(chain, make_feature())
        v2 = trace["versions"][2]
        self.assertEqual(v2["status"], "deleted")
        self.assertEqual(v2["positions"], ())
        self.assertEqual(v2["segments"], [])
        self.assertEqual(v2["deleted_source_positions"], (10, 11, 12, 1, 2, 3))
        v1 = trace["versions"][1]
        self.assertEqual(v1["status"], "partially_deleted")
        self.assertEqual(v1["positions"], (1, 2, 3))


class InsertedMemberPersistenceTest(unittest.TestCase):
    def _build(self):
        chain = versioning.create_chain("ref-x", REFERENCE, chain_id="chain-i")
        chain = versioning.rotate(chain, 4)          # v1
        chain = versioning.delete(chain, 8, 8)       # v2: origin position 12 deleted
        chain = versioning.insert(chain, 6, "GG")    # v3: internal, between origins 10 and 11
        chain = versioning.rotate(chain, 5)          # v4: rotate again
        chain = versioning.delete(chain, 2, 3)       # v5: deletes the two inserted bases
        return chain

    def test_inserted_bases_are_members_and_survive_rotation(self):
        trace = versioning.trace_feature(self._build(), make_feature())
        v3 = trace["versions"][3]
        self.assertEqual(v3["insertion_impact"]["location"], "internal")
        self.assertTrue(v3["insertion_impact"]["attributed"])
        self.assertEqual(v3["positions"], (6, 7, 8, 9, 10, 11, 12))
        inserted_sources = [
            source for source in v3["sources"] if source["kind"] == "inserted"
        ]
        self.assertEqual(inserted_sources, [
            {"kind": "inserted", "edit_version": 3, "offset": 0},
            {"kind": "inserted", "edit_version": 3, "offset": 1},
        ])
        self.assertEqual(
            [s["kind"] for s in v3["segments"][0]["sources"]],
            ["origin", "inserted", "inserted", "origin", "origin", "origin", "origin"],
        )
        self.assertEqual(v3["segments"][0]["source_positions"], [10, None, None, 11, 1, 2, 3])
        # the feature direction order is preserved across the insertion point
        self.assertEqual(
            [s.get("position") for s in v3["sources"]],
            [10, None, None, 11, 1, 2, 3],
        )

        v4 = trace["versions"][4]
        self.assertEqual(v4["positions"], (1, 2, 3, 4, 5, 6, 7))
        self.assertEqual(
            [s["kind"] for s in v4["sources"]],
            ["origin", "inserted", "inserted", "origin", "origin", "origin", "origin"],
        )
        inserted_after_rotation = [
            source for source in v4["sources"] if source["kind"] == "inserted"
        ]
        self.assertEqual(
            {s["edit_version"] for s in inserted_after_rotation}, {3}
        )

    def test_deleting_inserted_bases_reports_inserted_deleted_sources(self):
        trace = versioning.trace_feature(self._build(), make_feature())
        v5 = trace["versions"][5]
        self.assertEqual(v5["positions"], (1, 2, 3, 4, 5))
        self.assertEqual(v5["status"], "partially_deleted")
        self.assertEqual(v5["sources"][0], {"kind": "origin", "position": 10})
        self.assertEqual(v5["deleted_sources"], (
            {"kind": "inserted", "edit_version": 3, "offset": 0},
            {"kind": "inserted", "edit_version": 3, "offset": 1},
            {"kind": "origin", "position": 12},
        ))
        self.assertEqual(v5["deleted_source_positions"], (12,))
        self.assertEqual(tuple(v5["segments"][0]["sources"]), v5["sources"])

    def test_boundary_insertion_is_not_attributed(self):
        chain = versioning.create_chain("ref-x", REFERENCE, chain_id="chain-b")
        chain = versioning.rotate(chain, 4)
        chain = versioning.delete(chain, 8, 9)       # origins 12, 1 deleted
        chain = versioning.insert(chain, 9, "GG")    # boundary: only one side in feature
        trace = versioning.trace_feature(chain, make_feature())
        impact = trace["versions"][3]["insertion_impact"]
        self.assertEqual(impact["location"], "boundary")
        self.assertFalse(impact["attributed"])
        self.assertEqual(impact["inserted_sources"], ())
        v3 = trace["versions"][3]
        self.assertEqual(v3["positions"], (6, 7, 8, 9))
        self.assertTrue(
            all(source["kind"] == "origin" for source in v3["sources"])
        )
        # a later rotation still cannot pull the boundary bases into the feature
        chain = versioning.rotate(chain, 3)
        trace = versioning.trace_feature(chain, make_feature())
        self.assertEqual(trace["versions"][4]["positions"], (3, 4, 5, 6))
        self.assertTrue(
            all(s["kind"] == "origin" for s in trace["versions"][4]["sources"])
        )

    def test_outside_insertion_is_not_attributed(self):
        chain = versioning.create_chain("ref-x", REFERENCE, chain_id="chain-o")
        chain = versioning.insert(chain, 5, "GG")
        trace = versioning.trace_feature(chain, make_feature())
        impact = trace["versions"][1]["insertion_impact"]
        self.assertEqual(impact["location"], "outside")
        self.assertFalse(impact["attributed"])
        self.assertEqual(trace["versions"][1]["positions"], (12, 13, 14, 1, 2, 3))

    def test_atom_maps_cover_inserted_bases(self):
        chain = self._build()
        self.assertTrue(versioning.verify_chain(chain)["consistent"])
        for index in range(1, len(chain["versions"])):
            self.assertEqual(
                versioning.direct_atom_map(chain, index),
                versioning.composed_atom_map(chain, index),
            )
        atoms = versioning.composed_atom_map(chain, 5)
        self.assertNotIn(("inserted", 3, 0), atoms)
        self.assertNotIn(("inserted", 3, 1), atoms)
        atoms_v4 = versioning.composed_atom_map(chain, 4)
        self.assertEqual(atoms_v4[("inserted", 3, 0)], 2)
        self.assertEqual(atoms_v4[("inserted", 3, 1)], 3)
        atoms_v3 = versioning.composed_atom_map(chain, 3)
        self.assertEqual(atoms_v3[("inserted", 3, 0)], 7)
        self.assertEqual(atoms_v3[("inserted", 3, 1)], 8)

    def test_undo_drops_attributed_insertion_consistently(self):
        chain = self._build()
        full = versioning.trace_feature(chain, make_feature())
        chain = versioning.undo_last(chain)
        trace = versioning.trace_feature(chain, make_feature())
        self.assertEqual(trace["versions"], full["versions"][:5])
        self.assertTrue(versioning.verify_chain(chain)["consistent"])
        chain = versioning.undo_last(chain)
        trace = versioning.trace_feature(chain, make_feature())
        # undoing the rotation leaves the attributed insertion in place
        v3 = trace["versions"][3]
        self.assertEqual(
            [s["kind"] for s in v3["sources"]],
            ["origin", "inserted", "inserted", "origin", "origin", "origin", "origin"],
        )
        chain = versioning.undo_last(chain)              # undo the internal insertion
        trace = versioning.trace_feature(chain, make_feature())
        latest = trace["versions"][-1]
        self.assertTrue(
            all(s["kind"] == "origin" for s in latest["sources"])
        )
        self.assertEqual(latest["positions"], (6, 7, 8, 9, 10))


class VerifyTest(unittest.TestCase):
    def test_direct_matches_composed(self):
        chain = build_chain()
        report = versioning.verify_chain(chain)
        self.assertTrue(report["consistent"], report["mismatches"])
        self.assertEqual(report["checked_versions"], 3)
        for index in range(1, 4):
            self.assertEqual(
                versioning.direct_position_map(chain, index),
                versioning.composed_position_map(chain, index),
            )

    def test_verify_detects_tampering(self):
        chain = build_chain()
        chain["versions"][2]["provenance"][0] = {"kind": "carried", "position": 5}
        self.assertFalse(versioning.verify_chain(chain)["consistent"])

    def test_verify_detects_dropped_inserted_provenance_after_rotation(self):
        chain = versioning.create_chain("ref-x", REFERENCE, chain_id="chain-v")
        chain = versioning.insert(chain, 6, "G")
        chain = versioning.rotate(chain, 3)
        # the inserted base sat at position 7 of v1; rerouting its carried record
        # makes the composed atom map lose the inserted atom after rotation
        rotated_version = chain["versions"][2]
        self.assertEqual(rotated_version["provenance"][3], {"kind": "carried", "position": 7})
        rotated_version["provenance"][3] = {"kind": "carried", "position": 11}
        report = versioning.verify_chain(chain)
        self.assertFalse(report["consistent"])
        self.assertTrue(
            any(m.get("atom", {}).get("kind") == "inserted" for m in report["mismatches"])
        )


class PersistenceTest(unittest.TestCase):
    def test_save_reopen_preserves_trace_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = ProjectStore(Path(tmp) / "proj")
            store.initialize()
            reference = store.add_reference(REFERENCE, reference_id="ref-x")
            store.add_feature("ref-x", name="wrap", start=10, end=3, strand="-")
            chain = store.create_chain("ref-x", "chain-t")
            chain = store.update_chain("chain-t", "rotate", offset=4)
            chain = store.update_chain("chain-t", "delete", start=8, end=9)
            chain = store.update_chain("chain-t", "insert", after_position=9, bases="GG")

            reopened = ProjectStore(Path(tmp) / "proj")
            loaded = reopened.get_chain("chain-t")
            self.assertEqual(loaded, chain)
            inserted = [e for e in loaded["versions"][3]["provenance"] if e["kind"] == "inserted"]
            self.assertEqual(len(inserted), 2)
            self.assertEqual({e["edit_version"] for e in inserted}, {3})

            feature = reopened.get_features("ref-x")[0]
            self.assertEqual(
                versioning.trace_feature(loaded, feature),
                versioning.trace_feature(chain, feature),
            )
            self.assertTrue(versioning.verify_chain(loaded)["consistent"])

            undone = reopened.update_chain("chain-t", "undo")
            self.assertEqual(len(undone["versions"]), 3)
            self.assertEqual(reopened.get_chain("chain-t")["versions"], undone["versions"])

    def test_save_reopen_preserves_inserted_member_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = ProjectStore(Path(tmp) / "proj")
            store.initialize()
            store.add_reference(REFERENCE, reference_id="ref-x")
            store.add_feature("ref-x", name="wrap", start=10, end=3, strand="-")
            store.create_chain("ref-x", "chain-i")
            store.update_chain("chain-i", "rotate", offset=4)
            store.update_chain("chain-i", "delete", start=8, end=8)
            store.update_chain("chain-i", "insert", after_position=6, bases="GG")
            store.update_chain("chain-i", "rotate", offset=5)
            store.update_chain("chain-i", "delete", start=2, end=3)

            reopened = ProjectStore(Path(tmp) / "proj")
            feature = reopened.get_features("ref-x")[0]
            chain = reopened.get_chain("chain-i")
            trace = versioning.trace_feature(chain, feature)
            latest = trace["versions"][-1]
            self.assertEqual(latest["status"], "partially_deleted")
            self.assertEqual(
                [source["kind"] for source in latest["sources"]],
                ["origin"] * 5,
            )
            self.assertEqual(latest["deleted_sources"], (
                {"kind": "inserted", "edit_version": 3, "offset": 0},
                {"kind": "inserted", "edit_version": 3, "offset": 1},
                {"kind": "origin", "position": 12},
            ))
            self.assertTrue(versioning.verify_chain(chain)["consistent"])

            undone = reopened.update_chain("chain-i", "undo")
            trace = versioning.trace_feature(undone, feature)
            restored = trace["versions"][-1]
            self.assertEqual(restored["positions"], (1, 2, 3, 4, 5, 6, 7))
            self.assertEqual(
                [source["kind"] for source in restored["sources"]],
                ["origin", "inserted", "inserted", "origin", "origin", "origin", "origin"],
            )
            self.assertEqual(
                [source["edit_version"] for source in restored["sources"] if source["kind"] == "inserted"],
                [3, 3],
            )


if __name__ == "__main__":
    unittest.main()
