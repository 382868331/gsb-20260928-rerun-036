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


class InsertedFragmentTrackingTest(unittest.TestCase):
    def build_full_chain(self):
        chain = versioning.create_chain("ref-x", REFERENCE, chain_id="chain-i")
        chain = versioning.rotate(chain, 4)                    # v1
        chain = versioning.delete(chain, 8, 9)                 # v2: partial deletion
        chain = versioning.insert(chain, 7, "GG")              # v3: internal insertion
        chain = versioning.rotate(chain, 3)                    # v4: rotate again
        chain = versioning.delete(chain, 5, 6)                 # v5: delete inserted bases
        return chain

    def test_internal_insertion_joins_feature_and_rotates_along(self):
        trace = versioning.trace_feature(self.build_full_chain(), make_feature())
        v3 = trace["versions"][3]
        self.assertEqual(v3["insertion_impact"]["location"], "internal")
        self.assertEqual(v3["positions"], (6, 7, 8, 9, 10, 11))
        self.assertEqual(v3["inserted_positions"], (8, 9))
        provenance = [p for s in v3["segments"] for p in s["provenance"]]
        self.assertEqual(provenance[2:4], [
            {"kind": "inserted", "edit_version": 3, "offset": 0},
            {"kind": "inserted", "edit_version": 3, "offset": 1},
        ])
        self.assertEqual(provenance[0], {"kind": "origin", "position": 10})
        self.assertIsNone(v3["segments"][0]["source_positions"][2])

        v4 = trace["versions"][4]
        self.assertEqual(v4["status"], "partially_deleted")
        self.assertEqual(v4["deleted_source_positions"], (12, 1))
        self.assertEqual(v4["deleted_insertions"], ())
        self.assertEqual(v4["positions"], (3, 4, 5, 6, 7, 8))
        self.assertEqual(v4["inserted_positions"], (5, 6))
        provenance = [p for s in v4["segments"] for p in s["provenance"]]
        self.assertEqual(provenance[2], {"kind": "inserted", "edit_version": 3, "offset": 0})
        self.assertEqual([s["positions"] for s in v4["segments"]], [[3, 4, 5, 6, 7, 8]])

    def test_later_deletion_of_inserted_bases_is_traced(self):
        trace = versioning.trace_feature(self.build_full_chain(), make_feature())
        v5 = trace["versions"][5]
        self.assertEqual(v5["status"], "partially_deleted")
        self.assertEqual(v5["positions"], (3, 4, 5, 6))
        self.assertEqual(v5["inserted_positions"], ())
        self.assertEqual(v5["deleted_insertions"], (
            {"kind": "inserted", "edit_version": 3, "offset": 0},
            {"kind": "inserted", "edit_version": 3, "offset": 1},
        ))
        self.assertEqual(v5["deleted_source_positions"], (12, 1))

    def test_boundary_insertion_not_attributed(self):
        chain = versioning.create_chain("ref-x", REFERENCE)
        chain = versioning.rotate(chain, 4)
        chain = versioning.delete(chain, 8, 9)
        chain = versioning.insert(chain, 9, "GG")              # boundary
        chain = versioning.rotate(chain, 2)
        trace = versioning.trace_feature(chain, make_feature())
        v3 = trace["versions"][3]
        self.assertEqual(v3["insertion_impact"]["location"], "boundary")
        self.assertEqual(v3["inserted_positions"], ())
        self.assertEqual(v3["positions"], (6, 7, 8, 9))
        v4 = trace["versions"][4]
        self.assertEqual(v4["inserted_positions"], ())
        self.assertEqual(v4["positions"], (4, 5, 6, 7))
        provenance = [p for s in v4["segments"] for p in s["provenance"]]
        self.assertTrue(all(p["kind"] == "origin" for p in provenance))

    def test_verify_covers_inserted_provenance(self):
        chain = self.build_full_chain()
        report = versioning.verify_chain(chain)
        self.assertTrue(report["consistent"], report)
        chain["versions"][3]["provenance"][7] = {"kind": "inserted", "edit_version": 3, "offset": 1}
        report = versioning.verify_chain(chain)
        self.assertFalse(report["consistent"])
        self.assertTrue(report["insertion_mismatches"])

    def test_save_reopen_and_undo_preserve_insertion_trace(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = ProjectStore(Path(tmp) / "proj")
            store.initialize()
            store.add_reference(REFERENCE, reference_id="ref-x")
            store.add_feature("ref-x", name="wrap", start=10, end=3, strand="-")
            store.create_chain("ref-x", "chain-i")
            store.update_chain("chain-i", "rotate", offset=4)
            store.update_chain("chain-i", "delete", start=8, end=9)
            store.update_chain("chain-i", "insert", after_position=7, bases="GG")
            store.update_chain("chain-i", "rotate", offset=3)
            chain = store.update_chain("chain-i", "delete", start=5, end=6)

            reopened = ProjectStore(Path(tmp) / "proj")
            loaded = reopened.get_chain("chain-i")
            feature = reopened.get_features("ref-x")[0]
            self.assertEqual(
                versioning.trace_feature(loaded, feature),
                versioning.trace_feature(chain, feature),
            )
            trace = versioning.trace_feature(loaded, feature)
            self.assertEqual(trace["versions"][4]["inserted_positions"], (5, 6))
            self.assertTrue(versioning.verify_chain(loaded)["consistent"])

            undone = reopened.update_chain("chain-i", "undo")
            self.assertEqual(len(undone["versions"]), 5)
            restored = versioning.trace_feature(reopened.get_chain("chain-i"), feature)
            self.assertEqual(restored["versions"], trace["versions"][:5])
            self.assertEqual(restored["versions"][4]["inserted_positions"], (5, 6))
            self.assertEqual(restored["versions"][4]["deleted_insertions"], ())


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


if __name__ == "__main__":
    unittest.main()
