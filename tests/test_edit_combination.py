from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioseq_platform import combination, versioning
from bioseq_platform.cli import main as cli_main
from bioseq_platform.features import Feature
from bioseq_platform.storage import ProjectStore

REFERENCE = "ACGTACGTACGT"  # length 12


def make_feature(start=10, end=3, strand="-"):
    return Feature(
        id="feat-wrap", reference_id="ref-x", name="wrap",
        start=start, end=end, strand=strand,
    )


def base_chain():
    return versioning.create_chain("ref-x", REFERENCE, chain_id="chain-t")


class StudyTest(unittest.TestCase):
    def test_shared_deletion_ranges_merge(self):
        candidates = [
            {"id": "d1", "kind": "delete", "start": 3, "end": 5},
            {"id": "d2", "kind": "delete", "start": 4, "end": 7},
        ]
        study = combination.study_combinations(base_chain(), [], None, [], candidates)
        self.assertEqual(study["conflicts"], [])
        best = study["best"]
        self.assertEqual(best["candidate_ids"], ["d1", "d2"])
        self.assertEqual(best["deleted_sources"], 5)  # union 3..7, not 3+4
        result = combination.confirm_combination(
            base_chain(), [], None, [], candidates, ["d1", "d2"],
            new_chain_id="chain-shared",
        )
        self.assertEqual(result["final_sequence"], "ACTACGT")
        self.assertTrue(result["verification"]["consistent"])
        self.assertEqual(result["new_chain_id"], "chain-shared")
        self.assertEqual(result["base_chain_id"], "chain-t")

    def test_consecutive_insertions_same_anchor(self):
        candidates = [
            {"id": "c1", "kind": "insert", "after_position": 5, "bases": "AA"},
            {"id": "c2", "kind": "insert", "after_position": 5, "bases": "CC"},
        ]
        study = combination.study_combinations(base_chain(), [], None, [], candidates)
        self.assertEqual(study["conflicts"], [])
        result = combination.confirm_combination(
            base_chain(), [], None, [], candidates, ["c1", "c2"],
        )
        self.assertEqual(result["final_sequence"], "ACGTAAACCCGTACGT")
        self.assertTrue(result["verification"]["consistent"])

    def test_delete_insert_conflict_reports_actual_source(self):
        candidates = [
            {"id": "d1", "kind": "delete", "start": 4, "end": 6},
            {"id": "i1", "kind": "insert", "after_position": 5, "bases": "GG"},
        ]
        study = combination.study_combinations(base_chain(), [], None, [], candidates)
        self.assertEqual(len(study["conflicts"]), 1)
        conflict = study["conflicts"][0]
        self.assertEqual(conflict["reason"], "insert_anchor_deleted")
        self.assertEqual(conflict["affected_sources"], [{"kind": "origin", "position": 5}])
        for combo in study["combinations"]:
            self.assertNotEqual(set(combo["candidate_ids"]), {"d1", "i1"})
        with self.assertRaises(ValueError):
            combination.confirm_combination(
                base_chain(), [], None, [], candidates, ["d1", "i1"]
            )

    def test_delete_rotate_conflict_reports_actual_source(self):
        candidates = [
            {"id": "d1", "kind": "delete", "start": 4, "end": 6},
            {"id": "r1", "kind": "rotate", "offset": 3},
        ]
        study = combination.study_combinations(base_chain(), [], None, [], candidates)
        self.assertEqual(len(study["conflicts"]), 1)
        conflict = study["conflicts"][0]
        self.assertEqual(conflict["reason"], "rotation_start_deleted")
        self.assertEqual(conflict["affected_sources"], [{"kind": "origin", "position": 4}])

    def test_same_bases_from_different_versions_are_distinct_sources(self):
        chain = base_chain()
        chain = versioning.insert(chain, 2, "GG")   # v1 inserted atoms at positions 3,4
        chain = versioning.insert(chain, 6, "GG")   # v2 inserted atoms at positions 7,8
        candidates = [
            {"id": "d1", "kind": "delete", "start": 3, "end": 3},
            {"id": "d2", "kind": "delete", "start": 7, "end": 8},
        ]
        study = combination.study_combinations(chain, [], None, [], candidates)
        reports = {report["id"]: report for report in study["candidates"]}
        self.assertEqual(
            reports["d1"]["affected_sources"],
            [{"kind": "inserted", "edit_version": 1, "offset": 0}],
        )
        self.assertEqual(
            reports["d2"]["affected_sources"],
            [
                {"kind": "inserted", "edit_version": 2, "offset": 0},
                {"kind": "inserted", "edit_version": 2, "offset": 1},
            ],
        )
        result = combination.confirm_combination(
            chain, [], None, [], candidates, ["d1", "d2"]
        )
        self.assertTrue(result["verification"]["consistent"])
        self.assertEqual(result["final_sequence"], "ACGGTACGTACGT"[0:0] + "ACGGTACGTACGT")

    def test_kept_feature_violation_uses_feature_direction_order(self):
        feature = make_feature()
        candidates = [{"id": "d1", "kind": "delete", "start": 11, "end": 12}]
        study = combination.study_combinations(
            base_chain(), [feature], None, [], candidates
        )
        report = study["candidates"][0]
        self.assertEqual(len(report["violations"]), 1)
        violation = report["violations"][0]
        self.assertEqual(violation["feature_id"], "feat-wrap")
        self.assertEqual(violation["feature_strand"], "-")
        self.assertTrue(violation["feature_wraps_origin"])
        self.assertEqual(
            violation["affected_sources"],
            [{"kind": "origin", "position": 11}, {"kind": "origin", "position": 12}],
        )
        self.assertEqual(study["combinations"], [])
        self.assertIsNone(study["best"])

    def test_allowed_range_exempts_disruption(self):
        feature = make_feature()
        candidates = [{"id": "d1", "kind": "delete", "start": 11, "end": 12}]
        study = combination.study_combinations(
            base_chain(), [feature], None, [(11, 12)], candidates
        )
        self.assertEqual(study["candidates"][0]["violations"], [])
        self.assertEqual(study["best"]["candidate_ids"], ["d1"])

    def test_boundary_insert_on_kept_feature(self):
        feature = make_feature()
        candidates = [{"id": "i1", "kind": "insert", "after_position": 3, "bases": "GG"}]
        study = combination.study_combinations(
            base_chain(), [feature], None, [], candidates
        )
        violation = study["candidates"][0]["violations"][0]
        self.assertEqual(violation["location"], "boundary")
        allowed = combination.study_combinations(
            base_chain(), [feature], None, [(3, 3)], candidates
        )
        self.assertEqual(allowed["candidates"][0]["violations"], [])

    def test_minimal_change_ranking(self):
        candidates = [
            {"id": "d1", "kind": "delete", "start": 1, "end": 2},
            {"id": "d2", "kind": "delete", "start": 5, "end": 5},
            {"id": "i1", "kind": "insert", "after_position": 8, "bases": "A"},
        ]
        study = combination.study_combinations(base_chain(), [], None, [], candidates)
        best = study["best"]
        self.assertEqual(best["candidate_ids"], ["d1", "d2", "i1"])
        self.assertEqual(best["change_cost"], 4)
        two = [c for c in study["combinations"] if c["candidate_count"] == 2]
        self.assertEqual(two[0]["change_cost"], 2)
        self.assertEqual(two[-1]["change_cost"], 3)

    def test_confirm_builds_independent_chain_and_traces_features(self):
        feature = make_feature()
        chain = base_chain()
        candidates = [
            {"id": "d1", "kind": "delete", "start": 11, "end": 11},
            {"id": "r1", "kind": "rotate", "offset": 2},
        ]
        result = combination.confirm_combination(
            chain, [feature], None, [(11, 11)], candidates, ["d1", "r1"],
            new_chain_id="chain-branch",
        )
        self.assertEqual(result["final_sequence"], "GTACGTACTAC")
        self.assertTrue(result["verification"]["consistent"])
        kept = result["kept_features"][0]
        self.assertEqual(kept["feature_id"], "feat-wrap")
        self.assertEqual(kept["status"], "partially_deleted")
        member_positions = [s["position"] for s in kept["sources"] if s["kind"] == "origin"]
        self.assertEqual(member_positions, [10, 12, 1, 2, 3])
        # base chain untouched
        self.assertEqual(len(chain["versions"]), 1)
        trace = versioning.trace_feature(result["chain"], feature)
        self.assertEqual(trace["versions"][-1]["deleted_source_positions"], (11,))


class StoreFlowTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = ProjectStore(Path(self.tmp.name) / "proj")
        self.store.initialize()
        self.reference = self.store.add_reference(REFERENCE, reference_id="ref-x")
        self.feature = self.store.add_feature(
            "ref-x", name="wrap", start=10, end=3, strand="-"
        )
        self.chain = self.store.create_chain("ref-x", "chain-t")

    def tearDown(self):
        self.tmp.cleanup()

    def test_confirm_save_reload_verify(self):
        chain = self.store.get_chain("chain-t")
        features = self.store.get_features("ref-x")
        candidates = [
            {"id": "d1", "kind": "delete", "start": 11, "end": 12},
            {"id": "i1", "kind": "insert", "after_position": 2, "bases": "GG"},
        ]
        result = combination.confirm_combination(
            chain, features, None, [(11, 12), (2, 3)], candidates, ["d1", "i1"],
            new_chain_id="chain-branch",
        )
        self.store.save_chain(result["chain"])
        reloaded = self.store.get_chain("chain-branch")
        self.assertTrue(versioning.verify_chain(reloaded)["consistent"])
        trace = versioning.trace_feature(reloaded, self.feature)
        final = trace["versions"][-1]
        self.assertEqual(final["status"], "partially_deleted")
        self.assertEqual(final["deleted_source_positions"], (11, 12))
        self.assertEqual(
            [s["position"] for s in final["sources"] if s["kind"] == "origin"],
            [10, 1, 2, 3],
        )
        # original chain unchanged
        self.assertEqual(len(self.store.get_chain("chain-t")["versions"]), 1)

    def test_cli_study_and_confirm(self):
        project = str(Path(self.tmp.name) / "proj")
        candidates_path = Path(self.tmp.name) / "candidates.json"
        candidates_path.write_text(json.dumps([
            {"id": "d1", "kind": "delete", "start": 11, "end": 12},
            {"id": "i1", "kind": "insert", "after_position": 11, "bases": "GG"},
        ]))
        study_code = cli_main([
            "combo-study", project, "chain-t",
            "--candidates", str(candidates_path),
            "--keep", self.feature.id,
        ])
        self.assertEqual(study_code, 0)
        confirm_code = cli_main([
            "combo-confirm", project, "chain-t",
            "--candidates", str(candidates_path),
            "--select", "d1",
            "--keep", self.feature.id,
            "--allow", "11:12",
            "--new-chain-id", "chain-cli",
        ])
        self.assertEqual(confirm_code, 0)
        chain = self.store.get_chain("chain-cli")
        self.assertTrue(versioning.verify_chain(chain)["consistent"])
        self.assertEqual(chain["versions"][-1]["sequence"], "ACGTACGTAC")


if __name__ == "__main__":
    unittest.main()
