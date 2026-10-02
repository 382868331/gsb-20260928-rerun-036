from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioseq_platform.features import Feature, project_feature
from bioseq_platform.models import Alignment, AlignmentParameters
from bioseq_platform.pipeline import (
    project_run_record,
    run_from_files,
    select_alignment,
    selection_projection_diff,
)
from bioseq_platform.storage import ProjectStore


def rotate(sequence: str, offset: int) -> str:
    return sequence[offset:] + sequence[:offset]


def make_alignment(reference_coordinates, query_coordinates, variants=(), strand="forward", rotation_offset=0):
    return Alignment(
        score=0,
        reference_start=0,
        reference_end=0,
        query_start=0,
        query_end=0,
        aligned_reference="",
        aligned_query="",
        cigar="",
        strand=strand,
        rotation_offset=rotation_offset,
        variants=tuple(variants),
        reference_coordinates=tuple(reference_coordinates),
        query_coordinates=tuple(query_coordinates),
    )


def make_feature(start, end, strand="+", feature_id="feat-x"):
    return Feature(id=feature_id, reference_id="ref-x", name=feature_id, start=start, end=end, strand=strand)


class BoundaryRuleTest(unittest.TestCase):
    def test_boundary_insertion_not_internal_in_either_neighbor(self) -> None:
        feature_a = make_feature(2, 4, feature_id="feat-a")
        feature_b = make_feature(5, 7, feature_id="feat-b")
        insertion = {"type": "insertion", "after_position": 4, "inserted_bases": "G"}
        alignment = make_alignment([2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6], [insertion])
        impact_a = project_feature(feature_a, alignment, 8)["impacts"]
        impact_b = project_feature(feature_b, alignment, 8)["impacts"]
        self.assertEqual([item["location"] for item in impact_a], ["boundary"])
        self.assertEqual([item["location"] for item in impact_b], ["boundary"])

    def test_strictly_internal_insertion(self) -> None:
        feature = make_feature(2, 5)
        insertion = {"type": "insertion", "after_position": 3, "inserted_bases": "G"}
        alignment = make_alignment([2, 3, 4, 5], [1, 2, 3, 4], [insertion])
        impacts = project_feature(feature, alignment, 8)["impacts"]
        self.assertEqual([item["location"] for item in impacts], ["internal"])

    def test_wrap_feature_boundary_insertion_at_origin(self) -> None:
        feature = make_feature(7, 2)
        insertion = {"type": "insertion", "after_position": 2, "inserted_bases": "G"}
        alignment = make_alignment([7, 8, 1, 2], [1, 2, 3, 4], [insertion])
        impacts = project_feature(feature, alignment, 8)["impacts"]
        self.assertEqual([item["location"] for item in impacts], ["boundary"])
        internal = {"type": "insertion", "after_position": 8, "inserted_bases": "G"}
        alignment = make_alignment([7, 8, 1, 2], [1, 2, 3, 4], [internal])
        impacts = project_feature(feature, alignment, 8)["impacts"]
        self.assertEqual([item["location"] for item in impacts], ["internal"])

    def test_deletion_covering_whole_feature_has_no_range(self) -> None:
        feature = make_feature(3, 5)
        deletion = {"type": "deletion", "start": 3, "end": 5, "reference_bases": "CGT"}
        alignment = make_alignment([1, 2, 3, 4, 5, 6], [1, 2, 0, 0, 0, 3], [deletion])
        projection = project_feature(feature, alignment, 8)
        self.assertEqual(projection["status"], "deleted")
        self.assertIsNone(projection["query_start"])
        self.assertIsNone(projection["query_end"])
        self.assertTrue(projection["impacts"][0]["covers_feature"])

    def test_wrap_deletion_positions_cross_origin(self) -> None:
        feature = make_feature(8, 1)
        deletion = {"type": "deletion", "start": 7, "end": 2, "reference_bases": "TGCA"}
        alignment = make_alignment([7, 8, 1, 2], [0, 0, 0, 0], [deletion])
        projection = project_feature(feature, alignment, 8)
        self.assertEqual(projection["status"], "deleted")
        self.assertEqual(projection["impacts"][0]["positions_in_feature"], [1, 8])

    def test_partial_deletion_overlapping_feature_edge(self) -> None:
        feature = make_feature(3, 5)
        deletion = {"type": "deletion", "start": 5, "end": 7, "reference_bases": "TGC"}
        alignment = make_alignment([3, 4, 5, 6, 7], [1, 2, 0, 0, 0], [deletion])
        projection = project_feature(feature, alignment, 8)
        self.assertEqual(projection["status"], "projected")
        self.assertEqual(projection["impacts"][0]["location"], "partial")
        self.assertFalse(projection["impacts"][0]["covers_feature"])


class FeatureProjectionIntegrationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.project = Path(self.temporary.name) / "project"
        self.reference = self.project / "reference.fasta"
        self.sequences = self.project / "sequences.fasta"
        self.reference.parent.mkdir()
        self.store = ProjectStore(self.project)
        self.store.initialize()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _run(self, reference: str, query: str, feature_specs=(), select=None, mode="global", reference_id="ref"):
        self.reference.write_text(">" + reference_id + "\n" + reference + "\n", encoding="utf-8")
        self.sequences.write_text(">q1\n" + query + "\n", encoding="utf-8")
        parameters = AlignmentParameters(mode=mode, topology="circular")
        run_id, manifest, results = run_from_files(
            self.store,
            reference_path=str(self.reference),
            sequences_path=str(self.sequences),
            sequence_type="dna",
            parameters=parameters,
        )
        reference_id = manifest["reference_id"]
        for name, start, end, strand in feature_specs:
            self.store.add_feature(reference_id, name=name, start=start, end=end, strand=strand)
        if select is not None:
            select_alignment(self.store, run_id, "q1", select)
            results = self.store.get_run_results(run_id)
        return run_id, reference_id, results

    def test_wrap_feature_projects_as_single_query_range(self) -> None:
        run_id, _, results = self._run("AACGTTGC", "TGCAACGT", [("wrap", 7, 3, "+")])
        forward = next(item for item in results[0].alignments if item.strand == "forward")
        self.assertTrue(forward.wraps_reference)
        index = list(results[0].alignments).index(forward)
        projection = project_run_record(self.store, run_id, "q1", index)
        feature = projection["projections"][0]["features"][0]
        self.assertEqual(feature["status"], "projected")
        self.assertEqual(feature["query_start"], 2)
        self.assertEqual(feature["query_end"], 6)
        self.assertEqual(feature["reference_covered"], (1, 2, 3, 7, 8))

    def test_reverse_strand_projection_uses_original_query_coordinates(self) -> None:
        run_id, _, results = self._run("AACGTTGC", "GCAACGTT", [("gene", 2, 4, "+")])
        reverse = [item for item in results[0].alignments if item.strand == "reverse"]
        self.assertEqual(len(reverse), 1)
        index = list(results[0].alignments).index(reverse[0])
        projection = project_run_record(self.store, run_id, "q1", index)
        feature = projection["projections"][0]["features"][0]
        self.assertEqual(feature["query_strand"], "reverse")
        self.assertEqual(feature["query_start"], 5)
        self.assertEqual(feature["query_end"], 7)
        self.assertFalse(feature["query_wraps_origin"])

    def test_tied_rotations_keep_differing_projections(self) -> None:
        run_id, _, results = self._run("ACGTACGT", "ACGT", [("box", 1, 4, "+")], mode="local")
        forward = [item for item in results[0].alignments if item.strand == "forward"]
        self.assertGreaterEqual(len(forward), 2)
        projection = project_run_record(self.store, run_id, "q1")
        statuses = [
            entry["features"][0]["status"]
            for entry in projection["projections"]
            if entry["strand"] == "forward"
        ]
        self.assertIn("projected", statuses)
        self.assertIn("not_covered", statuses)

    def test_selection_persists_and_preview_reports_recompute(self) -> None:
        run_id, _, results = self._run("ACGTACGT", "ACGT", [("box", 1, 4, "+")])
        forward = [item for item in results[0].alignments if item.strand == "forward"]
        first = list(results[0].alignments).index(forward[0])
        second = list(results[0].alignments).index(forward[1])
        select_alignment(self.store, run_id, "q1", first)
        preview = selection_projection_diff(self.store, run_id, "q1", second)
        self.assertEqual(preview["previous_alignment_index"], first)
        self.assertEqual(len(preview["features_to_recompute"]), 1)
        select_alignment(self.store, run_id, "q1", second)
        reopened = ProjectStore(self.project)
        stored = reopened.get_run_results(run_id)[0]
        self.assertEqual(stored.selected_alignment_index, second)
        projection = project_run_record(reopened, run_id, "q1")
        selected = [entry for entry in projection["projections"] if entry["selected"]]
        self.assertEqual([entry["alignment_index"] for entry in selected], [second])

    def test_rotation_of_reference_and_query_preserves_feature_impact(self) -> None:
        reference = "AACGTTGC"
        query = "TGTAACGT"  # rotation 5 of reference with substitution at reference position 8
        run_id, _, results = self._run(reference, query, [("wrap", 7, 3, "+")])
        base_projection = project_run_record(self.store, run_id, "q1")
        base_entry = next(
            entry for entry in base_projection["projections"]
            if entry["strand"] == "forward" and entry["features"][0]["impacts"]
        )
        base_feature = base_entry["features"][0]
        base_impact = base_feature["impacts"][0]
        self.assertEqual(base_impact["type"], "substitution")
        self.assertEqual(base_impact["position"], 8)

        rotated_reference = rotate(reference, 3)
        rotated_query = rotate(query, 2)
        run_id2, _, results2 = self._run(
            rotated_reference, rotated_query, [("wrap", 4, 8, "+")], reference_id="ref-rot"
        )
        rotated_projection = project_run_record(self.store, run_id2, "q1")
        rotated_entry = next(
            entry for entry in rotated_projection["projections"]
            if entry["strand"] == "forward" and entry["features"][0]["impacts"]
        )
        rotated_feature = rotated_entry["features"][0]
        rotated_impact = rotated_feature["impacts"][0]
        self.assertEqual(rotated_impact["type"], "substitution")
        self.assertEqual(rotated_impact["position"], 5)  # 8 shifted by rotation 3
        self.assertEqual(rotated_feature["status"], base_feature["status"])
        expected_query_start = ((base_feature["query_start"] - 1 - 2) % 8) + 1
        self.assertEqual(rotated_feature["query_start"], expected_query_start)
        self.assertTrue(rotated_feature["query_wraps_origin"])


if __name__ == "__main__":
    unittest.main()
