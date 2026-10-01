from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioseq_platform.models import AlignmentParameters
from bioseq_platform.pipeline import run_from_files, select_alignment
from bioseq_platform.storage import ProjectStore


def rotate(sequence: str, offset: int) -> str:
    return sequence[offset:] + sequence[:offset]


class CircularAlignmentTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.project = Path(self.temporary.name) / "project"
        self.reference = self.project / "reference.fasta"
        self.sequences = self.project / "sequences.fasta"
        self.reference.parent.mkdir()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _run(self, reference: str, content: str, parameters: AlignmentParameters):
        self.reference.write_text(f">ref\n{reference}\n", encoding="utf-8")
        self.sequences.write_text(content, encoding="utf-8")
        store = ProjectStore(self.project)
        return run_from_files(
            store,
            reference_path=self.reference,
            sequences_path=self.sequences,
            sequence_type="dna",
            parameters=parameters,
        )

    def test_wrap_match_is_single_segment_not_split_fragments(self) -> None:
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, results = self._run("ACGT", ">rot\nGTAC\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 8)
        forward = [item for item in result.alignments if item.strand == "forward"]
        self.assertEqual(len(forward), 1)
        alignment = forward[0]
        self.assertEqual(alignment.cigar, "4=")
        self.assertEqual(alignment.rotation_offset, 2)
        self.assertTrue(alignment.wraps_reference)
        self.assertEqual(alignment.reference_coordinates, (3, 4, 1, 2))
        self.assertEqual(alignment.aligned_reference, "GTAC")
        self.assertEqual(alignment.aligned_query, "GTAC")

    def test_rotation_of_record_start_gives_equivalent_result(self) -> None:
        reference = "AACGTTGC"
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, first = self._run(reference, f">q\n{rotate(reference, 3)}\n", parameters)
        _, _, second = self._run(reference, f">q\n{rotate(reference, 5)}\n", parameters)
        first_forward = [item for item in first[0].alignments if item.strand == "forward"]
        second_forward = [item for item in second[0].alignments if item.strand == "forward"]
        self.assertEqual(len(first_forward), 1)
        self.assertEqual(len(second_forward), 1)
        self.assertEqual(first[0].score, second[0].score)
        self.assertEqual(first[0].score, 2 * len(reference))
        self.assertEqual(first_forward[0].cigar, second_forward[0].cigar)
        self.assertEqual(first_forward[0].cigar, "8=")
        self.assertEqual(first_forward[0].rotation_offset, 3)
        self.assertEqual(second_forward[0].rotation_offset, 5)
        self.assertEqual(first_forward[0].reference_coordinates, (4, 5, 6, 7, 8, 1, 2, 3))
        self.assertEqual(second_forward[0].reference_coordinates, (6, 7, 8, 1, 2, 3, 4, 5))

    def test_repeats_keep_all_tied_rotation_starts(self) -> None:
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, results = self._run("AATAAT", ">rep\nAATAAT\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 12)
        forward = [item for item in result.alignments if item.strand == "forward"]
        self.assertEqual(len(forward), 2)
        self.assertEqual([item.rotation_offset for item in forward], [0, 3])
        self.assertEqual(forward[0].reference_coordinates, (1, 2, 3, 4, 5, 6))
        self.assertEqual(forward[1].reference_coordinates, (4, 5, 6, 1, 2, 3))
        self.assertTrue(forward[1].wraps_reference)
        self.assertFalse(forward[0].wraps_reference)

    def test_wrap_crossing_substitution_uses_original_coordinates(self) -> None:
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, results = self._run("ACGT", ">sub\nGGAC\n", parameters)
        forward = [item for item in results[0].alignments if item.strand == "forward"]
        self.assertEqual(len(forward), 1)
        alignment = forward[0]
        self.assertEqual(alignment.cigar, "1=1X2=")
        self.assertEqual(alignment.reference_coordinates, (3, 4, 1, 2))
        self.assertEqual(alignment.variants, ({
            "type": "substitution",
            "position": 4,
            "reference_base": "T",
            "alternative_base": "G",
        },))

    def test_wrap_crossing_insertion_uses_original_coordinates(self) -> None:
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, results = self._run("ACGT", ">ins\nGTTAC\n", parameters)
        forward = [item for item in results[0].alignments if item.strand == "forward"]
        self.assertEqual(len(forward), 2)
        self.assertEqual({item.cigar for item in forward}, {"2=1I2=", "1=1I3="})
        boundary = [item for item in forward if item.cigar == "2=1I2="]
        self.assertEqual(len(boundary), 1)
        alignment = boundary[0]
        self.assertEqual(alignment.reference_coordinates, (3, 4, 0, 1, 2))
        self.assertEqual(alignment.variants, ({
            "type": "insertion",
            "after_position": 4,
            "inserted_bases": "T",
        },))
        self.assertEqual(
            {item.variants[0]["after_position"] for item in forward},
            {3, 4},
        )

    def test_short_query_is_not_reused_to_inflate_score(self) -> None:
        parameters = AlignmentParameters(mode="local", topology="circular")
        _, _, results = self._run("AAAAAA", ">short\nAAA\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 6)
        for alignment in result.alignments:
            used = sorted(coord for coord in alignment.query_coordinates if coord)
            self.assertEqual(used, [1, 2, 3])

    def test_long_query_cannot_wrap_reference_twice(self) -> None:
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, results = self._run("AAA", ">long\nAAAAAA\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 0)
        for alignment in result.alignments:
            reference_used = [coord for coord in alignment.reference_coordinates if coord]
            self.assertEqual(len(reference_used), len(set(reference_used)))

    def test_linear_and_circular_interpretations_are_comparable(self) -> None:
        linear = AlignmentParameters(mode="global", topology="linear")
        circular = AlignmentParameters(mode="global", topology="circular")
        linear_id, _, linear_results = self._run("ACGT", ">rot\nGTAC\n", linear)
        circular_id, _, circular_results = self._run("ACGT", ">rot\nGTAC\n", circular)
        self.assertNotEqual(linear_id, circular_id)
        self.assertGreater(circular_results[0].score, linear_results[0].score)
        self.assertEqual(circular_results[0].score, 8)
        self.assertTrue(all(not item.wraps_reference for item in linear_results[0].alignments))

    def test_selected_candidate_survives_save_and_reload(self) -> None:
        parameters = AlignmentParameters(mode="global", topology="circular")
        run_id, _, results = self._run("AATAAT", ">rep\nAATAAT\n", parameters)
        store = ProjectStore(self.project)
        select_alignment(store, run_id, "rep", 1)
        reloaded = ProjectStore(self.project).get_run_results(run_id)
        original = results[0].alignments
        again = reloaded[0].alignments
        self.assertEqual(reloaded[0].selected_alignment_index, 1)
        self.assertEqual(len(original), len(again))
        self.assertEqual(
            [(item.strand, item.rotation_offset, item.wraps_reference, item.reference_coordinates) for item in original],
            [(item.strand, item.rotation_offset, item.wraps_reference, item.reference_coordinates) for item in again],
        )
        self.assertEqual(again[1].rotation_offset, 3)

    def test_select_alignment_validates_record_and_index(self) -> None:
        parameters = AlignmentParameters(mode="global", topology="circular")
        run_id, _, _ = self._run("AATAAT", ">rep\nAATAAT\n", parameters)
        store = ProjectStore(self.project)
        with self.assertRaises(KeyError):
            select_alignment(store, run_id, "missing", 0)
        with self.assertRaises(IndexError):
            select_alignment(store, run_id, "rep", 999)


if __name__ == "__main__":
    unittest.main()
