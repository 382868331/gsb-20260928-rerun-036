from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioseq_platform.models import AlignmentParameters
from bioseq_platform.pipeline import run_from_files
from bioseq_platform.storage import ProjectStore
from bioseq_platform.alignment import reverse_complement


def rotate(sequence: str, k: int) -> str:
    return sequence[k:] + sequence[:k]


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

    def test_circular_global_reports_rotation_and_wrapped_coordinates(self) -> None:
        reference = "ACGTTGCA"
        query = rotate(reference, 3)
        self.assertEqual(query, "TTGCAACG")
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, results = self._run(reference, f">rot\n{query}\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 16)
        alignment = result.alignments[0]
        self.assertEqual(alignment.strand, "forward")
        self.assertEqual(alignment.query_rotation, 5)
        self.assertEqual(alignment.cigar, "8=")
        self.assertEqual(alignment.aligned_query, reference)
        self.assertEqual(alignment.query_coordinates, (6, 7, 8, 1, 2, 3, 4, 5))
        self.assertEqual((alignment.query_start, alignment.query_end), (6, 5))
        self.assertEqual(alignment.variants, ())

    def test_circular_global_reports_reverse_strand_and_start(self) -> None:
        reference = "ACCGTACG"
        query = rotate(reverse_complement(reference), 3)
        self.assertEqual(query, "ACGGTCGT")
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, results = self._run(reference, f">rev\n{query}\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 16)
        self.assertEqual(result.alignment_count, 1)
        alignment = result.alignments[0]
        self.assertEqual(alignment.strand, "reverse")
        self.assertEqual(alignment.query_rotation, 3)
        self.assertEqual(alignment.aligned_query, reference)
        self.assertEqual(alignment.query_coordinates, (5, 4, 3, 2, 1, 8, 7, 6))
        self.assertEqual((alignment.query_start, alignment.query_end), (5, 6))

    def test_circular_repeats_keep_all_tied_starts_and_selection_survives_reload(self) -> None:
        parameters = AlignmentParameters(mode="global", topology="circular")
        run_id, _, results = self._run("AAAA", ">rep\nAAAA\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 8)
        self.assertEqual(result.alignment_count, 4)
        self.assertEqual([item.query_rotation for item in result.alignments], [0, 1, 2, 3])
        self.assertTrue(all(item.strand == "forward" for item in result.alignments))
        self.assertEqual(
            [item.query_coordinates for item in result.alignments],
            [(1, 2, 3, 4), (2, 3, 4, 1), (3, 4, 1, 2), (4, 1, 2, 3)],
        )

        store = ProjectStore(self.project)
        chosen = store.select_alignment(run_id, "rep", 3)
        self.assertEqual(chosen.selected_alignment, 3)
        reloaded = store.get_run_results(run_id)[0]
        self.assertEqual(reloaded.selected_alignment, 3)
        self.assertEqual(
            [(item.strand, item.query_rotation, item.query_coordinates) for item in reloaded.alignments],
            [(item.strand, item.query_rotation, item.query_coordinates) for item in result.alignments],
        )
        self.assertEqual(reloaded.alignments[2].query_coordinates, (3, 4, 1, 2))

    def test_circular_mode_never_reuses_sequence_to_inflate_score(self) -> None:
        parameters = AlignmentParameters(mode="local", topology="circular")
        _, _, results = self._run("ACAC", ">short\nAC\n", parameters)
        self.assertEqual(results[0].score, 4)
        for alignment in results[0].alignments:
            self.assertLessEqual(len(alignment.aligned_query.replace("-", "")), 2)

        global_parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, global_results = self._run("ACAC", ">short\nAC\n", global_parameters)
        self.assertEqual(global_results[0].score, 0)
        self.assertEqual(global_results[0].alignments[0].cigar, "2=2D")

    def test_rotation_of_record_start_yields_equivalent_results(self) -> None:
        reference = "ACCGTACG"
        mutated = "ACCATACG"
        rotated = rotate(mutated, 3)
        self.assertEqual(rotated, "ATACGACC")
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, first = self._run(reference, f">a\n{mutated}\n", parameters)
        _, _, second = self._run(reference, f">b\n{rotated}\n", parameters)
        first_result, second_result = first[0], second[0]
        self.assertEqual(first_result.score, 13)
        self.assertEqual(second_result.score, 13)
        self.assertEqual(first_result.alignment_count, second_result.alignment_count)
        first_alignment = first_result.alignments[0]
        second_alignment = second_result.alignments[0]
        self.assertEqual(first_alignment.query_rotation, 0)
        self.assertEqual(second_alignment.query_rotation, 5)
        self.assertEqual(first_alignment.variants[0]["position"], 4)
        self.assertEqual(second_alignment.variants[0]["position"], 4)
        self.assertEqual(first_alignment.variants[0]["query_position"], 4)
        self.assertEqual(second_alignment.variants[0]["query_position"], 1)
        self.assertEqual(first_alignment.cigar, second_alignment.cigar)

    def test_rotation_invariance_with_repeats(self) -> None:
        reference = "AAA"
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, first = self._run(reference, ">a\nAAA\n", parameters)
        _, _, second = self._run(reference, ">b\nAAA\n", parameters)
        self.assertEqual(first[0].score, second[0].score)
        self.assertEqual(first[0].alignment_count, 3)
        self.assertEqual(second[0].alignment_count, 3)
        self.assertEqual(
            [item.query_coordinates for item in first[0].alignments],
            [item.query_coordinates for item in second[0].alignments],
        )

    def test_cross_boundary_variants_use_original_input_coordinates(self) -> None:
        reference = "ACCGTACG"
        inserted = reference + "T"
        inserted_query = rotate(inserted, 4)
        self.assertEqual(inserted_query, "TACGTACCG")
        deleted = "ACCGACG"
        deleted_query = rotate(deleted, 3)
        self.assertEqual(deleted_query, "GACGACC")
        parameters = AlignmentParameters(mode="global", topology="circular")
        _, _, results = self._run(
            reference,
            f">ins\n{inserted_query}\n>del\n{deleted_query}\n",
            parameters,
        )
        insertion, deletion = results
        ins_alignment = next(item for item in insertion.alignments if item.query_rotation == 5)
        self.assertEqual(ins_alignment.query_rotation, 5)
        self.assertEqual(ins_alignment.cigar, "8=1I")
        self.assertEqual(ins_alignment.variants, ({
            "type": "insertion",
            "after_position": 8,
            "inserted_bases": "T",
            "query_start": 5,
            "query_end": 5,
        },))
        self.assertEqual(ins_alignment.query_coordinates[-1], 5)

        del_alignment = next(item for item in deletion.alignments if item.query_rotation == 4)
        self.assertEqual(del_alignment.query_rotation, 4)
        self.assertEqual(del_alignment.cigar, "4=1D3=")
        self.assertEqual(del_alignment.variants, ({
            "type": "deletion",
            "start": 5,
            "end": 5,
            "reference_bases": "T",
            "query_after_position": 1,
        },))

    def test_linear_and_circular_interpretations_are_comparable(self) -> None:
        reference = "ACCGTACG"
        query = rotate(reference, 3)
        self.assertEqual(query, "GTACGACC")
        circular = AlignmentParameters(mode="global", topology="circular")
        linear = AlignmentParameters(mode="global", topology="linear")
        _, _, circular_results = self._run(reference, f">q\n{query}\n", circular)
        _, _, linear_results = self._run(reference, f">q\n{query}\n", linear)
        self.assertEqual(circular_results[0].score, 16)
        self.assertLess(linear_results[0].score, 16)

        circular_local = AlignmentParameters(mode="local", topology="circular")
        linear_local = AlignmentParameters(mode="local", topology="linear")
        _, _, circular_local_results = self._run(reference, f">q\n{query}\n", circular_local)
        _, _, linear_local_results = self._run(reference, f">q\n{query}\n", linear_local)
        self.assertEqual(circular_local_results[0].score, 16)
        self.assertLess(linear_local_results[0].score, 16)
        best = circular_local_results[0].alignments[0]
        self.assertEqual(best.query_coordinates, (6, 7, 8, 1, 2, 3, 4, 5))
        self.assertEqual(best.query_rotation, 5)

    def test_circular_batch_reload_preserves_candidate_order_strand_and_start(self) -> None:
        reference = "ACCGTACG"
        query = rotate(reverse_complement(reference), 3)
        parameters = AlignmentParameters(mode="global", topology="circular")
        run_id, _, results = self._run(reference, f">rev\n{query}\n", parameters)
        store = ProjectStore(self.project)
        reloaded = store.get_run_results(run_id)[0]
        original = results[0]
        self.assertEqual(
            [(item.strand, item.query_rotation, item.query_coordinates, item.score) for item in original.alignments],
            [(item.strand, item.query_rotation, item.query_coordinates, item.score) for item in reloaded.alignments],
        )
        self.assertEqual(reloaded.alignments[0].strand, "reverse")
        self.assertEqual(reloaded.alignments[0].query_rotation, 3)
        manifest = store.get_run_manifest(run_id)
        self.assertEqual(manifest["parameters"]["topology"], "circular")


if __name__ == "__main__":
    unittest.main()
