from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioseq_platform.models import AlignmentParameters, Scoring
from bioseq_platform.pipeline import run_from_files
from bioseq_platform.storage import ProjectStore
from bioseq_platform.alignment import reverse_complement


class AlignmentPipelineTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.project = Path(self.temporary.name) / "project"
        self.reference = self.project / "reference.fasta"
        self.sequences = self.project / "sequences.fasta"
        self.reference.parent.mkdir()
        self.reference.write_text(">ref\nAAAA\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _run(self, content: str, parameters: AlignmentParameters):
        self.sequences.write_text(content, encoding="utf-8")
        store = ProjectStore(self.project)
        results = run_from_files(
            store,
            reference_path=self.reference,
            sequences_path=self.sequences,
            sequence_type="dna",
            parameters=parameters,
        )
        return results

    def _run_against(self, reference: str, content: str, parameters: AlignmentParameters):
        self.reference.write_text(f">ref\n{reference}\n", encoding="utf-8")
        return self._run(content, parameters)

    def test_repeat_returns_all_tied_positions(self) -> None:
        parameters = AlignmentParameters()
        run_id, _, results = self._run(">repeat\nAAA\n", parameters)
        result = results[0]
        self.assertEqual(result.status, "aligned")
        self.assertEqual(result.score, 6)
        self.assertEqual(result.alignment_count, 2)
        self.assertEqual({item.reference_start for item in result.alignments}, {1, 2})
        self.assertEqual({item.reference_end for item in result.alignments}, {3, 4})
        starts = sorted(item.reference_start for item in result.alignments)
        self.assertTrue(all(item.strand == "forward" for item in result.alignments))
        self.assertTrue(all(item.reference_coordinates == tuple(range(start, start + 3)) for item, start in zip(result.alignments, starts)))
        self.assertTrue(all(item.aligned_reference == "AAA" for item in result.alignments))
        self.assertTrue(run_id.startswith("run-"))

    def test_repeat_with_short_indel_keeps_coordinate_identities(self) -> None:
        parameters = AlignmentParameters()
        _, _, results = self._run_against("AAAAAAAAAA", ">run\nAAAAAAAA\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 16)
        self.assertEqual(result.alignment_count, 3)
        expected = [(1, 8), (2, 9), (3, 10)]
        actual = [(item.reference_start, item.reference_end) for item in result.alignments]
        self.assertEqual(actual, expected)
        for alignment, (start, end) in zip(result.alignments, expected):
            self.assertEqual(alignment.aligned_reference, "AAAAAAAA")
            self.assertEqual(alignment.reference_coordinates, tuple(range(start, end + 1)))
            self.assertEqual(alignment.query_coordinates, (1, 2, 3, 4, 5, 6, 7, 8))

    def test_local_match_with_prefix_locates_original_interval(self) -> None:
        parameters = AlignmentParameters()
        _, _, results = self._run_against("TTTACGT", ">probe\nACGT\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 8)
        forward = [item for item in result.alignments if item.strand == "forward"]
        self.assertEqual(len(forward), 1)
        alignment = forward[0]
        self.assertEqual((alignment.reference_start, alignment.reference_end), (4, 7))
        self.assertEqual((alignment.query_start, alignment.query_end), (1, 4))
        self.assertEqual(alignment.reference_coordinates, (4, 5, 6, 7))
        self.assertEqual(alignment.cigar, "4=")

    def test_reverse_match_reports_strand_and_original_query_coordinates(self) -> None:
        parameters = AlignmentParameters()
        reference = "AAGCT"
        query = reverse_complement(reference)
        self.assertEqual(query, "AGCTT")
        _, _, results = self._run_against(reference, f">rev\n{query}\n", parameters)
        result = results[0]
        self.assertEqual(result.score, 10)
        reverse = [item for item in result.alignments if item.strand == "reverse"]
        self.assertEqual(len(reverse), 1)
        alignment = reverse[0]
        self.assertEqual((alignment.reference_start, alignment.reference_end), (1, 5))
        self.assertEqual((alignment.query_start, alignment.query_end), (5, 1))
        self.assertEqual(alignment.aligned_reference, alignment.aligned_query)
        self.assertEqual(alignment.aligned_reference, "AAGCT")
        self.assertEqual(alignment.query_coordinates, (5, 4, 3, 2, 1))
        self.assertEqual(alignment.trace_end_query, 1)
        self.assertEqual(alignment.trace_end_reference, 5)

    def test_substitution_score_and_variant(self) -> None:
        parameters = AlignmentParameters(mode="global")
        _, _, results = self._run_against("ACGTACGT", ">sub\nACGTTCGT\n", parameters)
        alignment = results[0].alignments[0]
        self.assertEqual(results[0].alignment_count, 1)
        self.assertEqual(alignment.score, 13)
        self.assertEqual(alignment.cigar, "4=1X3=")
        self.assertEqual(alignment.variants[0]["type"], "substitution")
        self.assertEqual(alignment.variants[0]["position"], 5)

    def test_deletion_and_insertion_gaps(self) -> None:
        parameters = AlignmentParameters(mode="global")
        _, _, results = self._run_against(
            "ACGTACGT",
            ">deletion\nACGACGT\n>insertion\nACGTATACGT\n",
            parameters,
        )
        deletion, insertion = results
        self.assertEqual(deletion.alignment_count, 1)
        self.assertEqual(deletion.alignments[0].cigar, "3=1D4=")
        self.assertEqual(deletion.alignments[0].variants, ({
            "type": "deletion",
            "start": 4,
            "end": 4,
            "reference_bases": "T",
        },))
        self.assertEqual(deletion.alignments[0].score, 12)
        self.assertEqual(insertion.alignment_count, 3)
        self.assertEqual({item.cigar for item in insertion.alignments}, {"5=2I3=", "4=2I4=", "3=2I5="})
        self.assertTrue(all(item.score == 12 for item in insertion.alignments))
        self.assertEqual([item.variants[0]["after_position"] for item in insertion.alignments], [5, 4, 3])
        self.assertTrue(all(item.variants[0]["type"] == "insertion" for item in insertion.alignments))
        inserted = {item.variants[0]["after_position"]: item.variants[0]["inserted_bases"] for item in insertion.alignments}
        self.assertEqual(inserted, {5: "TA", 4: "AT", 3: "TA"})

    def test_invalid_record_is_isolated(self) -> None:
        parameters = AlignmentParameters()
        run_id, manifest, results = self._run_against(
            "ACGTACGT",
            ">valid\nACG\n>bad\nACGTX\n",
            parameters,
        )
        self.assertEqual(len(results), 2)
        valid, invalid = results
        self.assertEqual(valid.status, "aligned")
        self.assertEqual(valid.score, 6)
        self.assertTrue(invalid.status.startswith("invalid:"))
        self.assertEqual(invalid.alignment_count, 0)
        self.assertEqual(manifest["valid_record_count"], 1)
        self.assertEqual(manifest["invalid_record_count"], 1)
        self.assertTrue((Path(self.project) / "runs" / run_id / "results.json").exists())

    def test_reopened_project_reproduces_same_batch(self) -> None:
        parameters = AlignmentParameters()
        first_id, _, first_results = self._run(">repeat\nAAA\n", parameters)
        reopened_store = ProjectStore(self.project)
        second_id, second_manifest, second_results = run_from_files(
            reopened_store,
            reference_path=self.reference,
            sequences_path=self.sequences,
            sequence_type="dna",
            parameters=parameters,
        )
        self.assertEqual(second_id, first_id)
        self.assertEqual(second_manifest["online_databases"], False)
        self.assertEqual(len(second_results), len(first_results))
        self.assertEqual({item.reference_start for item in second_results[0].alignments}, {1, 2})

    def test_saved_batch_reloads_original_coordinates_and_strand(self) -> None:
        parameters = AlignmentParameters()
        run_id, _, first_results = self._run(">repeat\nAAA\n", parameters)
        store = ProjectStore(self.project)
        reloaded_results = store.get_run_results(run_id)
        first = first_results[0].alignments
        reloaded = reloaded_results[0].alignments
        self.assertEqual(len(first), len(reloaded))
        self.assertEqual(
            [(item.reference_start, item.reference_end, item.query_start, item.query_end, item.strand, item.reference_coordinates, item.query_coordinates) for item in first],
            [(item.reference_start, item.reference_end, item.query_start, item.query_end, item.strand, item.reference_coordinates, item.query_coordinates) for item in reloaded],
        )
        manifest = store.get_run_manifest(run_id)
        self.assertEqual(manifest["reference"]["sequence"], "AAAA")


if __name__ == "__main__":
    unittest.main()
