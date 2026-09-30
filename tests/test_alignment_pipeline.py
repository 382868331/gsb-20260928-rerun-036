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

    def _run(
        self,
        content: str,
        parameters: AlignmentParameters,
        reference: str | None = None,
        sequence_type: str = "dna",
    ):
        if reference is not None:
            self.reference.write_text(f">ref\n{reference}\n", encoding="utf-8")
        self.sequences.write_text(content, encoding="utf-8")
        store = ProjectStore(self.project)
        return run_from_files(
            store,
            reference_path=self.reference,
            sequences_path=self.sequences,
            sequence_type=sequence_type,
            parameters=parameters,
        )

    def test_repeat_returns_all_tied_positions(self) -> None:
        parameters = AlignmentParameters()
        run_id, _, results = self._run(">repeat\nAAA\n", parameters)
        result = results[0]
        self.assertEqual(result.status, "aligned")
        self.assertEqual(result.score, 6)
        self.assertEqual(result.alignment_count, 2)
        self.assertEqual({item.orientation for item in result.alignments}, {"forward"})
        self.assertEqual({item.reference_start for item in result.alignments}, {1, 2})
        self.assertEqual({item.reference_end for item in result.alignments}, {3, 4})
        self.assertEqual({item.query_start for item in result.alignments}, {1})
        self.assertEqual({item.query_end for item in result.alignments}, {3})
        self.assertEqual({item.aligned_reference for item in result.alignments}, {"AAA"})
        self.assertEqual(len({item.candidate_id for item in result.alignments}), 2)
        self.assertEqual(result.alignments[0].coordinates, tuple(
            {"reference_position": index, "query_position": index} for index in range(1, 4)
        ))
        self.assertTrue(run_id.startswith("run-"))

    def test_substitution_score_and_variant(self) -> None:
        parameters = AlignmentParameters(mode="global")
        _, _, results = self._run(
            ">sub\nACGTTCGT\n",
            parameters,
            reference="ACGTACGT",
            sequence_type="protein",
        )
        alignment = results[0].alignments[0]
        self.assertEqual(results[0].alignment_count, 1)
        self.assertEqual(alignment.score, 13)
        self.assertEqual(alignment.cigar, "4=1X3=")
        self.assertEqual(alignment.variants[0]["type"], "substitution")
        self.assertEqual(alignment.variants[0]["position"], 5)

    def test_deletion_and_insertion_gaps(self) -> None:
        parameters = AlignmentParameters(mode="global")
        _, _, results = self._run(
            ">deletion\nACGACGT\n>insertion\nACGTAACGT\n",
            parameters,
            reference="ACGTACGT",
        )
        deletion, insertion = results
        self.assertEqual(deletion.alignments[0].cigar, "3=1D4=")
        self.assertEqual(deletion.alignments[0].variants, ({
            "type": "deletion",
            "start": 4,
            "end": 4,
            "reference_bases": "T",
        },))
        self.assertEqual(deletion.alignments[0].score, 12)
        self.assertEqual(insertion.alignments[0].cigar, "4=1I4=")
        self.assertEqual(insertion.alignments[0].variants[0]["type"], "insertion")
        self.assertEqual(insertion.alignments[0].variants[0]["after_position"], 4)
        self.assertEqual(insertion.alignments[0].variants[0]["inserted_bases"], "A")
        self.assertEqual(insertion.alignments[0].score, 14)

    def test_local_match_with_query_prefix_keeps_original_coordinates(self) -> None:
        parameters = AlignmentParameters()
        _, _, results = self._run(
            ">prefixed\nGGAAAA\n",
            parameters,
            reference="TTAAAATT",
        )
        alignment = results[0].alignments[0]
        self.assertEqual(results[0].alignment_count, 1)
        self.assertEqual(alignment.score, 8)
        self.assertEqual((alignment.reference_start, alignment.reference_end), (3, 6))
        self.assertEqual((alignment.query_start, alignment.query_end), (3, 6))
        self.assertEqual(alignment.aligned_reference, "AAAA")
        self.assertEqual(alignment.aligned_query, "AAAA")
        self.assertEqual(alignment.coordinates, tuple(
            {"reference_position": reference_position, "query_position": query_position}
            for reference_position, query_position in zip(range(3, 7), range(3, 7))
        ))

    def test_reverse_match_reports_strand_and_original_query_coordinates(self) -> None:
        parameters = AlignmentParameters(mode="global")
        _, _, results = self._run(
            ">reverse\nACGGTTTT\n",
            parameters,
            reference="AAAACCGT",
        )
        alignment = results[0].alignments[0]
        self.assertEqual(results[0].alignment_count, 1)
        self.assertEqual(alignment.orientation, "reverse")
        self.assertEqual((alignment.reference_strand, alignment.query_strand), ("+", "-"))
        self.assertEqual(alignment.score, 16)
        self.assertEqual((alignment.reference_start, alignment.reference_end), (1, 8))
        self.assertEqual((alignment.query_start, alignment.query_end), (1, 8))
        self.assertEqual(alignment.aligned_reference, "AAAACCGT")
        self.assertEqual(alignment.aligned_query, "AAAACCGT")
        self.assertEqual([position["query_position"] for position in alignment.coordinates], list(range(8, 0, -1)))

    def test_repeated_region_with_indel_keeps_distinct_origins(self) -> None:
        parameters = AlignmentParameters()
        _, _, results = self._run(
            ">repeat_indel\nACGAACG\n",
            parameters,
            reference="ACGACGACG",
        )
        alignments = results[0].alignments
        self.assertEqual(results[0].score, 10)
        self.assertEqual({(item.reference_start, item.reference_end) for item in alignments}, {(1, 6), (4, 9), (1, 9)})
        self.assertEqual({item.cigar for item in alignments}, {"4=1I2=", "3=1I3=", "4=2D3="})
        self.assertEqual(len({item.candidate_id for item in alignments}), 5)

    def test_saved_reverse_batch_can_be_reopened_to_original_coordinates(self) -> None:
        parameters = AlignmentParameters(mode="global")
        first_id, _, first_results = self._run(
            ">reverse\nACGGTTTT\n",
            parameters,
            reference="AAAACCGT",
        )
        reopened_store = ProjectStore(self.project)
        _, _, reopened_results = run_from_files(
            reopened_store,
            reference_path=self.reference,
            sequences_path=self.sequences,
            sequence_type="dna",
            parameters=parameters,
        )
        stored = reopened_store.get_run_results(first_id)[0].alignments[0]
        expected = first_results[0].alignments[0]
        self.assertEqual(reopened_results[0].alignments[0].candidate_id, expected.candidate_id)
        self.assertEqual(stored.orientation, "reverse")
        self.assertEqual((stored.reference_start, stored.reference_end), (1, 8))
        self.assertEqual((stored.query_start, stored.query_end), (1, 8))
        self.assertEqual([position["query_position"] for position in stored.coordinates], list(range(8, 0, -1)))

    def test_invalid_record_is_isolated(self) -> None:
        parameters = AlignmentParameters()
        run_id, manifest, results = self._run(
            ">valid\nACG\n>bad\nACGTX\n",
            parameters,
            reference="ACGTACGT",
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
        self.assertEqual({item.candidate_id for item in second_results[0].alignments},
                         {item.candidate_id for item in first_results[0].alignments})
        self.assertTrue(all(item.coordinates for item in second_results[0].alignments))


if __name__ == "__main__":
    unittest.main()
