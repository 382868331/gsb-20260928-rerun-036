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
        self.reference.write_text(">ref\nACGTACGT\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _run(self, content: str, parameters: AlignmentParameters):
        self.sequences.write_text(content, encoding="utf-8")
        store = ProjectStore(self.project)
        return run_from_files(
            store,
            reference_path=self.reference,
            sequences_path=self.sequences,
            sequence_type="dna",
            parameters=parameters,
        )

    def test_repeat_returns_all_tied_positions(self) -> None:
        parameters = AlignmentParameters()
        run_id, _, results = self._run(">repeat\nAAA\n", parameters)
        result = results[0]
        self.assertEqual(result.status, "aligned")
        self.assertEqual(result.score, 6)
        self.assertEqual(result.alignment_count, 2)
        self.assertEqual({item.reference_start for item in result.alignments}, {1, 2})
        self.assertEqual({item.reference_end for item in result.alignments}, {3, 4})
        self.assertTrue(run_id.startswith("run-"))

    def test_substitution_score_and_variant(self) -> None:
        parameters = AlignmentParameters(mode="global")
        _, _, results = self._run(">sub\nACGTTCGT\n", parameters)
        alignment = results[0].alignments[0]
        self.assertEqual(results[0].alignment_count, 1)
        self.assertEqual(alignment.score, 13)
        self.assertEqual(alignment.cigar, "4=1X3=")
        self.assertEqual(alignment.variants[0]["type"], "substitution")
        self.assertEqual(alignment.variants[0]["position"], 5)

    def test_deletion_and_insertion_gaps(self) -> None:
        parameters = AlignmentParameters(mode="global")
        _, _, results = self._run(
            ">deletion\nACGACGT\n>insertion\nACGTATACGT\n",
            parameters,
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
        self.assertEqual(insertion.alignments[0].cigar, "4=1I1=4=")
        self.assertEqual(insertion.alignments[0].variants[0]["type"], "insertion")
        self.assertEqual(insertion.alignments[0].variants[0]["after_position"], 4)
        self.assertEqual(insertion.alignments[0].variants[0]["inserted_bases"], "A")
        self.assertEqual(insertion.alignments[0].score, 14)

    def test_invalid_record_is_isolated(self) -> None:
        parameters = AlignmentParameters()
        run_id, manifest, results = self._run(
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


if __name__ == "__main__":
    unittest.main()
