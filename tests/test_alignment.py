import tempfile
import unittest
from pathlib import Path

from bioseq.alignment import AlignmentParameters, smith_waterman
from bioseq.runner import run_batch
from bioseq.storage import PlatformStore


PARAMETERS = AlignmentParameters()


class AlignmentTest(unittest.TestCase):
    def test_repeated_region_returns_each_tied_position(self) -> None:
        candidates = smith_waterman("AAT", "AATAAT", PARAMETERS)

        self.assertEqual(len(candidates), 2)
        self.assertEqual([candidate.reference_start for candidate in candidates], [0, 3])
        self.assertEqual([candidate.reference_end for candidate in candidates], [3, 6])
        self.assertTrue(all(candidate.score == 6 for candidate in candidates))
        self.assertTrue(
            all(candidate.aligned_query == "AAT" for candidate in candidates)
        )

    def test_gap_alignment_has_expected_score_and_coordinates(self) -> None:
        candidates = smith_waterman("ACGT", "ACTGT", PARAMETERS)

        self.assertEqual(len(candidates), 1)
        candidate = candidates[0]
        self.assertEqual(candidate.score, 6)
        self.assertEqual((candidate.query_start, candidate.query_end), (0, 4))
        self.assertEqual((candidate.reference_start, candidate.reference_end), (0, 5))
        self.assertEqual(candidate.aligned_query, "AC-GT")
        self.assertEqual(candidate.aligned_reference, "ACTGT")
        self.assertEqual(candidate.cigar, "2=1D2=")

    def test_batch_quarantines_invalid_input_and_is_reproducible(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            store = PlatformStore(Path(directory))
            reference_id = store.put_reference(
                {
                    "name": "repeat-ref",
                    "sequence_type": "dna",
                    "sequence": "AATAAT",
                }
            )
            dataset_id = store.put_dataset(
                {
                    "name": "fixture",
                    "sequence_type": "dna",
                    "sequences": [
                        {"id": "valid", "sequence": "AAT"},
                        {"id": "bad", "sequence": "AAN"},
                    ],
                }
            )

            first = run_batch(store, dataset_id, reference_id, PARAMETERS)
            second = run_batch(store, dataset_id, reference_id, PARAMETERS)

            self.assertEqual(first["batch"]["id"], second["batch"]["id"])
            self.assertEqual(first["results"], second["results"])
            self.assertEqual(first["batch"]["valid_sequence_count"], 1)
            self.assertEqual(first["batch"]["quarantined_sequence_count"], 1)
            self.assertEqual(first["quarantine"][0]["sequence_id"], "bad")
            self.assertIn("invalid characters", first["quarantine"][0]["errors"][0])
            self.assertEqual(
                [item["reference_start"] for item in first["results"][0]["candidates"]],
                [0, 3],
            )


if __name__ == "__main__":
    unittest.main()
