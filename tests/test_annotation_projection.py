from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioseq_platform.annotations import project_annotation, region_members
from bioseq_platform.models import Alignment, AlignmentParameters, FunctionalRegion
from bioseq_platform.pipeline import (
    annotation_from_file,
    build_projections,
    projection_preview,
    run_from_files,
    select_alignment,
)
from bioseq_platform.storage import ProjectStore


def rotate(sequence: str, offset: int) -> str:
    return sequence[offset:] + sequence[:offset]


class RegionModelTest(unittest.TestCase):
    def test_wrap_around_region_is_one_directed_interval(self) -> None:
        region = FunctionalRegion(id="R1", start=7, end=2)
        self.assertEqual(region_members(region, 8), [7, 8, 1, 2])
        self.assertTrue(region.end < region.start)

    def test_reverse_region_reverses_member_order(self) -> None:
        region = FunctionalRegion(id="R2", start=7, end=2, strand="reverse")
        self.assertEqual(region_members(region, 8), [2, 1, 8, 7])

    def test_region_validates_against_reference_length(self) -> None:
        with self.assertRaises(ValueError):
            FunctionalRegion(id="R3", start=1, end=9).validate(8)
        with self.assertRaises(ValueError):
            FunctionalRegion(id="R4", start=0, end=2).validate(8)


class AnnotationProjectionTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.project = Path(self.temporary.name) / "project"
        self.reference = self.project / "reference.fasta"
        self.sequences = self.project / "sequences.fasta"
        self.annotation = self.project / "regions.tsv"
        self.reference.parent.mkdir()
        self.store = ProjectStore(self.project)
        self.store.initialize()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _run(self, reference: str, content: str, parameters: AlignmentParameters | None = None):
        parameters = parameters or AlignmentParameters(mode="global", topology="circular")
        self.reference.write_text(f">ref\n{reference}\n", encoding="utf-8")
        self.sequences.write_text(content, encoding="utf-8")
        return run_from_files(
            self.store,
            reference_path=self.reference,
            sequences_path=self.sequences,
            sequence_type="dna",
            parameters=parameters,
        )

    def _annotate(self, rows: list[str]) -> str:
        self.annotation.write_text(
            "\n".join(rows) + "\n",
            encoding="utf-8",
        )
        return annotation_from_file(
            self.store,
            annotation_path=self.annotation,
            reference_id="ref",
        )

    def _alignment(self, run_id: str, record_id: str, index: int = 0) -> Alignment:
        return self.store.get_run_results(run_id)[0].alignments[index]

    def test_crossing_region_projects_as_single_query_range(self) -> None:
        run_id, _, _ = self._run("ACGT", ">q\nGTAC\n")
        aid = self._annotate(["wrap	3	2	+	wrap-feature"])
        payload = build_projections(
            self.store, run_id=run_id, record_id="q",
            annotation_id=aid,
        )
        region = payload["candidates"][0]["regions"][0]
        self.assertEqual(region["region_id"], "wrap")
        self.assertTrue(region["region_wraps"])
        self.assertEqual(region["status"], "covered")
        self.assertEqual(region["query_range"]["start"], 1)
        self.assertEqual(region["query_range"]["end"], 4)
        self.assertEqual(
            region["projected_reference_segments"],
            [{"start": 3, "end": 2, "wraps": True}],
        )

    def test_boundary_insertion_counted_once_for_upstream_region(self) -> None:
        run_id, _, _ = self._run("ACGTACGT", ">ins\nACGTAACGT\n")
        aid = self._annotate(["A	1	4	+", "B	5	8	+"])
        payload = build_projections(
            self.store, run_id=run_id, record_id="ins", annotation_id=aid
        )
        anchored = [
            candidate for candidate in payload["candidates"]
            if candidate["strand"] == "forward"
            and candidate["regions"][0]["variants"]["insertions"]
            and candidate["regions"][0]["variants"]["insertions"][0]["after_reference_position"] == 4
        ]
        self.assertEqual(len(anchored), 1)
        regions = {item["region_id"]: item for item in anchored[0]["regions"]}
        self.assertEqual(
            regions["A"]["variants"]["insertions"][0]["location"], "end_boundary"
        )
        self.assertEqual(
            regions["A"]["variants"]["insertions"][0]["region_offset"], 4
        )
        self.assertEqual(regions["B"]["variants"]["insertions"], [])
        total_insertions = sum(
            len(item["variants"]["insertions"])
            for candidate in payload["candidates"]
            for item in candidate["regions"]
        )
        self.assertGreaterEqual(total_insertions, 1)

    def test_deletion_covering_whole_region_has_no_query_range(self) -> None:
        run_id, _, _ = self._run("ACGTACGT", ">del\nACGTACGT\n")
        # Query equals reference; instead force a full deletion explicitly.
        alignment = Alignment(
            score=8,
            reference_start=3,
            reference_end=6,
            query_start=1,
            query_end=4,
            aligned_reference="GTAC",
            aligned_query="G--C",
            cigar="1=2D1=",
            strand="forward",
            rotation_offset=0,
            wraps_reference=False,
            variants=({
                "type": "deletion",
                "start": 4,
                "end": 5,
                "reference_bases": "TA",
            },),
            reference_coordinates=(3, 4, 5, 6),
            query_coordinates=(1, 0, 0, 2),
        )
        region = FunctionalRegion(id="middle", start=4, end=5)
        projection = project_annotation([region], alignment, 8)[0]
        self.assertEqual(projection["status"], "deleted")
        self.assertIsNone(projection["query_range"])
        self.assertEqual(
            projection["deleted_reference_segments"],
            [{"start": 4, "end": 5, "wraps": False}],
        )
        self.assertTrue(projection["variants"]["deletions"][0]["removes_region"])

    def test_reverse_complement_returns_to_original_query_coordinates(self) -> None:
        reference = "AACGTTGC"
        from bioseq_platform.alignment import reverse_complement
        query = reverse_complement(rotate(reference, 3))
        run_id, _, _ = self._run(reference, f">rev\n{query}\n")
        aid = self._annotate(["fwd	1	4	+", "rev	5	8	-"])
        payload = build_projections(
            self.store, run_id=run_id, record_id="rev", annotation_id=aid
        )
        reverse_candidates = [c for c in payload["candidates"] if c["strand"] == "reverse"]
        self.assertTrue(reverse_candidates)
        candidate = reverse_candidates[0]
        for region in candidate["regions"]:
            if region["status"] == "covered":
                self.assertEqual(region["query_range"]["strand"], "reverse")
                self.assertGreaterEqual(region["query_range"]["start"], region["query_range"]["end"])

    def test_repeat_tie_keeps_distinct_projections(self) -> None:
        run_id, _, _ = self._run("AATAAT", ">rep\nAATAAT\n")
        aid = self._annotate(["first	1	3	+", "second	4	6	+"])
        payload = build_projections(
            self.store, run_id=run_id, record_id="rep", annotation_id=aid
        )
        forward = [c for c in payload["candidates"] if c["strand"] == "forward"]
        offsets = {c["rotation_offset"] for c in forward}
        self.assertEqual(offsets, {0, 3})
        region_sets = [
            tuple((r["region_id"], r["query_range"]) for r in c["regions"])
            for c in forward
        ]
        self.assertEqual(len({repr(item) for item in region_sets}), 2)

    def test_reselect_previews_regions_that_need_recompute(self) -> None:
        run_id, _, _ = self._run("AATAAT", ">rep\nAATAAT\n")
        aid = self._annotate(["first	1	3	+", "second	4	6	+"])
        build_projections(
            self.store, run_id=run_id, record_id="rep", annotation_id=aid
        )
        select_alignment(self.store, run_id, "rep", 0)
        build_projections(
            self.store, run_id=run_id, record_id="rep", annotation_id=aid
        )
        preview = projection_preview(
            self.store,
            run_id=run_id,
            record_id="rep",
            annotation_id=aid,
            target_alignment_index=1,
        )
        self.assertTrue(preview["recompute_required"])
        self.assertGreater(preview["changed_region_count"], 0)
        self.assertTrue(all(item["needs_recompute"] for item in preview["changed_regions"]))

    def test_projection_and_selection_survive_reopen(self) -> None:
        run_id, _, _ = self._run("AATAAT", ">rep\nAATAAT\n")
        aid = self._annotate(["first	1	3	+"])
        build_projections(
            self.store, run_id=run_id, record_id="rep", annotation_id=aid
        )
        select_alignment(self.store, run_id, "rep", 1)
        reopened_store = ProjectStore(self.project)
        reloaded = reopened_store.get_run_results(run_id)
        self.assertEqual(reloaded[0].selected_alignment_index, 1)
        projection = reopened_store.load_projection(run_id, "rep", aid)
        self.assertIsNotNone(projection)
        self.assertEqual(len(projection["candidates"]), 2)

    def test_rotation_equivalence_both_molecules(self) -> None:
        reference = "AACGTTGC"
        from bioseq_platform.alignment import reverse_complement
        query = reverse_complement(rotate(reference, 3))
        run_id_a, _, _ = self._run(
            reference,
            f">q\n{query}\n",
        )
        # Rotate both reference and query by the same amount; the biological
        # event is identical, only the record start moves.
        shift = 2
        run_id_b, _, _ = self._run(
            rotate(reference, shift),
            f">q\n{rotate(query, shift)}\n",
        )

        def rotations_with_changes(run_id: str):
            results = self.store.get_run_results(run_id)
            summary = []
            for candidate in results[0].alignments:
                if candidate.strand != "reverse":
                    continue
                summary.append((
                    candidate.rotation_offset,
                    candidate.cigar,
                    tuple(
                        (r, q)
                        for r, q in zip(candidate.reference_coordinates, candidate.query_coordinates)
                        if r and q
                    ),
                ))
            return sorted(summary)

        self.assertTrue(rotations_with_changes(run_id_a))
        self.assertEqual(
            [item[1] for item in rotations_with_changes(run_id_a)],
            [item[1] for item in rotations_with_changes(run_id_b)],
        )


def _rotate_regions(regions: list[FunctionalRegion], shift: int, reference_length: int):
    shifted = []
    for region in regions:
        members = region_members(region, reference_length)
        rotated = [((position - 1 + shift) % reference_length) + 1 for position in members]
        shifted.append(FunctionalRegion(
            id=region.id,
            start=rotated[0],
            end=rotated[-1],
            strand=region.strand,
            name=region.name,
        ))
    return shifted


class RotationEquivalenceTest(unittest.TestCase):
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

    def test_functional_impact_equivalent_when_both_molecules_rotate(self) -> None:
        from bioseq_platform.alignment import reverse_complement

        reference = "AACGTTGC"
        shift = 3
        rotated_query = rotate(reverse_complement(reference), shift)
        regions = [FunctionalRegion(id="f", start=2, end=5), FunctionalRegion(id="r", start=6, end=1, strand="reverse")]

        def project(reference_sequence: str, query_sequence: str, annotation_regions):
            self.reference.write_text(f">ref\n{reference_sequence}\n", encoding="utf-8")
            self.sequences.write_text(f">q\n{query_sequence}\n", encoding="utf-8")
            run_id, _, _ = run_from_files(
                self.store,
                reference_path=self.reference,
                sequences_path=self.sequences,
                sequence_type="dna",
                parameters=AlignmentParameters(mode="global", topology="circular"),
            )
            aid = self.store.add_annotation(annotation_regions, reference_id="ref")
            payload = build_projections(
                self.store, run_id=run_id, record_id="q", annotation_id=aid
            )
            selected = next(
                candidate for candidate in payload["candidates"]
                if candidate["strand"] == "reverse"
            )
            return {region["region_id"]: region for region in selected["regions"]}

        before = project(reference, rotated_query, regions)
        after = project(
            rotate(reference, shift),
            rotate(rotated_query, shift),
            _rotate_regions(regions, shift, len(reference)),
        )
        for region_id in ("f", "r"):
            self.assertEqual(before[region_id]["status"], after[region_id]["status"])
            before_range = before[region_id]["query_range"]
            after_range = after[region_id]["query_range"]
            self.assertEqual(before_range["strand"], after_range["strand"])
            self.assertEqual(
                len(before_range["segments"]),
                len(after_range["segments"]),
            )
            self.assertEqual(
                before[region_id]["variants"],
                after[region_id]["variants"],
            )


if __name__ == "__main__":
    unittest.main()
