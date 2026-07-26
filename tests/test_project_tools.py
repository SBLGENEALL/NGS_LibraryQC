#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path

from project_tools.cleanup_project import (
    discover_fastq_groups,
    discover_metadata,
    discover_references,
    discover_result_dirs,
    guess_fastq_role,
)
from project_tools.configure_project import reference_columns


class ProjectToolTests(unittest.TestCase):
    def test_inventory_excludes_pipeline_test_data(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "03_NGS"
            package = root / "NGS_LibraryQC-1.3"
            package.mkdir(parents=True)

            raw = root / "01_5UTR_Plasmid" / "02_RawData" / "2026-07-23"
            raw.mkdir(parents=True)
            (raw / "pDNA_R1_001.fastq.gz").write_bytes(b"raw")

            subset = root / "01_5UTR_Plasmid" / "fastq_1pct"
            subset.mkdir(parents=True)
            (subset / "pDNA_1pct_R1.fastq.gz").write_bytes(b"subset")

            reference = root / "01_5UTR_Plasmid" / "01_Reference"
            reference.mkdir(parents=True)
            (reference / "5UTR_reference.tsv").write_text(
                "Variant_ID\tUTR_sequence\nv1\tACGT\n",
                encoding="utf-8",
            )

            metadata = raw / "SampleSheet.csv"
            metadata.write_text(
                "[Data]\nSample_ID,index,index2\npDNA,ACGT,TGCA\n",
                encoding="utf-8",
            )

            result = root / "01_5UTR_Plasmid" / "03_Analysis" / "v1.2" / "results"
            combined = result / "combined"
            combined.mkdir(parents=True)
            (combined / "ALL_ASSIGNED_variant_counts.tsv").write_text(
                "variant_ids\ttotal_count\nv1\t1\n",
                encoding="utf-8",
            )

            synthetic = (
                root
                / "00_Tools"
                / "NGS_LibraryQC-v1.2"
                / "tests"
                / "synthetic"
                / "fastq"
            )
            synthetic.mkdir(parents=True)
            (synthetic / "toy_R1.fastq.gz").write_bytes(b"toy")

            results = discover_result_dirs(root, package)
            groups = discover_fastq_groups(root, package, results)
            references = discover_references(root, package, results)
            metadata_files = discover_metadata(root, package, results)

            self.assertEqual(results, [result])
            self.assertEqual({parent for parent, _ in groups}, {raw, subset})
            self.assertEqual(guess_fastq_role(raw), "raw")
            self.assertEqual(guess_fastq_role(subset), "subset")
            self.assertEqual(references, [reference / "5UTR_reference.tsv"])
            self.assertEqual(metadata_files, [metadata])

    def test_reference_column_detection(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "5UTR_reference.tsv"
            path.write_text(
                "Variant_ID\tUTR_sequence\nv1\tACGT\n",
                encoding="utf-8",
            )
            self.assertEqual(reference_columns(path), ("Variant_ID", "UTR_sequence"))


if __name__ == "__main__":
    unittest.main()
