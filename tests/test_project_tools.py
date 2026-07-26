#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from amplicon_qc import read_reference
from project_tools.configure_project import reference_columns, select_optional


class ProjectToolTests(unittest.TestCase):
    def test_reference_column_detection_for_tsv_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "5UTR_reference.tsv"
            path.write_text(
                "Variant_ID\tUTR_sequence\nv1\tACGT\n",
                encoding="utf-8",
            )
            self.assertEqual(reference_columns(path), ("Variant_ID", "UTR_sequence"))

    def test_correct_candidate_sequence_header_is_detected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "5UTR_reference.csv"
            path.write_text(
                "ID,Left Homology,5UTR candidate sequence,Right Homology,"
                "Final_oligo,UTR_Length,Oligo_Length\n"
                "v1,AAAA,ACGT,TTTT,AAAAACGTTTTT,4,12\n",
                encoding="utf-8",
            )
            self.assertEqual(
                reference_columns(path),
                ("ID", "5UTR candidate sequence"),
            )
            library = read_reference(
                path,
                "Variant_ID",
                "UTR_sequence",
                "CTATAAAAGAGCTCACAACCCCTCA",
                "GGAGGCCACACCCGCCACTCACCTG",
                "auto",
            )
            self.assertEqual(library.id_to_sequence, {"v1": "ACGT"})

    def test_optional_candidate_can_be_declined(self):
        candidate = Path("/tmp/Top_Unknown_Barcodes.csv")
        with patch("builtins.input", return_value="no"):
            self.assertIsNone(
                select_optional([candidate], "Top Unknown Barcodes file")
            )


if __name__ == "__main__":
    unittest.main()
