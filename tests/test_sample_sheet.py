#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path

from amplicon_qc import parse_sample_sheet


class SampleSheetTests(unittest.TestCase):
    def parse(self, content: str):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "SampleSheet.csv"
            path.write_text(content, encoding="utf-8")
            return parse_sample_sheet(path)

    def test_illumina_sectioned_sample_sheet(self):
        rows = self.parse(
            "[Header]\n"
            "FileFormatVersion,2\n\n"
            "[Reads]\n151\n151\n\n"
            "[Settings]\nAdapter,AGATCGGAAGAGC\n\n"
            "[Data]\n"
            "Sample_ID,index,index2\n"
            "pDNA,ACGTACGT,TGCATGCA\n"
        )
        self.assertEqual(
            rows,
            [{"sample": "PDNA", "i7": "ACGTACGT", "i5": "TGCATGCA"}],
        )

    def test_surplus_trailing_columns_do_not_raise(self):
        rows = self.parse(
            "[Header]\n"
            "FileFormatVersion,2\n\n"
            "[Data]\n"
            "Sample_ID,index,index2\n"
            "pDNA,ACGTACGT,TGCATGCA,,\n"
        )
        self.assertEqual(
            rows,
            [{"sample": "PDNA", "i7": "ACGTACGT", "i5": "TGCATGCA"}],
        )

    def test_flat_csv_remains_supported(self):
        rows = self.parse(
            "Sample_ID,index,index2\n"
            "pDNA,ACGTACGT,TGCATGCA\n"
        )
        self.assertEqual(
            rows,
            [{"sample": "PDNA", "i7": "ACGTACGT", "i5": "TGCATGCA"}],
        )


if __name__ == "__main__":
    unittest.main()
