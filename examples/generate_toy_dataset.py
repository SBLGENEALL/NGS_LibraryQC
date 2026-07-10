#!/usr/bin/env python3
"""Generate a small paired-end-like toy FASTQ dataset for NGS_LibraryQC.

The toy dataset includes:
- 10 expected reference inserts
- one dropout reference (UTR_010)
- deliberately biased counts
- staggered 1-10 nt random bases before the left anchor
- both forward and reverse-complement reads
- a few non-reference inserts
- a few reads without detectable anchors
"""

from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path

from Bio.Seq import Seq

LEFT = "CTATAAAAGAGCTCACAACCCCTCA"
RIGHT = "GGAGGCCACACCCGCCACTCACCTG"

REFERENCES = {
    "UTR_001": "ACGTTGCAACGTTGCAACGTTGCA",
    "UTR_002": "TGCATGCATGCATGCATGCATGCA",
    "UTR_003": "GGGAAACCCGGGAAACCCGGGAAA",
    "UTR_004": "ATATCGCGATATCGCGATATCGCG",
    "UTR_005": "CAGTCAGTCAGTCAGTCAGTCAGT",
    "UTR_006": "GCTAGCTAGCTAGCTAGCTAGCTA",
    "UTR_007": "AACCGGTTAACCGGTTAACCGGTT",
    "UTR_008": "TTGGAACCTTGGAACCTTGGAACC",
    "UTR_009": "AGCTTGCAAGCTTGCAAGCTTGCA",
    "UTR_010": "CGATCGATCGATCGATCGATCGAT",
}

# UTR_010 is intentionally absent.
COUNTS = {
    "UTR_001": 120,
    "UTR_002": 80,
    "UTR_003": 55,
    "UTR_004": 35,
    "UTR_005": 25,
    "UTR_006": 20,
    "UTR_007": 15,
    "UTR_008": 10,
    "UTR_009": 5,
}

NON_REFERENCE = {
    "NONREF_A": ("ACGTTGCAACGTTGCAACGTTGCT", 7),
    "NONREF_B": ("CCCCAAAATTTTGGGGCCCCAAAA", 4),
}


def random_stagger() -> str:
    length = random.randint(1, 10)
    return "".join(random.choice("ACGT") for _ in range(length))


def write_fastq_record(handle, name: str, seq: str) -> None:
    quality = "I" * len(seq)
    handle.write(f"@{name}\n{seq}\n+\n{quality}\n")


def main() -> None:
    random.seed(20260710)
    root = Path(__file__).resolve().parents[1]
    out_dir = root / "examples" / "toy_data"
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_path = out_dir / "toy_reference.csv"
    fastq_path = out_dir / "toy_reads.fastq.gz"

    with ref_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["utr_id", "utr_seq"])
        for utr_id, seq in REFERENCES.items():
            writer.writerow([utr_id, seq])

    records = []
    read_index = 1

    for utr_id, count in COUNTS.items():
        insert = REFERENCES[utr_id]
        for _ in range(count):
            seq = random_stagger() + LEFT + insert + RIGHT + random_stagger()
            if read_index % 2 == 0:
                seq = str(Seq(seq).reverse_complement())
            records.append((f"toy_{read_index}_{utr_id}", seq))
            read_index += 1

    for label, (insert, count) in NON_REFERENCE.items():
        for _ in range(count):
            seq = random_stagger() + LEFT + insert + RIGHT + random_stagger()
            if read_index % 2 == 0:
                seq = str(Seq(seq).reverse_complement())
            records.append((f"toy_{read_index}_{label}", seq))
            read_index += 1

    for _ in range(6):
        seq = "".join(random.choice("ACGT") for _ in range(100))
        records.append((f"toy_{read_index}_NOANCHOR", seq))
        read_index += 1

    random.shuffle(records)
    with gzip.open(fastq_path, "wt") as handle:
        for name, seq in records:
            write_fastq_record(handle, name, seq)

    print(f"Created: {ref_path}")
    print(f"Created: {fastq_path}")
    print(f"Total reads: {len(records)}")
    print("Expected detected references: 9 / 10")
    print("Expected dropout: UTR_010")
    print("Expected non-reference reads: 11")
    print("Expected anchor-missing reads: 6")


if __name__ == "__main__":
    main()
