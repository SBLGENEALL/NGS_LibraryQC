#!/usr/bin/env python3
"""Create a small deterministic paired-end test dataset."""

from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path


LEFT = "CTATAAAAGAGCTCACAACCCCTCA"
RIGHT_PRIMER = "GGAGGCCACACCCGCCACTCACCTG"


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


RIGHT_FORWARD = rc(RIGHT_PRIMER)
ADAPTER_PAD = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" * 5


def fastq_record(name: str, sequence: str, index: str, read_number: int) -> str:
    sequence = sequence[:150]
    quality = "I" * len(sequence)
    return f"@TEST:1:FLOW:1:1:1:{name} {read_number}:N:0:{index}\n{sequence}\n+\n{quality}\n"


def main() -> None:
    random.seed(124)
    base = Path(__file__).resolve().parent / "synthetic"
    fastq = base / "fastq"
    fastq.mkdir(parents=True, exist_ok=True)

    references = []
    used = set()
    for i in range(20):
        while True:
            length = 80 + (i * 7) % 51
            sequence = "".join(random.choice("ACGT") for _ in range(length))
            if sequence not in used:
                used.add(sequence)
                break
        references.append((f"target_{i + 1:04d}", sequence))

    with (base / "design.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Variant_ID", "target_sequence"])
        writer.writerows(references)

    with (base / "design_oligo.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["Variant_ID", "target_sequence"])
        writer.writerows(
            (variant_id, LEFT[-15:] + sequence + RIGHT_FORWARD[:15])
            for variant_id, sequence in references
        )

    assigned_counts = [200, 160, 130, 100, 80, 60, 50, 40, 30, 20] + [10] * 8 + [0, 0]
    undetermined_counts = [20, 16, 13, 10, 8] + [2] * 13 + [0, 0]

    def write_pair(prefix: str, counts, index: str, introduce_error: bool) -> None:
        r1_path = fastq / f"{prefix}_S1_L001_R1_001.fastq.gz"
        r2_path = fastq / f"{prefix}_S1_L001_R2_001.fastq.gz"
        with gzip.open(r1_path, "wt", encoding="utf-8") as out1, gzip.open(
            r2_path, "wt", encoding="utf-8"
        ) as out2:
            number = 0
            for (_, target), count in zip(references, counts):
                for copy in range(count):
                    number += 1
                    spacer_left = "".join(
                        random.choice("ACGT") for _ in range(copy % 6)
                    )
                    spacer_right = "".join(
                        random.choice("ACGT") for _ in range((copy + 2) % 6)
                    )
                    fragment = spacer_left + LEFT + target + RIGHT_FORWARD + spacer_right
                    r1 = fragment + ADAPTER_PAD
                    r2 = rc(fragment) + ADAPTER_PAD
                    if introduce_error and copy == 0 and len(target) > 30:
                        # One recoverable target mismatch in R1 only. R2 consensus is correct.
                        pos = len(spacer_left) + len(LEFT) + 20
                        replacement = "A" if r1[pos] != "A" else "C"
                        r1 = r1[:pos] + replacement + r1[pos + 1 :]
                    out1.write(fastq_record(str(number), r1, index, 1))
                    out2.write(fastq_record(str(number), r2, index, 2))

    write_pair("pDNA", assigned_counts, "ACGTACGT+TGCATGCA", True)
    write_pair(
        "Undetermined",
        undetermined_counts,
        "ACGTACNT+TGCATGCA",
        False,
    )

    phix_path = (
        Path(__file__).resolve().parents[1]
        / "references"
        / "NC_001422.1_phiX174.fasta"
    )
    phix = "".join(
        line.strip()
        for line in phix_path.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )
    circular_phix = phix + phix[:299]
    undetermined_r1 = fastq / "Undetermined_S1_L001_R1_001.fastq.gz"
    undetermined_r2 = fastq / "Undetermined_S1_L001_R2_001.fastq.gz"
    with gzip.open(undetermined_r1, "at", encoding="utf-8") as out1, gzip.open(
        undetermined_r2, "at", encoding="utf-8"
    ) as out2:
        for number in range(50):
            start = (number * 97) % len(phix)
            fragment = circular_phix[start : start + 300]
            out1.write(fastq_record(f"PHIX{number}", fragment, "NNNNNNNN", 1))
            out2.write(fastq_record(f"PHIX{number}", rc(fragment), "NNNNNNNN", 2))

    # Raw index reads permit direct cycle-level Q/base composition checks.
    i1_path = fastq / "Run_S1_L001_I1_001.fastq.gz"
    with gzip.open(i1_path, "wt", encoding="utf-8") as out:
        for i in range(200):
            seq = "ACGTACGT" if i < 180 else "ACGTACNT"
            qual = "I" * 6 + ("!" if "N" in seq else "I") + "I"
            out.write(f"@IDX:{i} 1:N:0:{seq}\n{seq}\n+\n{qual}\n")

    with (base / "SampleSheet.csv").open("w", newline="", encoding="utf-8") as handle:
        handle.write("[Header]\nFileFormatVersion,2\n\n[Data]\n")
        handle.write("Sample_ID,index,index2\n")
        handle.write("pDNA,ACGTACGT,TGCATGCA\n")

    with (base / "Top_Unknown_Barcodes.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "count"])
        writer.writerow(["ACGTACNT+TGCATGCA", sum(undetermined_counts)])
        writer.writerow(["ACGTACAT+TGCATGCA", 7])

    config = f"""[project]
target_name = synthetic_barcode

[input]
input_dir = {fastq}
reference = {base / "design.csv"}
outdir = {base / "results"}
sample_sheet = {base / "SampleSheet.csv"}
top_unknown = {base / "Top_Unknown_Barcodes.csv"}

[reference]
id_column = Variant_ID
sequence_column = target_sequence
reference_mode = auto

[amplicon]
left_primer = {LEFT}
right_primer = {RIGHT_PRIMER}

[matching]
anchor_mismatches = 1
signature_length = 20
max_mismatches = 2

[index_qc]
max_index_distance = 2

[control_qc]
phix_enabled = yes
expected_phix_percent =
phix_kmer_length = 27
phix_min_kmer_hits = 1

[runtime]
max_reads =
write_observed_sequences = yes
overwrite = yes
"""
    (base / "config.ini").write_text(config, encoding="utf-8")

    print(base)


if __name__ == "__main__":
    main()
