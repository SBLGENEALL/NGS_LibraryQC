#!/usr/bin/env python3
"""Create the project INI by discovery and yes/no selection."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


LEFT_PRIMER = "CTATAAAAGAGCTCACAACCCCTCA"
RIGHT_PRIMER = "GGAGGCCACACCCGCCACTCACCTG"


def ask(question: str) -> bool:
    while True:
        answer = input(f"{question} [yes/no]: ").strip().lower()
        if answer in {"yes", "y"}:
            return True
        if answer in {"no", "n"}:
            return False
        print("Please answer yes or no.")


def select_one(paths, label):
    if not paths:
        raise RuntimeError(f"No {label} candidate was found.")
    if len(paths) == 1:
        path = paths[0]
        if ask(f"Use this {label}?\n  {path}"):
            return path
        raise RuntimeError(f"No {label} was selected.")
    for path in paths:
        if ask(f"Use this {label}?\n  {path}"):
            return path
    raise RuntimeError(f"No {label} was selected.")


def select_optional(paths, label):
    for path in paths:
        if ask(f"Use this optional {label}?\n  {path}"):
            return path
    print(f"No {label} selected; continuing without it.")
    return None


def relative_to_config(path: Path, config_dir: Path) -> str:
    return str(Path("..") / path.relative_to(config_dir.parent))


def reference_columns(path: Path):
    if path.suffix.lower() in {".fa", ".fasta", ".fna"}:
        return "Variant_ID", "target_sequence"
    first = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()[0]
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    header = next(csv.reader([first], delimiter=delimiter))
    lookup = {re.sub(r"[^a-z0-9]", "", item.lower()): item.strip() for item in header}
    id_column = next(
        (
            lookup[key]
            for key in ("variantid", "sequenceid", "designid", "id", "name")
            if key in lookup
        ),
        "Variant_ID",
    )
    sequence_column = next(
        (
            lookup[key]
            for key in (
                "utrsequence",
                "5utrsequence",
                "5utrcandidatesequence",
                "targetsequence",
                "sequence",
                "seq",
                "finaloligo",
            )
            if key in lookup
        ),
        "UTR_sequence",
    )
    return id_column, sequence_column


def main() -> int:
    pipeline = Path(__file__).resolve().parents[1]
    project = pipeline.parent
    config_dir = project / "config"
    config_path = config_dir / "libraryqc.ini"
    reference_dir = project / "reference"

    references = sorted(
        path
        for path in reference_dir.iterdir()
        if path.is_file()
        and path.suffix.lower() in {".csv", ".tsv", ".txt", ".fa", ".fasta", ".fna"}
    )
    reference = select_one(references, "reference")
    sample_sheets = sorted(
        path
        for path in config_dir.iterdir()
        if path.is_file()
        and "samplesheet" in re.sub(r"[^a-z0-9]", "", path.name.lower())
    )
    sample_sheet = select_optional(sample_sheets, "SampleSheet") if sample_sheets else None
    top_unknowns = sorted(
        path
        for path in config_dir.iterdir()
        if path.is_file()
        and "topunknownbarcodes" in re.sub(r"[^a-z0-9]", "", path.name.lower())
    )
    top_unknown = (
        select_optional(top_unknowns, "Top Unknown Barcodes file")
        if top_unknowns
        else None
    )
    id_column, sequence_column = reference_columns(reference)

    if config_path.exists() and not ask(f"Replace the existing configuration?\n  {config_path}"):
        print("Existing configuration kept.")
        return 0

    sample_value = relative_to_config(sample_sheet, config_dir) if sample_sheet else ""
    unknown_value = relative_to_config(top_unknown, config_dir) if top_unknown else ""
    content = f"""[project]
target_name = 5UTR

[input]
# run.sh overrides input_dir and outdir for 1pct/full runs.
input_dir = ../raw_data/1pct
reference = {relative_to_config(reference, config_dir)}
outdir = ../results/1pct
sample_sheet = {sample_value}
top_unknown = {unknown_value}

[reference]
id_column = {id_column}
sequence_column = {sequence_column}
reference_mode = auto

[amplicon]
left_primer = {LEFT_PRIMER}
right_primer = {RIGHT_PRIMER}

[matching]
anchor_mismatches = 1
signature_length = 20
max_mismatches = 2

[index_qc]
max_index_distance = 2

[control_qc]
phix_enabled = yes
phix_reference =
expected_phix_percent = 10
phix_kmer_length = 27
phix_min_kmer_hits = 1

[runtime]
max_reads =
write_observed_sequences = yes
overwrite = no
"""
    config_path.write_text(content, encoding="utf-8")
    print(f"Configuration created: {config_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, IndexError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
