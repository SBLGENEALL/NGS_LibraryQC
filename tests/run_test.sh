#!/usr/bin/env bash
set -euo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd -- "${TEST_DIR}/.." && pwd)"

python3 "${TEST_DIR}/make_synthetic_data.py"

python3 "${PIPELINE_DIR}/amplicon_qc.py" \
  --config "${TEST_DIR}/synthetic/config.ini"

python3 - "${TEST_DIR}/synthetic/results" <<'PY'
import csv
import sys
from pathlib import Path

result_dir = Path(sys.argv[1])
pipeline_dir = result_dir.parents[2]
sys.path.insert(0, str(pipeline_dir))
from amplicon_qc import read_reference

LEFT_PRIMER = "CTATAAAAGAGCTCACAACCCCTCA"
RIGHT_PRIMER = "GGAGGCCACACCCGCCACTCACCTG"

def read_metrics(path):
    with path.open() as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}

assigned = read_metrics(
    result_dir / "combined" / "ALL_ASSIGNED_library_metrics.tsv"
)
with_undetermined = read_metrics(
    result_dir / "combined" / "ALL_WITH_UNDETERMINED_library_metrics.tsv"
)

assert int(assigned["unique_design_sequences"]) == 20
assert int(assigned["detected_1x"]) == 18
assert int(assigned["dropout_sequences"]) == 2
assert int(with_undetermined["mapped_reads"]) > int(assigned["mapped_reads"])
assert int(assigned["phix_reads"]) == 0
assert int(with_undetermined["phix_reads"]) == 50
assert float(with_undetermined["mapped_fraction_non_phix_percent"]) > float(
    with_undetermined["mapped_fraction_percent"]
)
assert (result_dir / "RESULTS_TO_SHARE.txt").is_file()
assert (result_dir / "report.html").is_file()
assert (result_dir / "observed_target_sequences.tsv.gz").is_file()
assert (result_dir / "combined" / "variant_count_matrix.tsv").is_file()
assert (result_dir / "combined" / "variant_rpm_matrix.tsv").is_file()
assert (result_dir / "combined" / "all_sample_metrics.tsv").is_file()

oligo_reference = read_reference(
    result_dir.parent / "design_oligo.csv",
    "Variant_ID",
    "target_sequence",
    LEFT_PRIMER,
    RIGHT_PRIMER,
    "auto",
)
plain_reference = read_reference(
    result_dir.parent / "design.csv",
    "Variant_ID",
    "target_sequence",
    LEFT_PRIMER,
    RIGHT_PRIMER,
    "auto",
)
assert {x.sequence for x in oligo_reference.unique} == {
    x.sequence for x in plain_reference.unique
}

print("Synthetic integration test passed")
PY
