#!/usr/bin/env bash
set -euo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd -- "${TEST_DIR}/.." && pwd)"

cd "${PIPELINE_DIR}"
python3 -m unittest tests/test_sample_sheet.py tests/test_project_tools.py

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
        return {row["metric"]: row["value"] for row in csv.DictReader(handle)}

assigned = read_metrics(
    result_dir / "combined" / "ALL_ASSIGNED_library_metrics.csv"
)
with_undetermined = read_metrics(
    result_dir / "combined" / "ALL_WITH_UNDETERMINED_library_metrics.csv"
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
assert (result_dir / "observed_target_sequences.csv.gz").is_file()
assert (result_dir / "combined" / "variant_count_matrix.csv").is_file()
assert (result_dir / "combined" / "variant_rpm_matrix.csv").is_file()
assert (result_dir / "combined" / "all_sample_metrics.csv").is_file()
assert not list(result_dir.rglob("*.tsv"))

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

if command -v Rscript >/dev/null 2>&1; then
    Rscript "${PIPELINE_DIR}/scripts/plot_library_qc.R" \
        "${TEST_DIR}/synthetic/results" \
        "${TEST_DIR}/synthetic/results/figures"
    test -f "${TEST_DIR}/synthetic/results/figures/01_coverage_classes.png"
    test -f "${TEST_DIR}/synthetic/results/figures/library_qc_figures.pdf"

    python3 - "${TEST_DIR}/synthetic/plot_fixture" <<'PY'
import csv
import sys
from pathlib import Path

result = Path(sys.argv[1])
combined = result / "combined"
combined.mkdir(parents=True, exist_ok=True)
groups = [(11, 0, 0), (113, 1, 9), (864, 10, 99), (1013, 100, 599)]
rows = []
number = 0
for size, low, high in groups:
    span = high - low + 1
    for index in range(size):
        number += 1
        total = low + (index % span)
        rows.append(
            {
                "variant_ids": f"UTR_{number:04d}",
                "length": 50 + (number % 81),
                "gc_percent": 30 + (number % 41),
                "exact_count": total,
                "near_count": 0,
                "total_count": total,
            }
        )
with (combined / "ALL_ASSIGNED_variant_counts.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0])
    writer.writeheader()
    writer.writerows(rows)
PY

    Rscript "${PIPELINE_DIR}/scripts/plot_library_qc.R" \
        "${TEST_DIR}/synthetic/plot_fixture" \
        "${TEST_DIR}/synthetic/plot_fixture/figures"

    python3 - "${TEST_DIR}/synthetic/plot_fixture/figures" <<'PY'
import csv
import sys
from pathlib import Path

figures = Path(sys.argv[1])
with (figures / "coverage_class_summary.csv").open(encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
assert [int(row["variants"]) for row in rows] == [11, 113, 864, 1013]
assert (figures / "01_coverage_classes.png").is_file()
assert (figures / "library_qc_figures.pdf").is_file()
print("R plotting fixture passed")
PY
fi
