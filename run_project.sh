#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
MODE="${1:-}"
CONFIG="${PROJECT_DIR}/config/libraryqc.ini"
VERSION="$(tr -d '[:space:]' < "${SCRIPT_DIR}/VERSION")"

if [[ "${MODE}" != "1pct" && "${MODE}" != "full" ]]; then
    echo "Usage: bash pipeline/run_project.sh 1pct|full" >&2
    exit 2
fi

if [[ ! -f "${CONFIG}" ]]; then
    python3 "${SCRIPT_DIR}/project_tools/configure_project.py"
fi

if [[ "${MODE}" == "1pct" ]]; then
    INPUT_DIR="${PROJECT_DIR}/subset_1pct"
else
    INPUT_DIR="${PROJECT_DIR}/raw_data"
fi

if ! find "${INPUT_DIR}" -type f \( -iname '*.fastq' -o -iname '*.fastq.gz' -o -iname '*.fq' -o -iname '*.fq.gz' \) -print -quit | grep -q .; then
    echo "No FASTQ files found under: ${INPUT_DIR}" >&2
    exit 2
fi

DATE="$(date +%Y%m%d)"
BASE_OUT="${PROJECT_DIR}/results/${DATE}_${MODE}"
OUTDIR="${BASE_OUT}"
RERUN=2
while [[ -e "${OUTDIR}" ]]; do
    OUTDIR="${BASE_OUT}_run${RERUN}"
    ((RERUN += 1))
done

TMPDIR_RUN="$(mktemp -d)"
trap 'rm -rf -- "${TMPDIR_RUN}"' EXIT

echo "Mode: ${MODE}"
echo "Pipeline version: ${VERSION}"
echo "Input: ${INPUT_DIR}"
echo "Output: ${OUTDIR}"
echo

set +e
python3 "${SCRIPT_DIR}/project_tools/run_with_resources.py" \
    --resource-file "${TMPDIR_RUN}/resource_usage.txt" \
    --log-file "${TMPDIR_RUN}/run.log" \
    -- \
    python3 "${SCRIPT_DIR}/amplicon_qc.py" \
        --config "${CONFIG}" \
        --input-dir "${INPUT_DIR}" \
        --outdir "${OUTDIR}"
STATUS="$?"
set -e

mkdir -p "${OUTDIR}"
mv "${TMPDIR_RUN}/resource_usage.txt" "${OUTDIR}/resource_usage.txt"
mv "${TMPDIR_RUN}/run.log" "${OUTDIR}/run.log"

if [[ "${STATUS}" -ne 0 ]]; then
    echo "Analysis failed with exit code ${STATUS}." >&2
    echo "Check: ${OUTDIR}/run.log" >&2
    exit "${STATUS}"
fi

if ! grep -Fq '[7/7] Creating HTML and shareable summary' "${OUTDIR}/run.log" ||
   [[ ! -f "${OUTDIR}/manifest.json" ]] ||
   [[ ! -f "${OUTDIR}/report.html" ]] ||
   [[ ! -f "${OUTDIR}/RESULTS_TO_SHARE.txt" ]]; then
    echo "Analysis did not pass the [7/7] completion check." >&2
    exit 1
fi

if command -v Rscript >/dev/null 2>&1; then
    Rscript "${SCRIPT_DIR}/scripts/plot_library_qc.R" "${OUTDIR}" "${OUTDIR}/figures" \
        2>&1 | tee -a "${OUTDIR}/run.log"
else
    echo "Rscript was not found; R figures were skipped." | tee -a "${OUTDIR}/run.log"
fi

if [[ "${MODE}" == "1pct" ]] &&
   find "${PROJECT_DIR}/raw_data" -type f \( -iname '*.fastq' -o -iname '*.fastq.gz' -o -iname '*.fq' -o -iname '*.fq.gz' \) -print -quit | grep -q .; then
    python3 "${SCRIPT_DIR}/project_tools/estimate_full_run.py" \
        "${OUTDIR}" "${PROJECT_DIR}/subset_1pct" "${PROJECT_DIR}/raw_data" \
        | tee -a "${OUTDIR}/run.log"
fi

echo
echo "Completed [7/7]: ${OUTDIR}"
echo "Keep the cleanup quarantine until this result has been reviewed."
