#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
MODE="${1:-}"
OPTION="${2:-}"
CONFIG="${PROJECT_DIR}/config/libraryqc.ini"

usage() {
    echo "Usage:" >&2
    echo "  bash scripts/run.sh plot" >&2
    echo "  bash scripts/run.sh 1pct [--replace]" >&2
    echo "  bash scripts/run.sh full [--replace]" >&2
}

plot_result() {
    local result_dir="$1"
    if [[ ! -d "${result_dir}" ]]; then
        echo "Result folder not found: ${result_dir}" >&2
        exit 2
    fi
    if ! command -v Rscript >/dev/null 2>&1; then
        echo "Rscript is required to create the figures." >&2
        exit 2
    fi
    Rscript \
        "${SCRIPT_DIR}/plot_library_qc.R" \
        "${result_dir}" \
        "${result_dir}/figures"
    echo "Figures: ${result_dir}/figures"
}

if [[ "${MODE}" == "plot" ]]; then
    plot_result "${PROJECT_DIR}/results/1pct"
    exit 0
fi

if [[ "${MODE}" != "1pct" && "${MODE}" != "full" ]]; then
    usage
    exit 2
fi
if [[ -n "${OPTION}" && "${OPTION}" != "--replace" ]]; then
    usage
    exit 2
fi

if ! command -v flock >/dev/null 2>&1; then
    echo "Required command not found: flock" >&2
    exit 2
fi
LOCK_FILE="${PROJECT_DIR}/config/.analysis.lock"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
    echo "Another 1pct or full analysis is already running for this project." >&2
    exit 2
fi

if [[ ! -f "${CONFIG}" ]]; then
    python3 "${SCRIPT_DIR}/project_tools/configure_project.py"
fi

INPUT_DIR="${PROJECT_DIR}/raw_data/${MODE}"
OUTDIR="${PROJECT_DIR}/results/${MODE}"

if ! find "${INPUT_DIR}" -type f \
    \( -iname '*.fastq' -o -iname '*.fastq.gz' -o -iname '*.fq' -o -iname '*.fq.gz' \) \
    -print -quit | grep -q .; then
    echo "No FASTQ files found under: ${INPUT_DIR}" >&2
    exit 2
fi

if [[ -e "${OUTDIR}" ]]; then
    if [[ "${OPTION}" != "--replace" ]]; then
        echo "Result already exists: ${OUTDIR}" >&2
        echo "Use 'plot' to redraw it, or add --replace to rerun the analysis." >&2
        exit 2
    fi
    ARCHIVE_DIR="${PROJECT_DIR}/../archive"
    mkdir -p "${ARCHIVE_DIR}"
    BACKUP_BASE="${ARCHIVE_DIR}/previous_${MODE}_result"
    BACKUP="${BACKUP_BASE}"
    BACKUP_NUMBER=2
    while [[ -e "${BACKUP}" ]]; do
        BACKUP="${BACKUP_BASE}_${BACKUP_NUMBER}"
        ((BACKUP_NUMBER += 1))
    done
    mv -- "${OUTDIR}" "${BACKUP}"
    echo "Previous result moved to: ${BACKUP}"
fi

TMPDIR_RUN="$(mktemp -d)"
trap 'rm -rf -- "${TMPDIR_RUN}"' EXIT

echo "Mode: ${MODE}"
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

plot_result "${OUTDIR}" 2>&1 | tee -a "${OUTDIR}/run.log"

echo
echo "Completed [7/7]: ${OUTDIR}"
