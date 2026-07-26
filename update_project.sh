#!/usr/bin/env bash
set -Eeuo pipefail

PACKAGE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${NGS_LIBRARYQC_PROJECT:-/data/user/MCET03/03_NGS/01_5UTR_Plasmid}"
TARGET_DIR="${PROJECT_DIR}/pipeline"
NGS_ROOT="$(cd -- "${PROJECT_DIR}/.." && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
STAGING_DIR="${PROJECT_DIR}/pipeline_update_staging_${STAMP}"
BACKUP_DIR="${NGS_ROOT}_pipeline_backup_${STAMP}"

if [[ ! -f "${PACKAGE_DIR}/amplicon_qc.py" ]] ||
   [[ ! -f "${PACKAGE_DIR}/run_project.sh" ]] ||
   [[ ! -f "${PACKAGE_DIR}/VERSION" ]]; then
    echo "This update package is incomplete." >&2
    exit 2
fi

if [[ ! -d "${PROJECT_DIR}" ]] || [[ ! -d "${TARGET_DIR}" ]]; then
    echo "Existing project pipeline was not found: ${TARGET_DIR}" >&2
    exit 2
fi

if pgrep -f '[a]mplicon_qc.py' >/dev/null 2>&1; then
    echo "amplicon_qc.py is running. Wait for it to stop before updating." >&2
    exit 2
fi

echo "Update source: ${PACKAGE_DIR}"
echo "Project: ${PROJECT_DIR}"
echo "Previous pipeline backup: ${BACKUP_DIR}"
echo "Raw data, reference, config, and results will not be changed."

while true; do
    read -r -p "Install this PR development update? [yes/no]: " ANSWER
    case "${ANSWER,,}" in
        yes|y) break ;;
        no|n)
            echo "Cancelled. No files were changed."
            exit 0
            ;;
        *) echo "Please answer yes or no." ;;
    esac
done

mkdir -- "${STAGING_DIR}"
cp -a -- "${PACKAGE_DIR}/." "${STAGING_DIR}/"

rollback() {
    if [[ ! -e "${TARGET_DIR}" && -d "${BACKUP_DIR}" ]]; then
        mv -- "${BACKUP_DIR}" "${TARGET_DIR}"
        echo "Update failed; the previous pipeline was restored." >&2
    fi
}
trap rollback ERR

mv -- "${TARGET_DIR}" "${BACKUP_DIR}"
mv -- "${STAGING_DIR}" "${TARGET_DIR}"
trap - ERR

echo "PR development update installed."
echo "Backup kept at: ${BACKUP_DIR}"
echo "Next: bash ${TARGET_DIR}/run_project.sh 1pct"
