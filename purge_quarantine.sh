#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
NGS_ROOT="$(cd -- "${PROJECT_DIR}/.." && pwd)"
PARENT_DIR="$(cd -- "${NGS_ROOT}/.." && pwd)"

mapfile -d '' QUARANTINES < <(
    find "${PARENT_DIR}" -maxdepth 1 -mindepth 1 -type d \
        -name "$(basename "${NGS_ROOT}")_cleanup_quarantine_*" -print0 | sort -z
)

if [[ "${#QUARANTINES[@]}" -eq 0 ]]; then
    echo "No cleanup quarantine was found."
    exit 0
fi

for TARGET in "${QUARANTINES[@]}"; do
    case "${TARGET}" in
        "${PARENT_DIR}/$(basename "${NGS_ROOT}")_cleanup_quarantine_"*) ;;
        *)
            echo "Refusing unexpected target: ${TARGET}" >&2
            exit 1
            ;;
    esac
    echo
    echo "Quarantine: ${TARGET}"
    du -sh -- "${TARGET}"
    read -r -p "Permanently delete this quarantine? [yes/no]: " ANSWER
    if [[ "${ANSWER,,}" == "yes" || "${ANSWER,,}" == "y" ]]; then
        rm -rf -- "${TARGET}"
        echo "Deleted: ${TARGET}"
    else
        echo "Kept: ${TARGET}"
    fi
done
