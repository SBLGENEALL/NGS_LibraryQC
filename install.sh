#!/usr/bin/env bash
set -Eeuo pipefail

PACKAGE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${NGS_LIBRARYQC_PROJECT:-/data/user/MCET03/03_NGS/01_5UTR_Plasmid}"
NGS_ROOT="$(cd -- "${PROJECT_DIR}/.." && pwd)"
ARCHIVE_ROOT="${NGS_ROOT}/archive"

if [[ ! -d "${PROJECT_DIR}" ]]; then
    echo "Project folder not found: ${PROJECT_DIR}" >&2
    exit 2
fi
if pgrep -f '[a]mplicon_qc.py' >/dev/null 2>&1; then
    echo "An analysis is running. Install after it finishes." >&2
    exit 2
fi
for required in \
    amplicon_qc.py \
    run_project.sh \
    scripts/plot_library_qc.R \
    project_tools/configure_project.py \
    project_tools/run_with_resources.py \
    references/NC_001422.1_phiX174.fasta; do
    if [[ ! -f "${PACKAGE_DIR}/${required}" ]]; then
        echo "Package file missing: ${required}" >&2
        exit 2
    fi
done

echo "Project: ${PROJECT_DIR}"
echo "Final folders: raw_data, reference, scripts, config, results"
echo "Old extra items will be moved to: ${ARCHIVE_ROOT}"
read -r -p "Install the clean pipeline? [yes/no]: " ANSWER
case "${ANSWER,,}" in
    yes|y) ;;
    *)
        echo "Cancelled. No files were changed."
        exit 0
        ;;
esac

mkdir -p \
    "${PROJECT_DIR}/raw_data/1pct" \
    "${PROJECT_DIR}/raw_data/full" \
    "${PROJECT_DIR}/reference" \
    "${PROJECT_DIR}/config" \
    "${PROJECT_DIR}/results" \
    "${ARCHIVE_ROOT}"

unique_path() {
    local base="$1"
    local candidate="${base}"
    local number=2
    while [[ -e "${candidate}" ]]; do
        candidate="${base}_${number}"
        ((number += 1))
    done
    printf '%s\n' "${candidate}"
}

move_children() {
    local source="$1"
    local destination="$2"
    [[ -d "${source}" ]] || return 0
    mkdir -p "${destination}"
    while IFS= read -r -d '' item; do
        local target="${destination}/$(basename -- "${item}")"
        if [[ -e "${target}" ]]; then
            echo "Name conflict; installation stopped before overwriting: ${target}" >&2
            exit 2
        fi
        mv -- "${item}" "${target}"
    done < <(find "${source}" -mindepth 1 -maxdepth 1 -print0)
}

# Put the extracted 1% FASTQs under raw_data/1pct.
if [[ -d "${PROJECT_DIR}/subset_1pct" ]]; then
    move_children "${PROJECT_DIR}/subset_1pct" "${PROJECT_DIR}/raw_data/1pct"
    rmdir "${PROJECT_DIR}/subset_1pct"
fi

# Put all original FASTQ folders/files under raw_data/full.
while IFS= read -r -d '' item; do
    target="${PROJECT_DIR}/raw_data/full/$(basename -- "${item}")"
    if [[ -e "${target}" ]]; then
        echo "Name conflict; installation stopped before overwriting: ${target}" >&2
        exit 2
    fi
    mv -- "${item}" "${target}"
done < <(
    find "${PROJECT_DIR}/raw_data" -mindepth 1 -maxdepth 1 \
        ! -name 1pct ! -name full -print0
)

# Preserve the newest completed legacy 1% result at the fixed final path.
if [[ ! -e "${PROJECT_DIR}/results/1pct" ]]; then
    LEGACY_RESULT="$(
        find "${PROJECT_DIR}/results" -mindepth 1 -maxdepth 1 -type d \
            ! -name 1pct ! -name full \
            \( -name '*1pct*' -o -name '*1percent*' -o -name '*1_percent*' \) \
            -printf '%T@ %p\n' |
        sort -nr |
        while read -r _ path; do
            if [[ -f "${path}/combined/ALL_ASSIGNED_variant_counts.csv" ]] ||
               [[ -f "${path}/combined/ALL_ASSIGNED_variant_counts.tsv" ]]; then
                printf '%s\n' "${path}"
                break
            fi
        done
    )"
    if [[ -n "${LEGACY_RESULT}" ]]; then
        mv -- "${LEGACY_RESULT}" "${PROJECT_DIR}/results/1pct"
    fi
fi

OLD_RESULTS="$(unique_path "${ARCHIVE_ROOT}/old_results")"
mkdir -p "${OLD_RESULTS}"
while IFS= read -r -d '' item; do
    mv -- "${item}" "${OLD_RESULTS}/"
done < <(
    find "${PROJECT_DIR}/results" -mindepth 1 -maxdepth 1 \
        ! -name 1pct ! -name full -print0
)
rmdir "${OLD_RESULTS}" 2>/dev/null || true

# Keep only the five requested top-level project folders.
OLD_ITEMS="$(unique_path "${ARCHIVE_ROOT}/old_project_items")"
mkdir -p "${OLD_ITEMS}"
while IFS= read -r -d '' item; do
    mv -- "${item}" "${OLD_ITEMS}/"
done < <(
    find "${PROJECT_DIR}" -mindepth 1 -maxdepth 1 \
        ! -name raw_data \
        ! -name reference \
        ! -name scripts \
        ! -name config \
        ! -name results \
        -print0
)

if [[ -d "${PROJECT_DIR}/scripts" ]]; then
    mv -- "${PROJECT_DIR}/scripts" "${OLD_ITEMS}/scripts"
fi
rmdir "${OLD_ITEMS}" 2>/dev/null || true

NEW_SCRIPTS="${PROJECT_DIR}/scripts"
mkdir -p \
    "${NEW_SCRIPTS}/project_tools" \
    "${NEW_SCRIPTS}/references"

install -m 0644 "${PACKAGE_DIR}/amplicon_qc.py" "${NEW_SCRIPTS}/amplicon_qc.py"
install -m 0755 "${PACKAGE_DIR}/run_project.sh" "${NEW_SCRIPTS}/run.sh"
install -m 0755 "${PACKAGE_DIR}/scripts/plot_library_qc.R" \
    "${NEW_SCRIPTS}/plot_library_qc.R"
install -m 0644 "${PACKAGE_DIR}/project_tools/configure_project.py" \
    "${NEW_SCRIPTS}/project_tools/configure_project.py"
install -m 0644 "${PACKAGE_DIR}/project_tools/run_with_resources.py" \
    "${NEW_SCRIPTS}/project_tools/run_with_resources.py"
install -m 0644 "${PACKAGE_DIR}/references/NC_001422.1_phiX174.fasta" \
    "${NEW_SCRIPTS}/references/NC_001422.1_phiX174.fasta"

echo
echo "Clean pipeline installed."
echo "Replot current 1% result:"
echo "  bash ${NEW_SCRIPTS}/run.sh plot"
echo "Run the final pipeline on 1% FASTQs:"
echo "  bash ${NEW_SCRIPTS}/run.sh 1pct --replace"
