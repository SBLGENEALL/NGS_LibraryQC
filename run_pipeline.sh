#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/config.ini}"

if [[ ! -f "${CONFIG}" ]]; then
  echo "Configuration file not found: ${CONFIG}" >&2
  echo "Copy config/config_template.ini to config.ini and edit it." >&2
  exit 2
fi

python3 "${SCRIPT_DIR}/amplicon_qc.py" --config "${CONFIG}"

echo
echo "Analysis complete"
echo "Open the configured output directory and check:"
echo "  RESULTS_TO_SHARE.txt"
echo "  report.html"
