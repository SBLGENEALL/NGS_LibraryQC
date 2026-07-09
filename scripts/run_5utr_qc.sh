#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root:
#   bash scripts/run_5utr_qc.sh

python ngs_libraryqc.py --config configs/example_5utr_config.json
