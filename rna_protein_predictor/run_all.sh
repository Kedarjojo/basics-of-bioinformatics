#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"

skip_install="${RNA_PROTEIN_SKIP_INSTALL:-0}"
for argument in "$@"; do
    if [[ "$argument" == "--dry-run" ]]; then
        skip_install="1"
    fi
done

if [[ "$skip_install" != "1" ]]; then
    python3 -m pip install -r requirements.txt
fi
python3 scripts/run_pipeline.py "$@"
