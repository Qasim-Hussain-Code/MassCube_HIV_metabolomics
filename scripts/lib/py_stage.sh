#!/usr/bin/env bash
# py_stage.sh <stage> <script.py> [--force]: stamp and resource row around a Python script
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
STAGE="$1"
SCRIPT="$2"
FORCE=0
[[ "${3:-}" == "--force" ]] && FORCE=1
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"
run_python_measured "$STAGE" "$(native_path "$ROOT/$SCRIPT")" --root "$(native_path "$ROOT")"
stage_finish "$STAGE" "$MEASURED_RSS" "$SCRIPT"
