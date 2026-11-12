#!/usr/bin/env bash
# Run one R stage with the same bookkeeping as the shell stages: skip when
# its stamp exists, measure elapsed time, peak resident memory and disk,
# and write logs/<stage>.resources.tsv and the stamp.
#
# Usage: scripts/lib/r_stage.sh <stage> <script.R> [--force]
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
run_r_measured "$STAGE" "$(native_path "$ROOT/$SCRIPT")"
stage_finish "$STAGE" "$MEASURED_RSS" "$SCRIPT"
