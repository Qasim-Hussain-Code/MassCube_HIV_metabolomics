#!/usr/bin/env bash
# Apply the pilot's m/z tolerance rule, project the full mode, and delete
# the pilot's raw files (preregistration.md, sections 7 and 16).
#
# Usage: scripts/05_pilot_rules.sh --mode pos|neg [--force] [--help]
#
# If the rule widens the tolerances, the parameter file is rewritten, the
# pilot's outputs are cleared and the pilot is processed again from its
# kept raw files, so no second download is needed; the rule is then
# recorded once more. Writes results/<mode>_mz_rule.tsv,
# results/<mode>_pilot_internal_standard.tsv and
# results/<mode>_pilot_projection.tsv.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=lib/common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
MODE=""
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        --help|-h) sed -n '2,13p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done
case "$MODE" in pos|neg) ;; *) echo "05_pilot_rules: --mode pos or --mode neg is required" >&2; exit 2 ;; esac
STAGE="05_rules_$MODE"
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"

RULES="$(native_path "$ROOT/scripts/lib/pilot_rules.py")"
set +e
run_python_measured "$STAGE" "$RULES" --root "$(native_path "$ROOT")" --mode "$MODE"
rc=$?
set -e
if [[ $rc -eq 10 ]]; then
    echo "$STAGE: tolerance widened; processing the pilot again from its kept raw files"
    any_python "$RULES" --root "$(native_path "$ROOT")" --mode "$MODE" --reset-pilot
    bash "$ROOT/scripts/04_process_mode.sh" --mode "$MODE" --pilot --force
    run_python_measured "$STAGE" "$RULES" --root "$(native_path "$ROOT")" --mode "$MODE"
elif [[ $rc -ne 0 ]]; then
    exit "$rc"
fi

stage_finish "$STAGE" "$MEASURED_RSS" "$(awk -F'\t' 'NR == 2 {print $NF}' "$ROOT/results/${MODE}_mz_rule.tsv")"
