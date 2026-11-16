#!/usr/bin/env bash
# Align one mode once all its files are done, then free its scan caches.
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
        --help|-h) echo "usage: scripts/06_align_mode.sh --mode pos|neg [--force]"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done
case "$MODE" in pos|neg) ;; *) echo "06_align_mode: --mode pos or --mode neg is required" >&2; exit 2 ;; esac
STAGE="06_align_$MODE"
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"

run_python_measured "$STAGE" "$(native_path "$ROOT/scripts/lib/align_mode.py")" --root "$(native_path "$ROOT")" --mode "$MODE"

stage_finish "$STAGE" "$MEASURED_RSS" "$(awk -F'\t' 'NR == 2 {printf "%s aligned features from %s files", $3, $2}' "$ROOT/results/${MODE}_alignment_summary.tsv")"
