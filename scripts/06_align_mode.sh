#!/usr/bin/env bash
# Align one ion mode with MassCube once all of its files are processed:
# run order, the retention time and drift rules, alignment with gap
# filling, per-injection metrics, then deletion of the scan caches and
# compression of the per-file tables so the other mode fits on disk.
#
# Usage: scripts/06_align_mode.sh --mode pos|neg [--force] [--help]
#
# Writes results/<mode>_run_order.tsv, <mode>_rt_rule.tsv,
# <mode>_drift_rule.tsv, <mode>_injection_metrics.tsv and
# <mode>_alignment_summary.tsv; fills injection_order and batch in
# config/sample_table_<mode>.csv; may switch correct_rt and
# signal_normalization on in config/masscube_parameters_<mode>.csv; and
# leaves the aligned features in work/<mode>/project_files/.

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
        --help|-h) sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done
case "$MODE" in pos|neg) ;; *) echo "06_align_mode: --mode pos or --mode neg is required" >&2; exit 2 ;; esac
STAGE="06_align_$MODE"
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"

run_python_measured "$STAGE" "$(native_path "$ROOT/scripts/lib/align_mode.py")" --root "$(native_path "$ROOT")" --mode "$MODE"

stage_finish "$STAGE" "$MEASURED_RSS" "$(awk -F'\t' 'NR == 2 {printf "%s aligned features from %s files", $3, $2}' "$ROOT/results/${MODE}_alignment_summary.tsv")"
