#!/usr/bin/env bash
# Write the MassCube sample tables for both ion modes from the deposited
# ISA-Tab, and stop if any study sample lacks a group label.
#
# Usage: scripts/03_sample_table.sh [--force] [--help]
#
# Writes config/sample_table_pos.csv, config/sample_table_neg.csv and
# results/sample_counts.tsv. Injection order and batch are left empty;
# the conversion stage fills them from each file's acquisition time.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=lib/common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
STAGE="03_sample_table"
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --help|-h) sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"

run_python_measured "$STAGE" "$(native_path "$ROOT/scripts/lib/build_sample_tables.py")" --root "$(native_path "$ROOT")"

stage_finish "$STAGE" "$MEASURED_RSS" "$(awk -F'\t' 'NR > 1 {s += $4} END {print s}' "$ROOT/results/sample_counts.tsv") files across both modes"
