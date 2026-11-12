#!/usr/bin/env bash
# Injection QC, the run order test, the feature filters, the collapse to
# one feature per MassCube group, the antiretroviral positive control and
# the PCA with QCs (preregistration.md, sections 5, 6, 8 and 9). Writes the
# log2 matrices the statistics stage reads.
#
# Usage: scripts/08_qc_filters.sh [--force] [--help]
#
# Outputs are listed in scripts/lib/qc_filters.py.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=lib/common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --help|-h) sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done
STAGE="08_qc_filters"
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"

run_python_measured "$STAGE" "$(native_path "$ROOT/scripts/lib/qc_filters.py")" --root "$(native_path "$ROOT")"

stage_finish "$STAGE" "$MEASURED_RSS" "$(awk -F'\t' 'NR == 2 {print "positive control passes: " $1}' "$ROOT/results/antiretroviral_control_verdict.tsv")"
