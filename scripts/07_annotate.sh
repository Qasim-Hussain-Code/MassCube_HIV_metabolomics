#!/usr/bin/env bash
# MS/MS annotation of both modes, then normalised tables. The 3.3 GB database
# only comes down after both modes' caches are gone (prereg section 14).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=lib/common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --help|-h) echo "usage: scripts/07_annotate.sh [--force]"; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done
STAGE="07_annotate"
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"

run_python_measured "$STAGE" "$(native_path "$ROOT/scripts/lib/annotate.py")" --root "$(native_path "$ROOT")"

stage_finish "$STAGE" "$MEASURED_RSS" "$(awk -F'\t' 'NR == 2 {printf "%s spectra, loaded in %s s, +%s MB", $2, $3, $6}' "$ROOT/results/database_load.tsv")"
