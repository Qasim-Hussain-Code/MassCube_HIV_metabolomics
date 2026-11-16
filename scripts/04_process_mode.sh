#!/usr/bin/env bash
# download -> mzML -> MassCube, one raw file at a time; resumes where it stopped
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=lib/common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
MODE=""
PILOT=0
LIMIT=0
CONNECTIONS=8
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --pilot) PILOT=1; shift ;;
        --limit) LIMIT="$2"; shift 2 ;;
        --connections) CONNECTIONS="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        --help|-h) echo "usage: scripts/04_process_mode.sh --mode pos|neg [--pilot] [--limit N] [--connections K] [--force]"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done
case "$MODE" in pos|neg) ;; *) echo "04_process_mode: --mode pos or --mode neg is required" >&2; exit 2 ;; esac

if [[ $PILOT -eq 1 ]]; then STAGE="04_pilot_$MODE"; else STAGE="04_process_$MODE"; fi
# the full mode needs the tolerances the pilot settles (prereg section 7)
if [[ $PILOT -eq 0 && $LIMIT -eq 0 && ! -f "$(stamp_file "05_rules_$MODE")" ]]; then
    echo "04_process_mode: the pilot rules for $MODE have not been applied; run scripts/05_pilot_rules.sh --mode $MODE first" >&2
    exit 1
fi
# --limit is for testing: no stamp read or written
if [[ $LIMIT -eq 0 ]]; then stage_skip_if_done "$STAGE" "$FORCE"; fi
stage_start "$STAGE"

ARGS=(--root "$(native_path "$ROOT")" --mode "$MODE" --connections "$CONNECTIONS")
[[ $PILOT -eq 1 ]] && ARGS+=(--pilot)
[[ $LIMIT -gt 0 ]] && ARGS+=(--limit "$LIMIT")
run_python_measured "$STAGE" "$(native_path "$ROOT/scripts/lib/process_mode.py")" "${ARGS[@]}"

n="$(awk -F'\t' -v m="$([[ $MODE == pos ]] && echo positive || echo negative)" '$2 == m' "$ROOT/results/conversion_manifest.tsv" | wc -l)"
if [[ $LIMIT -gt 0 ]]; then
    echo "$STAGE: test run of $LIMIT file(s) finished; $n file(s) of this mode in the manifest; no stamp written"
    exit 0
fi
stage_finish "$STAGE" "$MEASURED_RSS" "$n files of this mode in the manifest"
