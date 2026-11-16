#!/usr/bin/env bash
# MTBLS9264 ISA-Tab, maf files, EBI checksums and raw file sizes. No spectra.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=lib/common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
STAGE="02_download_metadata"
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --help|-h) echo "usage: scripts/02_download_metadata.sh [--force]"; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"
# under 1 MB, so this only guards the free-space floor
require_space $((5 * 1024 * 1024))

run_python_measured "$STAGE" "$(native_path "$ROOT/scripts/lib/fetch_metadata.py")" --root "$(native_path "$ROOT")"

stage_finish "$STAGE" "$MEASURED_RSS" "$(awk -F'\t' 'NR == 2 {printf "%d raw files, %d bytes", $8, $9}' "$ROOT/results/mtbls9264_study.tsv")"
