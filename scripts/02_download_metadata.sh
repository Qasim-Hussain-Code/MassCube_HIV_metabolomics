#!/usr/bin/env bash
# Fetch the MTBLS9264 study metadata: the ISA-Tab investigation, study and
# assay tables, the authors' metabolite assignment files, the published
# checksums, and the size of every raw file. No spectrum is downloaded.
#
# Usage: scripts/02_download_metadata.sh [--force] [--help]
#
# Writes data/metadata/ (untracked), results/mtbls9264_files.tsv,
# results/mtbls9264_metadata_files.tsv and results/mtbls9264_study.tsv.

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
        --help|-h) sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"
# The metadata come to under 1 MB; the check is for the drive floor.
require_space $((5 * 1024 * 1024))

# Standard library only, so it also runs before the environment exists;
# memory is measured when the environment is there to measure it.
run_python_measured "$STAGE" "$(native_path "$ROOT/scripts/lib/fetch_metadata.py")" --root "$(native_path "$ROOT")"

stage_finish "$STAGE" "$MEASURED_RSS" "$(awk -F'\t' 'NR == 2 {printf "%d raw files, %d bytes", $8, $9}' "$ROOT/results/mtbls9264_study.tsv")"
