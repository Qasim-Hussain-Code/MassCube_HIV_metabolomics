#!/usr/bin/env bash
# Annotate both ion modes with the MassCube MS/MS database and the
# antiretroviral spectra of the GNPS Drug Library, export the feature
# tables, and normalise them. Runs only once both modes are aligned and
# their scan caches deleted, so the database never shares the disk with
# them (preregistration.md, sections 9 and 15).
#
# Usage: scripts/07_annotate.sh [--force] [--help]
#
# Writes results/database_version.tsv, database_load.tsv,
# drug_library_subset.tsv, and per mode <mode>_annotation.tsv.gz,
# <mode>_drug_hits.tsv, <mode>_features.tsv.gz, <mode>_detected.tsv.gz,
# <mode>_height_raw.tsv.gz, <mode>_height_normalised.tsv.gz,
# <mode>_area_raw.tsv.gz, <mode>_area_normalised.tsv.gz and, when under
# 20 MB, <mode>_features.msp.gz. The 3.3 GB database is deleted at the end.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=lib/common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --help|-h) sed -n '2,17p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done
STAGE="07_annotate"
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"

run_python_measured "$STAGE" "$(native_path "$ROOT/scripts/lib/annotate.py")" --root "$(native_path "$ROOT")"

stage_finish "$STAGE" "$MEASURED_RSS" "$(awk -F'\t' 'NR == 2 {printf "%s spectra, loaded in %s s, +%s MB", $2, $3, $6}' "$ROOT/results/database_load.tsv")"
