#!/usr/bin/env bash
# Run the analysis in order. Every stage skips itself when its stamp in
# logs/ exists, so the script can be re-run after a failure and picks
# up where it stopped.
#
# Usage:
#   bash run_all.sh [--from <stage>] [--force] [--help]
#
#   --from <stage>   start at a stage number, for example --from 04
#   --force          re-run stages even if their stamp exists
#
# Each stage writes its elapsed time, peak resident memory and disk use
# to logs/<stage>.resources.tsv; the collated table is
# logs/resources_all.tsv.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FROM="00"
FORCE=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        --force) FORCE=(--force); shift ;;
        --help|-h) sed -n '2,14p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

FROM_NUM="${FROM:0:2}"
# Base 10 explicitly: bash reads "08" and "09" as invalid octal.
runs_from() { [[ $((10#$1)) -ge $((10#$FROM_NUM)) ]]; }

mkdir -p "$ROOT/logs" "$ROOT/results" "$ROOT/figures" "$ROOT/data"
T_ALL=$(date +%s)

stage() {
    # stage <number> <description> <command...>
    local num="$1" desc="$2"; shift 2
    if ! runs_from "$num"; then return 0; fi
    echo "==> [$num] $desc"
    "$@" 2>&1 | tee "$ROOT/logs/${num}_console.log"
    local rc=${PIPESTATUS[0]}
    if [[ $rc -ne 0 ]]; then echo "stage $num failed with exit code $rc" >&2; exit "$rc"; fi
}

stage 00 "configure" bash "$ROOT/scripts/00_configure.sh" --yes
stage 01 "install" bash "$ROOT/scripts/01_install.sh" "${FORCE[@]+"${FORCE[@]}"}"
stage 02 "study metadata" bash "$ROOT/scripts/02_download_metadata.sh" "${FORCE[@]+"${FORCE[@]}"}"
stage 03 "sample tables" bash "$ROOT/scripts/03_sample_table.sh" "${FORCE[@]+"${FORCE[@]}"}"

{
    printf "stage\telapsed_s\tpeak_rss_mb\tdisk_before_gb\tdisk_after_gb\tnote\n"
    for f in "$ROOT"/logs/*.resources.tsv; do
        [[ -f "$f" ]] && tail -n +2 "$f" | cut -f1-6
    done
} > "$ROOT/logs/resources_all.tsv"
echo "==> finished in $(( $(date +%s) - T_ALL )) s; per-stage cost in logs/resources_all.tsv"
