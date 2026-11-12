#!/usr/bin/env bash
# Run the analysis in order. Every stage skips itself when its stamp in
# logs/ exists, so the script can be re-run after a failure and picks
# up where it stopped.
#
# Usage:
#   bash run_all.sh [--from <stage>] [--force] [--help]
#
#   --from <stage>   start at the first stage numbered <stage> or later and
#                    run everything after it, for example --from 05
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
        --help|-h) sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

FROM_NUM="${FROM:0:2}"
# --from N starts at the first stage numbered N or later and then runs
# everything after it in order; the ion modes repeat stage numbers, so a
# plain comparison would skip the second mode's early stages. Base 10
# explicitly: bash reads "08" and "09" as invalid octal.
STARTED=0
runs_from() {
    if [[ $STARTED -eq 1 ]] || [[ $((10#$1)) -ge $((10#$FROM_NUM)) ]]; then STARTED=1; return 0; fi
    return 1
}

mkdir -p "$ROOT/logs" "$ROOT/results" "$ROOT/figures" "$ROOT/data"
T_ALL=$(date +%s)

stage() {
    # stage <number> <description> <command...>
    local num="$1" desc="$2"; shift 2
    if ! runs_from "$num"; then return 0; fi
    echo "==> [$num] $desc"
    local log="${num}_${desc//[^A-Za-z0-9]/_}"
    "$@" 2>&1 | tee "$ROOT/logs/${log}_console.log"
    local rc=${PIPESTATUS[0]}
    if [[ $rc -ne 0 ]]; then echo "stage $num failed with exit code $rc" >&2; exit "$rc"; fi
}

stage 00 "configure" bash "$ROOT/scripts/00_configure.sh" --yes
stage 01 "install" bash "$ROOT/scripts/01_install.sh" "${FORCE[@]+"${FORCE[@]}"}"
stage 02 "study metadata" bash "$ROOT/scripts/02_download_metadata.sh" "${FORCE[@]+"${FORCE[@]}"}"
stage 03 "sample tables" bash "$ROOT/scripts/03_sample_table.sh" "${FORCE[@]+"${FORCE[@]}"}"
# One ion mode at a time, so that no more than one mode's scan caches are
# ever on disk: pilot, pilot rules, every file, alignment and clean-up,
# then the same for the other mode. The remaining files share number 05
# with the pilot rules they depend on.
for mode in pos neg; do
    stage 04 "$mode mode, pilot" bash "$ROOT/scripts/04_process_mode.sh" --mode "$mode" --pilot "${FORCE[@]+"${FORCE[@]}"}"
    stage 05 "$mode mode, pilot rules" bash "$ROOT/scripts/05_pilot_rules.sh" --mode "$mode" "${FORCE[@]+"${FORCE[@]}"}"
    stage 05 "$mode mode, all files" bash "$ROOT/scripts/04_process_mode.sh" --mode "$mode" "${FORCE[@]+"${FORCE[@]}"}"
    stage 06 "$mode mode, alignment" bash "$ROOT/scripts/06_align_mode.sh" --mode "$mode" "${FORCE[@]+"${FORCE[@]}"}"
done
stage 07 "annotation and normalisation" bash "$ROOT/scripts/07_annotate.sh" "${FORCE[@]+"${FORCE[@]}"}"

{
    printf "stage\telapsed_s\tpeak_rss_mb\tdisk_before_gb\tdisk_after_gb\tnote\n"
    for f in "$ROOT"/logs/*.resources.tsv; do
        [[ -f "$f" ]] && tail -n +2 "$f" | cut -f1-6
    done
} > "$ROOT/logs/resources_all.tsv"
echo "==> finished in $(( $(date +%s) - T_ALL )) s; per-stage cost in logs/resources_all.tsv"
