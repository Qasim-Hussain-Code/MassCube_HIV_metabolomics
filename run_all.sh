#!/usr/bin/env bash
# Whole analysis, stage by stage. Finished stages leave a stamp in logs/ and get skipped on a rerun.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FROM="00"
FORCE=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        --force) FORCE=(--force); shift ;;
        --help|-h) echo "usage: bash run_all.sh [--from NN] [--force]   (--from 05 starts at stage 05 and runs the rest)"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

FROM_NUM="${FROM:0:2}"
# The two ion modes reuse stage numbers, so once --from has matched, everything after it runs.
# 10# because bash reads 08 and 09 as octal.
STARTED=0
runs_from() {
    if [[ $STARTED -eq 1 ]] || [[ $((10#$1)) -ge $((10#$FROM_NUM)) ]]; then STARTED=1; return 0; fi
    return 1
}

mkdir -p "$ROOT/logs" "$ROOT/results" "$ROOT/figures" "$ROOT/data"
T_ALL=$(date +%s)

stage() {   # stage NN "description" command...
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
# One mode at a time: both modes' scan caches together don't fit in 10 GB.
for mode in pos neg; do
    stage 04 "$mode mode, pilot" bash "$ROOT/scripts/04_process_mode.sh" --mode "$mode" --pilot "${FORCE[@]+"${FORCE[@]}"}"
    stage 05 "$mode mode, pilot rules" bash "$ROOT/scripts/05_pilot_rules.sh" --mode "$mode" "${FORCE[@]+"${FORCE[@]}"}"
    stage 05 "$mode mode, all files" bash "$ROOT/scripts/04_process_mode.sh" --mode "$mode" "${FORCE[@]+"${FORCE[@]}"}"
    stage 06 "$mode mode, alignment" bash "$ROOT/scripts/06_align_mode.sh" --mode "$mode" "${FORCE[@]+"${FORCE[@]}"}"
done
stage 07 "annotation and normalisation" bash "$ROOT/scripts/07_annotate.sh" "${FORCE[@]+"${FORCE[@]}"}"
stage 08 "injection QC, filters, positive control" bash "$ROOT/scripts/08_qc_filters.sh" "${FORCE[@]+"${FORCE[@]}"}"
stage 09 "statistics" bash "$ROOT/scripts/lib/r_stage.sh" 09_statistics scripts/09_statistics.R "${FORCE[@]+"${FORCE[@]}"}"
# exploratory, preregistration section 12
stage 10 "pathway over-representation" bash "$ROOT/scripts/lib/py_stage.sh" 10_pathways scripts/lib/pathways.py "${FORCE[@]+"${FORCE[@]}"}"
stage 10 "classifier" bash "$ROOT/scripts/lib/py_stage.sh" 10_classifier scripts/lib/classifier.py "${FORCE[@]+"${FORCE[@]}"}"
stage 10 "PERMANOVA" bash "$ROOT/scripts/lib/r_stage.sh" 10_permanova scripts/10_permanova.R "${FORCE[@]+"${FORCE[@]}"}"
stage 11 "figures" bash "$ROOT/scripts/lib/r_stage.sh" 11_figures scripts/11_figures.R "${FORCE[@]+"${FORCE[@]}"}"

{
    printf "stage\telapsed_s\tpeak_rss_mb\tdisk_before_gb\tdisk_after_gb\tnote\n"
    for f in "$ROOT"/logs/*.resources.tsv; do
        [[ -f "$f" ]] && tail -n +2 "$f" | cut -f1-6
    done
} > "$ROOT/logs/resources_all.tsv"
echo "==> finished in $(( $(date +%s) - T_ALL )) s; per-stage cost in logs/resources_all.tsv"
