#!/usr/bin/env bash
# Shared shell helpers, sourced by every stage: configuration, stage
# stamps, the disk budget and resource logging.
#
# A stage sources this file, calls stage_skip_if_done, does its work,
# and finishes with stage_finish. Peak resident memory comes from
# scripts/lib/measure.py, which wraps the stage's heavy command; GNU
# time does not exist on Windows, and conda-forge's time package has no
# Windows build, so the measurement is done the same way everywhere.

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
CONF="$ROOT/project.conf"

load_conf() {
    if [[ ! -f "$CONF" ]]; then
        echo "project.conf not found; run scripts/00_configure.sh first" >&2
        exit 1
    fi
    # shellcheck source=/dev/null
    source "$CONF"
}

set_conf() {
    # Replace or append KEY="value" in project.conf.
    local key="$1" val="$2"
    if grep -q "^$key=" "$CONF"; then
        sed -i.bak "s|^$key=.*|$key=\"$val\"|" "$CONF" && rm -f "$CONF.bak"
    else
        echo "$key=\"$val\"" >> "$CONF"
    fi
}

stamp_file() { echo "$ROOT/logs/$1.done"; }

stage_skip_if_done() {
    # stage_skip_if_done <stage> <force 0|1>; exits 0 when the stamp exists
    local stage="$1" force="${2:-0}" f
    f="$(stamp_file "$stage")"
    if [[ -f "$f" && "$force" -eq 0 ]]; then
        echo "$stage: already completed on $(head -1 "$f"); skipping (use --force to redo)."
        exit 0
    fi
}

stage_mark_done() {
    mkdir -p "$ROOT/logs"
    date -u +%Y-%m-%dT%H:%M:%SZ > "$(stamp_file "$1")"
}

# ---- disk budget -----------------------------------------------------

any_python() {
    # The environment's Python once it exists; before that, the first
    # interpreter that actually runs ("python3" on Windows can be a Store
    # alias that is on PATH and does nothing).
    if [[ -x "$ROOT/env/python.exe" || -x "$ROOT/env/bin/python" ]]; then
        "$ROOT/scripts/lib/in_env.sh" python "$@"
        return
    fi
    local p
    for p in python python3; do
        if "$p" -c "import sys" >/dev/null 2>&1; then "$p" "$@"; return; fi
    done
    echo "no working Python found" >&2
    return 1
}

native_path() {
    case "$(uname -s)" in MINGW*|MSYS*|CYGWIN*) cygpath -w "$1" ;; *) echo "$1" ;; esac
}

project_bytes() {
    # Everything this project wrote, environment and caches included.
    # MSYS du opens every file, which over the environment took more than
    # half an hour on the build machine with the virus scanner running;
    # disk_usage.py reads sizes from the directory listings in seconds.
    any_python -I "$ROOT/scripts/lib/disk_usage.py" "$(native_path "$ROOT")"
}

project_bytes_du() {
    # The literal `du -sb .`, for the final check and the README. du
    # exits non-zero on a file locked by another process; the total is
    # still printed and still usable.
    { du -sb "$ROOT" 2>/dev/null || true; } | awk '{print $1 + 0}'
}

free_bytes() {
    df -Pk "$ROOT" | awk 'NR==2 {printf "%.0f", $4 * 1024}'
}

gb() { awk -v b="$1" 'BEGIN {printf "%.2f", b / 1024^3}'; }

require_space() {
    # require_space <incoming_bytes> [<current_project_bytes>]
    # Two limits, both checked before anything is fetched or written:
    # the project ceiling from project.conf, and a floor of free space
    # on the drive itself. The floor exists because C: on the build
    # machine is shared with a pagefile that Windows grows on its own.
    local incoming="$1" current="${2:-}" free ceiling floor
    [[ -n "$current" ]] || current="$(project_bytes)"
    free="$(free_bytes)"
    ceiling="$(awk -v g="$DISK_CEILING_GB" 'BEGIN {printf "%.0f", g * 1024^3}')"
    floor="$(awk -v g="$FREE_FLOOR_GB" 'BEGIN {printf "%.0f", g * 1024^3}')"
    if (( current + incoming > ceiling )); then
        echo "disk budget: project is $(gb "$current") GB and the next item adds $(gb "$incoming") GB, above the ${DISK_CEILING_GB} GB ceiling; stopping" >&2
        return 1
    fi
    if (( free - incoming < floor )); then
        echo "disk budget: $(gb "$free") GB free on the drive and the next item needs $(gb "$incoming") GB, which would leave less than ${FREE_FLOOR_GB} GB; stopping" >&2
        return 1
    fi
}

# ---- resource logging --------------------------------------------------

STAGE_T0=""
STAGE_DISK0=""

stage_start() {
    STAGE_T0="$(date +%s)"
    STAGE_DISK0="$(project_bytes)"
    echo "$1: project at $(gb "$STAGE_DISK0") GB, $(gb "$(free_bytes)") GB free on the drive"
}

stage_finish() {
    # stage_finish <stage> <peak_rss_mb or NA> <note>
    local stage="$1" rss="$2" note="${3:-}" t1 disk1
    t1="$(date +%s)"
    disk1="$(project_bytes)"
    mkdir -p "$ROOT/logs"
    printf "stage\telapsed_s\tpeak_rss_mb\tdisk_before_gb\tdisk_after_gb\tnote\tfinished\n%s\t%d\t%s\t%s\t%s\t%s\t%s\n" \
        "$stage" "$((t1 - STAGE_T0))" "$rss" "$(gb "$STAGE_DISK0")" "$(gb "$disk1")" "$note" \
        "$(date +%Y-%m-%dT%H:%M:%S%z)" > "$ROOT/logs/$stage.resources.tsv"
    stage_mark_done "$stage"
    echo "$stage: done in $((t1 - STAGE_T0)) s, peak RSS ${rss} MB, project now $(gb "$disk1") GB"
}

env_python_exe() {
    local p
    for p in "$ROOT/env/python.exe" "$ROOT/env/bin/python"; do
        if [[ -x "$p" ]]; then echo "$p"; return; fi
    done
}

env_rscript_exe() {
    local p
    for p in "$ROOT/env/Scripts/Rscript.exe" "$ROOT/env/Lib/R/bin/Rscript.exe" "$ROOT/env/bin/Rscript"; do
        if [[ -x "$p" ]]; then echo "$p"; return; fi
    done
}

MEASURED_RSS="NA"
run_r_measured() {
    # run_r_measured <label> <script.R> [args...]: as run_python_measured,
    # for an R script run by the environment's Rscript.
    local label="$1"; shift
    local rs rc=0
    rs="$(env_rscript_exe)"
    if [[ -z "$rs" ]]; then echo "no Rscript in env/; run scripts/01_install.sh" >&2; return 1; fi
    "$ROOT/scripts/lib/in_env.sh" python "$ROOT/scripts/lib/measure.py" \
        --label "$label" --log "$(native_path "$ROOT/logs/measure.tsv")" -- "$(native_path "$rs")" "$@" || rc=$?
    MEASURED_RSS="$(awk -F'\t' -v l="$label" '$1 == l {v = $4} END {print (v == "" ? "NA" : v)}' "$ROOT/logs/measure.tsv")"
    return $rc
}

run_python_measured() {
    # run_python_measured <label> <script> [args...]: run a Python script
    # from the environment under measure.py and leave its peak resident
    # memory (MB) in MEASURED_RSS. measure.py is handed the interpreter
    # itself, not in_env.sh, because Windows cannot start a shell script
    # as a process; in_env.sh around measure.py sets the PATH both inherit.
    local label="$1"; shift
    local py
    py="$(env_python_exe)"
    MEASURED_RSS="NA"
    if [[ -z "$py" ]]; then
        any_python "$@"
        return
    fi
    local rc=0
    "$ROOT/scripts/lib/in_env.sh" python "$ROOT/scripts/lib/measure.py" \
        --label "$label" --log "$(native_path "$ROOT/logs/measure.tsv")" -- "$(native_path "$py")" "$@" || rc=$?
    MEASURED_RSS="$(awk -F'\t' -v l="$label" '$1 == l {v = $4} END {print (v == "" ? "NA" : v)}' "$ROOT/logs/measure.tsv")"
    # The script's own exit status, so a caller can act on a specific code.
    return $rc
}
