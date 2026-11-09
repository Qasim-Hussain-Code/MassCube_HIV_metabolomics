#!/usr/bin/env bash
# Detect or accept the machine's threads, RAM and disk budget and write
# project.conf at the repository root. Every later stage sources that
# file, so the same scripts run on another machine without edits.
#
# Usage:
#   scripts/00_configure.sh [--threads N] [--ram GB] [--disk GB] [--yes] [--force] [--help]
#
# --disk is the ceiling for everything the project writes, environment
# and caches included; stages stop 0.5 GB below it. Without --yes the
# detected values are shown and must be confirmed. An existing
# project.conf is kept unless values or --force are given.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONF="$ROOT/project.conf"
THREADS=""
RAM_GB=""
DISK_GB=""
YES=0
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --ram)     RAM_GB="$2"; shift 2 ;;
        --disk)    DISK_GB="$2"; shift 2 ;;
        --yes)     YES=1; shift ;;
        --force)   FORCE=1; shift ;;
        --help|-h) sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [[ -f "$CONF" && $FORCE -eq 0 && -z "$THREADS$RAM_GB$DISK_GB" ]]; then
    echo "00_configure: project.conf already exists; skipping (use --force to rewrite)."
    exit 0
fi

detect_threads() {
    if command -v nproc >/dev/null 2>&1; then nproc
    elif command -v sysctl >/dev/null 2>&1; then sysctl -n hw.ncpu
    else echo 2; fi
}

detect_ram_gb() {
    # Git Bash has no /proc/meminfo, so Windows is asked directly.
    if [[ -r /proc/meminfo ]]; then
        awk '/MemTotal/ {printf "%d", $2 / 1024 / 1024 + 0.5}' /proc/meminfo
    elif command -v powershell.exe >/dev/null 2>&1; then
        powershell.exe -NoProfile -Command \
            "[int][math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB)" | tr -d '\r'
    elif command -v sysctl >/dev/null 2>&1; then
        echo $(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))
    else
        echo 8
    fi
}

detect_disk_gb() {
    df -Pk "$ROOT" | awk 'NR==2 {printf "%d", $4 / 1024 / 1024}'
}

DET_THREADS="$(detect_threads)"
DET_RAM="$(detect_ram_gb)"
DET_DISK="$(detect_disk_gb)"
THREADS="${THREADS:-$DET_THREADS}"
RAM_GB="${RAM_GB:-$DET_RAM}"
DISK_GB="${DISK_GB:-$DET_DISK}"

if (( DISK_GB > DET_DISK )); then
    echo "00_configure: --disk $DISK_GB GB exceeds the $DET_DISK GB free on this drive" >&2
    exit 1
fi

# Peak resident memory is planned against 14 GB on a 16 GB machine and
# proportionally below that elsewhere; installed RAM is not free RAM.
PLAN_RAM_GB="$(awk -v r="$RAM_GB" 'BEGIN {v = r - 2; if (v < 2) v = 2; printf "%d", v}')"
DISK_CEILING_GB="$(awk -v d="$DISK_GB" 'BEGIN {printf "%.1f", d - 0.5}')"

cat <<EOF
00_configure: detected
  threads : $DET_THREADS (using $THREADS)
  RAM     : $DET_RAM GB installed (using $RAM_GB, planning peak memory against $PLAN_RAM_GB)
  disk    : $DET_DISK GB free on the drive (budget $DISK_GB GB, stages stop at $DISK_CEILING_GB)
EOF

if [[ $YES -eq 0 ]]; then
    read -r -p "Write these values to project.conf? [y/N] " ans
    case "$ans" in y|Y|yes|YES) ;; *) echo "aborted"; exit 1 ;; esac
fi

cat > "$CONF" <<EOF
# Written by scripts/00_configure.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ)
# Sourced by every stage; re-run the script rather than editing by hand.
PROJECT_ROOT="$ROOT"
THREADS=$THREADS
RAM_GB=$RAM_GB
PLAN_RAM_GB=$PLAN_RAM_GB
DISK_GB=$DISK_GB
DISK_CEILING_GB=$DISK_CEILING_GB
FREE_FLOOR_GB=2
SEED=20261109
# Filled in by scripts/01_install.sh
CONDA_EXE_PATH=""
ENV_PREFIX=""
PYTHON=""
RSCRIPT=""
TRFP=""
EOF

echo "00_configure: wrote $CONF"
