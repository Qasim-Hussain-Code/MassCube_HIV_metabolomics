#!/usr/bin/env bash
# Build the software environment inside the project and record what was
# installed.
#
# The conda environment goes to env/ and every cache (conda packages,
# pip, R's downloads) to .cache/, all untracked, so that the disk budget
# in project.conf measures them. An existing conda is used if one is on
# PATH or in the usual places; otherwise Miniforge is installed into
# miniforge/. ThermoRawFileParser 1.4.5 comes from its GitHub release on
# Windows, where bioconda has no build, and from bioconda elsewhere. The
# R packages are installed by BiocManager inside the conda R. Caches are
# emptied at the end.
#
# Usage: scripts/01_install.sh [--force] [--help]
#
# Writes config/python_packages.tsv, config/r_packages.tsv,
# config/environment.lock.yml, config/tool_versions.tsv,
# logs/sessionInfo_install.txt, and on first retrieval the
# ThermoRawFileParser checksum in config/external_checksums.tsv.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=lib/common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
STAGE="01_install"
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --help|-h) sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done
stage_skip_if_done "$STAGE" "$FORCE"
stage_start "$STAGE"

LOG="$ROOT/logs/01_install.log"
ENV="$ROOT/env"
TRFP_VERSION="1.4.5"
TRFP_ZIP_URL="https://github.com/CompOmics/ThermoRawFileParser/releases/download/v${TRFP_VERSION}/ThermoRawFileParser${TRFP_VERSION}.zip"
mkdir -p "$ROOT/logs" "$ROOT/.cache/tmp"

case "$(uname -s)" in MINGW*|MSYS*|CYGWIN*) WINDOWS=1 ;; *) WINDOWS=0 ;; esac
winpath() { if [[ $WINDOWS -eq 1 ]]; then cygpath -w "$1"; else echo "$1"; fi; }

CONDA_PKGS_DIRS="$(winpath "$ROOT/.cache/conda_pkgs")"
PIP_CACHE_DIR="$(winpath "$ROOT/.cache/pip")"
export CONDA_PKGS_DIRS PIP_CACHE_DIR
export CONDA_NOTIFY_OUTDATED_CONDA=false

find_conda() {
    local c
    for c in "$(command -v mamba || true)" "$(command -v conda || true)" \
             "$HOME/miniforge3/condabin/conda" "$HOME/miniforge3/Scripts/conda.exe" \
             "$HOME/miniconda3/condabin/conda" "$HOME/miniconda3/Scripts/conda.exe" \
             "$HOME/anaconda3/Scripts/conda.exe" "$ROOT/miniforge/condabin/conda" \
             "$ROOT/miniforge/Scripts/conda.exe"; do
        if [[ -n "$c" && -x "$c" ]]; then echo "$c"; return; fi
    done
}

CONDA="$(find_conda)"
if [[ -z "$CONDA" ]]; then
    # About 90 MB of installer; Miniforge itself is about 0.4 GB once unpacked.
    require_space $((600 * 1024 * 1024))
    echo "01_install: no conda found; installing Miniforge into miniforge/" | tee -a "$LOG"
    if [[ $WINDOWS -eq 1 ]]; then
        curl -fL --retry 5 -o "$ROOT/.cache/Miniforge3.exe" \
            "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe"
        "$ROOT/.cache/Miniforge3.exe" /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S "/D=$(winpath "$ROOT/miniforge")"
        CONDA="$ROOT/miniforge/Scripts/conda.exe"
    else
        curl -fL --retry 5 -o "$ROOT/.cache/Miniforge3.sh" \
            "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
        bash "$ROOT/.cache/Miniforge3.sh" -b -p "$ROOT/miniforge"
        CONDA="$ROOT/miniforge/bin/conda"
    fi
    rm -f "$ROOT/.cache/Miniforge3.exe" "$ROOT/.cache/Miniforge3.sh"
fi
echo "01_install: using $CONDA ($("$CONDA" --version 2>&1 | tr -d '\r'))" | tee -a "$LOG"

if [[ ! -d "$ENV/conda-meta" ]]; then
    # Room for the environment and its package cache together, until the
    # cache is cleaned below; the measured size is in the resource log.
    require_space $((5 * 1024 * 1024 * 1024))
    echo "01_install: creating the environment in env/ (downloads about 1 GB)" | tee -a "$LOG"
    "$CONDA" env create -y -p "$(winpath "$ENV")" -f "$(winpath "$ROOT/config/environment.yml")" 2>&1 | tee -a "$LOG"
else
    echo "01_install: env/ already exists" | tee -a "$LOG"
fi
PY=("$ROOT/scripts/lib/in_env.sh" python)
RS=("$ROOT/scripts/lib/in_env.sh" Rscript)
"${PY[@]}" -c "import masscube, ms_entropy, psutil" 2>&1 | tee -a "$LOG"

# ThermoRawFileParser. The 1.4.5 zip on GitHub publishes no checksum, so
# the sha256 seen at first retrieval is recorded and every later
# retrieval is checked against it.
CHECKSUMS="$ROOT/config/external_checksums.tsv"
[[ -f "$CHECKSUMS" ]] || printf "file\tsha256\tsize\tsource\trecorded\n" > "$CHECKSUMS"
if [[ $WINDOWS -eq 1 ]]; then
    TRFP_DIR="$ENV/opt/ThermoRawFileParser"
    TRFP="$TRFP_DIR/ThermoRawFileParser.exe"
    TRFP_SOURCE="GitHub release v$TRFP_VERSION"
    zip="$ROOT/.cache/ThermoRawFileParser${TRFP_VERSION}.zip"
    if [[ ! -x "$TRFP" ]]; then
        curl -fL --retry 5 --retry-delay 10 -o "$zip" "$TRFP_ZIP_URL"
        sha="$(sha256sum "$zip" | cut -d' ' -f1)"
        known="$(awk -F'\t' -v f="$(basename "$zip")" '$1 == f {print $2}' "$CHECKSUMS")"
        if [[ -z "$known" ]]; then
            printf "%s\t%s\t%s\t%s\t%s\n" "$(basename "$zip")" "$sha" "$(stat -c %s "$zip")" "$TRFP_ZIP_URL" \
                "$(date +%Y-%m-%dT%H:%M:%S%z)" >> "$CHECKSUMS"
        elif [[ "$known" != "$sha" ]]; then
            echo "01_install: ThermoRawFileParser zip sha256 $sha differs from the recorded $known" >&2
            exit 1
        fi
        mkdir -p "$TRFP_DIR"
        "${PY[@]}" -I -m zipfile -e "$(winpath "$zip")" "$(winpath "$TRFP_DIR")"
        rm -f "$zip"
    fi
    TRFP_SHA="sha256:$(awk -F'\t' -v f="ThermoRawFileParser${TRFP_VERSION}.zip" '$1 == f {print $2}' "$CHECKSUMS")"
else
    TRFP="$ENV/bin/thermorawfileparser"
    TRFP_SOURCE="bioconda thermorawfileparser=$TRFP_VERSION"
    TRFP_SHA="NA"
    [[ -x "$TRFP" ]] || "$CONDA" install -y -p "$ENV" -c conda-forge -c bioconda "thermorawfileparser=$TRFP_VERSION" 2>&1 | tee -a "$LOG"
fi
"$TRFP" --version 2>&1 | tr -d '\r' | tee -a "$LOG" || true

echo "01_install: installing R packages" | tee -a "$LOG"
TMP="$(winpath "$ROOT/.cache/tmp")" TEMP="$(winpath "$ROOT/.cache/tmp")" TMPDIR="$(winpath "$ROOT/.cache/tmp")" \
    "${RS[@]}" "$ROOT/scripts/lib/install_r_packages.R" \
    "$(winpath "$ROOT/config/r_packages.tsv")" "$(winpath "$ROOT/logs/sessionInfo_install.txt")" 2>&1 | tee -a "$LOG"
R_VERSION="$("${RS[@]}" -e 'cat(R.version.string)' 2>/dev/null | tr -d '\r')"

"${PY[@]}" "$ROOT/scripts/lib/record_versions.py" --root "$(winpath "$ROOT")" \
    --trfp-version "$TRFP_VERSION" --trfp-source "$TRFP_SOURCE" --trfp-sha256 "$TRFP_SHA" --r-version "$R_VERSION" 2>&1 | tee -a "$LOG"
# The lock file names exact builds; the local prefix line is dropped
# because it is a path on the build machine.
"$CONDA" env export -p "$(winpath "$ENV")" 2>/dev/null | tr -d '\r' | grep -v '^prefix:' > "$ROOT/config/environment.lock.yml"

set_conf CONDA_EXE_PATH "$CONDA"
set_conf ENV_PREFIX "$ENV"
set_conf PYTHON "$ROOT/scripts/lib/in_env.sh python"
set_conf RSCRIPT "$ROOT/scripts/lib/in_env.sh Rscript"
set_conf TRFP "$TRFP"

echo "01_install: cleaning caches" | tee -a "$LOG"
"$CONDA" clean -a -y 2>&1 | tail -3 | tee -a "$LOG"
"${PY[@]}" -m pip cache purge 2>&1 | tee -a "$LOG" || true
rm -rf "$ROOT/.cache/tmp" "$ROOT/.cache/pip"
mkdir -p "$ROOT/.cache/tmp"

# Measured once here; later stages add it to a walk of the rest of the
# tree instead of walking env/ again (see project_bytes in common.sh).
ENV_BYTES="$({ du -sb "$ENV" 2>/dev/null || true; } | awk '{print $1 + 0}')"
set_conf ENV_BYTES "$ENV_BYTES"
echo "01_install: env/ is $(gb "$ENV_BYTES") GB after cleaning" | tee -a "$LOG"

stage_finish "$STAGE" NA "masscube $("${PY[@]}" -c 'import importlib.metadata as m; print(m.version("masscube"))'), $R_VERSION, ThermoRawFileParser $TRFP_VERSION"
