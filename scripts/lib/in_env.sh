#!/usr/bin/env bash
# Run a program from the project's conda environment without activating
# it. A conda Python or R on Windows finds its DLLs only when the
# environment's Library directories are on PATH, which activation would
# normally arrange.
#
# Usage: scripts/lib/in_env.sh <python|Rscript|program> [args...]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENV="$ROOT/env"
if [[ ! -d "$ENV" ]]; then
    echo "in_env.sh: no environment at $ENV; run scripts/01_install.sh" >&2
    exit 1
fi
export PATH="$ENV:$ENV/Library/mingw-w64/bin:$ENV/Library/usr/bin:$ENV/Library/bin:$ENV/Scripts:$ENV/bin:$ENV/Lib/R/bin/x64:$PATH"
export CONDA_PREFIX="$ENV"
# Matplotlib and numba write caches under the home directory unless told
# otherwise; inside the project they are counted and cleaned with the rest.
export MPLCONFIGDIR="$ROOT/.cache/matplotlib"
export NUMBA_CACHE_DIR="$ROOT/.cache/numba"
export PYTHONDONTWRITEBYTECODE=1
prog="$1"; shift
case "$prog" in
    python)
        for p in "$ENV/python.exe" "$ENV/bin/python"; do [[ -x "$p" ]] && exec "$p" "$@"; done ;;
    Rscript)
        for p in "$ENV/Scripts/Rscript.exe" "$ENV/Lib/R/bin/Rscript.exe" "$ENV/bin/Rscript"; do [[ -x "$p" ]] && exec "$p" "$@"; done ;;
    *)
        exec "$prog" "$@" ;;
esac
echo "in_env.sh: $prog not found in $ENV" >&2
exit 1
