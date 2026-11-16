#!/usr/bin/env bash
# in_env.sh python|Rscript|prog args... -- runs from env/ without conda activate
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENV="$ROOT/env"
if [[ ! -d "$ENV" ]]; then
    echo "in_env.sh: no environment at $ENV; run scripts/01_install.sh" >&2
    exit 1
fi
# conda Python and R on Windows only find their DLLs with Library/ on PATH
export PATH="$ENV:$ENV/Library/mingw-w64/bin:$ENV/Library/usr/bin:$ENV/Library/bin:$ENV/Scripts:$ENV/bin:$ENV/Lib/R/bin/x64:$PATH"
export CONDA_PREFIX="$ENV"
# keep matplotlib and numba caches out of ~ so the disk budget counts them
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
