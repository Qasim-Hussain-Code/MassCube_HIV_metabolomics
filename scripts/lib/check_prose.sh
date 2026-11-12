#!/usr/bin/env bash
# Search the tracked prose and code comments for the characters and words
# the writing rules exclude. Exit 1 if anything is found.
#
# Characters (any tracked text file): en dash U+2013, em dash U+2014,
# curly quotes U+2018, U+2019, U+201C, U+201D, and emoji and pictographs
# (U+1F000 and above, U+2600 to U+27BF).
# Words and phrases (README.md, preregistration.md, data/README.md, and
# every comment line in the code): the list below. "robust" is allowed
# only in its statistical sense, which a pattern cannot judge, so its
# uses are printed for reading rather than failed.
#
# Usage: scripts/lib/check_prose.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
# grep -P refuses to run outside a UTF-8 locale, and a refusal (exit 2)
# must not read as "nothing found".
export LC_ALL=C.UTF-8 LANG=C.UTF-8

mapfile -t TEXT < <(git ls-files | grep -vE '\.(png|pdf|gz|jpg|jpeg|tiff|pkl|zip)$')
PROSE=(README.md preregistration.md data/README.md)
status=0

echo "== dashes, curly quotes, emoji =="
rc=0
grep -nP '[\x{2013}\x{2014}\x{2018}\x{2019}\x{201C}\x{201D}\x{2600}-\x{27BF}\x{1F000}-\x{1FFFF}]' "${TEXT[@]}" || rc=$?
case $rc in
    0) status=1 ;;
    1) echo "none" ;;
    *) echo "the character check itself failed (grep exit $rc)"; status=1 ;;
esac

# Comment lines only: '#' in shell, Python and R, plus Python docstrings,
# which are written as prose. This file is excluded because it holds the list.
comments() {
    local f
    for f in $(git ls-files '*.sh' '*.py' '*.R' | grep -v 'check_prose.sh'); do
        case "$f" in
            *.py) awk -v f="$f" '/"""/ {n = gsub(/"""/, "&"); if (n == 1) inside = !inside; print f ":" NR ":" $0; next}
                                 inside || /^[[:space:]]*#/ || /[[:space:]]#[[:space:]]/ {print f ":" NR ":" $0}' "$f" ;;
            *) awk -v f="$f" '/^[[:space:]]*#/ || /[[:space:]]#[[:space:]]/ {print f ":" NR ":" $0}' "$f" ;;
        esac
    done
}

WORDS='delve|tapestry|testament|crucial|intricate|meticulous|seamless|leverag|showcase|foster|transformative|groundbreaking|comprehensive|holistic|multifaceted|realm|underscore|pivotal|plays a vital role|not just|not only|in conclusion|moreover|furthermore|additionally|highlight'
echo "== banned words and phrases =="
hits="$( { grep -niE -- "$WORDS" "${PROSE[@]}" 2>/dev/null || true; comments | grep -iE -- "$WORDS" || true; } )"
if [[ -n "$hits" ]]; then echo "$hits"; status=1; else echo "none"; fi

echo "== uses of 'robust', to be read (not failed) =="
{ grep -niw -- "robust" "${PROSE[@]}" 2>/dev/null || true; comments | grep -iw -- "robust" || true; } | sed 's/^/  /'

exit $status
