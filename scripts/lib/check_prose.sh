#!/usr/bin/env bash
# Dashes, curly quotes, emoji and the banned words in tracked prose and comments; exit 1 on any hit.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
# grep -P won't run outside a UTF-8 locale, and its exit 2 must not pass as "none found"
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

# comment lines and docstrings; this file is left out since it holds the list
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

# only the statistical sense of robust is allowed, and grep can't tell, so these get read by eye
echo "== uses of 'robust', to be read (not failed) =="
{ grep -niw -- "robust" "${PROSE[@]}" 2>/dev/null || true; comments | grep -iw -- "robust" || true; } | sed 's/^/  /'

exit $status
