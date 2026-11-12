#!/usr/bin/env bash
# The checks run before the repository is handed over. Each prints what it
# found; the exit status is 1 if any check fails.
#
# Usage: scripts/lib/final_checks.sh
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT" || exit 1
# shellcheck source=common.sh
source "$ROOT/scripts/lib/common.sh"
load_conf
status=0

echo "== 1. shellcheck on every shell script =="
sc=""
for c in "$ROOT/env/bin/shellcheck.exe" "$ROOT/env/Library/bin/shellcheck.exe" "$ROOT/env/bin/shellcheck" "$(command -v shellcheck || true)"; do
    [[ -n "$c" && -x "$c" ]] && { sc="$c"; break; }
done
if [[ -n "$sc" ]]; then
    if "$sc" -x -e SC1091 scripts/*.sh scripts/lib/*.sh run_all.sh; then echo "clean"; else status=1; fi
else
    echo "shellcheck not found"; status=1
fi

echo "== 2. every figure the README shows exists and a script writes it =="
while read -r fig; do
    if [[ ! -s "$fig" ]]; then echo "MISSING: $fig"; status=1; continue; fi
    stem="$(basename "$fig" .png)"
    # Mirror plots are named per feature at run time (ms2_<mode>_<id>).
    [[ "$stem" == ms2_* ]] && stem="ms2_%s_%d"
    grep -qF -- "$stem" scripts/*.R || { echo "NOT WRITTEN BY A SCRIPT: $fig"; status=1; }
done < <(grep -oE 'figures/[A-Za-z0-9_./-]+\.png' README.md | sort -u)
echo "checked $(grep -oE 'figures/[A-Za-z0-9_./-]+\.png' README.md | sort -u | wc -l) figures"

echo "== 3. prose: dashes, curly quotes, emoji, banned words =="
if bash scripts/lib/check_prose.sh > "$ROOT/logs/prose_check.log" 2>&1; then echo "clean"; else cat "$ROOT/logs/prose_check.log"; status=1; fi

echo "== 4. no tracked file over 50 MB, repository under 200 MB =="
big=$(git ls-files -z | xargs -0 stat -c '%s %n' 2>/dev/null | awk '$1 > 52428800')
if [[ -n "$big" ]]; then echo "$big"; status=1; else echo "largest tracked file: $(git ls-files -z | xargs -0 stat -c '%s %n' | sort -n | tail -1)"; fi
pack=$(git count-objects -v | awk '/size-pack/ {print $2}')
echo "pack size ${pack} KiB"
(( pack < 200 * 1024 )) || status=1

echo "== 5. every stage has its stamp and resource row =="
for st in 01_install 02_download_metadata 03_sample_table 04_pilot_pos 05_rules_pos 04_process_pos 06_align_pos \
          04_pilot_neg 05_rules_neg 04_process_neg 06_align_neg 07_annotate 08_qc_filters 09_statistics \
          10_pathways 10_classifier 10_permanova 11_replication 12_figures; do
    if [[ -f "logs/$st.done" && -f "logs/$st.resources.tsv" ]]; then
        printf "%-22s %s\n" "$st" "$(tail -1 "logs/$st.resources.tsv" | cut -f2-5)"
    else
        echo "$st: no stamp or no resource row"; status=1
    fi
done

echo "== 6. disk: the literal du -sb of the project against the budget =="
bytes="$(project_bytes_du)"
echo "du -sb: $(gb "$bytes") GB against ${DISK_GB} GB"
(( bytes < DISK_GB * 1024 * 1024 * 1024 )) || status=1

echo "== 7. authorship and trace check of the pre-push hook, on HEAD =="
if [[ -x .git/hooks/pre-push ]]; then
    printf 'refs/heads/main %s refs/heads/main %s\n' "$(git rev-parse HEAD)" 0000000000000000000000000000000000000000 \
        | bash .git/hooks/pre-push origin manual && echo "passes" || status=1
else
    echo "no pre-push hook installed in this clone"
fi

exit $status
