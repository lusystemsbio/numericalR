#!/usr/bin/env bash
# Open a chapter and its original .Rmd side by side in VS Code for proofreading.
#
#   scripts/proofread.sh 2a          # new 02a-*.qmd (edit) + archive/02A.Rmd (reference)
#   scripts/proofread.sh 2a --diff   # open them in the side-by-side diff view instead
#
# The old file is found from the .qmd's own `old_id:` frontmatter, so split
# chapters (2e and 2f both map to 02E.Rmd) resolve correctly.
set -euo pipefail
cd "$(dirname "$0")/.."

code_id="${1:-}"; mode="${2:-}"
[ -z "$code_id" ] && { echo "usage: $0 <chapter e.g. 2a> [--diff]"; exit 1; }

part="$(printf '%02d' "${code_id%%[a-z]*}")"      # 2a -> 02
letter="${code_id##*[0-9]}"                       # 2a -> a
qmd="$(ls "${part}"-*/"${part}${letter}"-*.qmd 2>/dev/null | grep -vE -- '-(R|py)\.qmd$' | head -1 || true)"
[ -z "$qmd" ] && { echo "no .qmd found for '$code_id'"; exit 1; }

old_id="$(grep -m1 '^old_id:' "$qmd" | sed 's/old_id://; s/[" ]//g')"
rmd="archive/${old_id}.Rmd"

echo "new: $qmd"
echo "old: $rmd  (old_id=$old_id)"
if [ ! -f "$rmd" ]; then echo "  (no archived .Rmd; this chapter is newly authored)"; code "$qmd"; exit 0; fi

if [ "$mode" = "--diff" ]; then
  code --diff "$rmd" "$qmd"                        # side-by-side diff (old | new)
else
  code "$qmd" "$rmd"                               # two tabs; press Cmd-\ to split into panes
  echo "  tip: press Cmd-\\ to put them side by side (edit the .qmd, keep the .Rmd as reference)"
fi
