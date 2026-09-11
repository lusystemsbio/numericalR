#!/usr/bin/env bash
# Open a chapter and its original .Rmd side by side in VS Code for proofreading.
#
#   scripts/proofread.sh 2a          # new 02a-*.qmd (edit) + archive/02A.Rmd (reference)
#   scripts/proofread.sh 2a --diff   # open them in the side-by-side diff view instead
#
# The old file is found from the .qmd's own `old_id:` frontmatter, so split
# chapters (2e and 2f both map to 02E.Rmd) resolve correctly.
#
# archive/ has since been deleted from the working tree, so the reference .Rmd is
# recovered from git history into a temp file. Nothing to restore by hand.
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

# archive/ is no longer in the working tree; pull the reference out of git history.
if [ ! -f "$rmd" ]; then
  last="$(git rev-list -1 HEAD -- "$rmd" 2>/dev/null || true)"
  if [ -n "$last" ]; then
    tmp="$(mktemp -t "proofread-${old_id}")-${old_id}.Rmd"
    # the newest version of the path, else its state just before the commit that removed it
    if git show "${last}:${rmd}" > "$tmp" 2>/dev/null \
       || git show "${last}^:${rmd}" > "$tmp" 2>/dev/null; then
      rmd="$tmp"
    else
      rm -f "$tmp"
    fi
  fi
fi

echo "new: $qmd"
echo "old: $rmd  (old_id=$old_id)"
if [ ! -f "$rmd" ]; then echo "  (no archived .Rmd; this chapter is newly authored)"; code "$qmd"; exit 0; fi

if [ "$mode" = "--diff" ]; then
  code --diff "$rmd" "$qmd"                        # side-by-side diff (old | new)
else
  code "$qmd" "$rmd"                               # two tabs; press Cmd-\ to split into panes
  echo "  tip: press Cmd-\\ to put them side by side (edit the .qmd, keep the .Rmd as reference)"
fi
