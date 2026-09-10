#!/bin/sh
# Build the print editions (PDF and EPUB) of the book.
#
# A plain `quarto render` builds HTML only, because the pdf/epub formats are
# not declared in `_quarto.yml`; they live in `_quarto-print.yml` and are
# applied by Quarto's `print` profile. Run this script when you want them.
#
# Equivalent to running, from the repo root:
#     quarto render --profile print --to pdf
#     quarto render --profile print --to epub
#
# The PDF needs TinyTeX (`quarto install tinytex`, one time). Both are whole-book
# merged documents, so always render the whole project, never a single chapter.
set -eu
cd "$(dirname "$0")/.."

target="${1:-all}"

case "$target" in
  pdf|epub)
    quarto render --profile print --to "$target"
    ;;
  all)
    quarto render --profile print --to pdf
    quarto render --profile print --to epub
    ;;
  *)
    echo "usage: $0 [pdf|epub|all]" >&2
    exit 1
    ;;
esac

echo
echo "Output in _book/:"
ls -1 _book/*.pdf _book/*.epub 2>/dev/null || echo "  (nothing found)"
