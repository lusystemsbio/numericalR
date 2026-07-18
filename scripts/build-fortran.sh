#!/bin/sh
#
# Compile the Fortran sources used by some chapters into shared libraries.
#
# Runs automatically as a Quarto pre-render step (see `pre-render` in
# _quarto.yml), so `quarto render` rebuilds the .so files for the current
# platform before the chapters that load them execute. The compiled .so files
# are build artifacts and are git-ignored; only the .f90 sources are tracked.
#
# Requirement: gfortran on PATH (macOS: `brew install gcc`; Debian/Ubuntu:
# `apt-get install gfortran`).
#
# POSIX sh only (macOS ships bash 3.2, and Quarto may invoke via /bin/sh).
set -eu

# Only active book chapters (NN-*/src). Skip archive/ and the legacy extra/
# materials, which are not part of the rendered book.
sources=$(find . -name '*.f90' -not -path './archive/*' -not -path './extra/*' | sort)

if [ -z "$sources" ]; then
  echo "build-fortran: no .f90 sources found; nothing to do."
  exit 0
fi

if ! command -v gfortran >/dev/null 2>&1; then
  echo "build-fortran: ERROR - gfortran not found on PATH, but Fortran sources exist:" >&2
  echo "$sources" | sed 's/^/  /' >&2
  echo "Install gfortran (macOS: 'brew install gcc'; Ubuntu: 'apt-get install gfortran')." >&2
  exit 1
fi

echo "$sources" | while IFS= read -r f90; do
  so="${f90%.f90}.so"
  echo "build-fortran: $f90 -> $so"
  gfortran -fpic -shared "$f90" -o "$so" || exit 1
done
