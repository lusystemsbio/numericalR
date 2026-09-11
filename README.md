<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/lockup-dark.png">
  <img src="assets/lockup-light.png" alt="Numerical ReciPy — Computational Systems Biology and Biophysics in R and Python" width="620">
</picture>

#### [Mingyang Lu](https://lusystemsbio.northeastern.edu) — Lu Lab for Computational Systems Biology, Northeastern University

This book is a hands-on introduction to the numerical methods and algorithms used
in computational systems biology and biophysics: modeling gene-regulatory circuits,
solving ordinary, stochastic, and partial differential equations, running molecular
dynamics and Monte Carlo simulations, performing global optimization, and analyzing
high-dimensional data. Rather than treating solvers as black boxes, each topic
develops the underlying theory and then implements the algorithm from scratch, in
the context of real problems in bioengineering, biomedical engineering, and data
science.

It is written for graduate and advanced undergraduate students, and for researchers
who want to build and run their own models. The mathematics assumes calculus and a
first course in ordinary differential equations. The programming assumes some prior
experience, though not necessarily in R or Python, since Part 1 sets up both.

Every method is presented independently of language and then implemented in both
**R** and **Python**, so the same concept can be studied and run in either language.
The book is built with [Quarto](https://quarto.org) and rendered to HTML, PDF, and
EPUB from a single source. Quarto executes both languages in the same document, so
every result and figure is produced by the code printed beside it.

Each major method is also stated as a **recipe**: a short box giving the objective,
the model, the numerical method, a concrete test case, the output to produce, and a
verification check, written fully enough that the method can be rebuilt without
reading the implementation. A recipe doubles as a structured prompt for an AI
assistant, and the HTML edition puts a **Copy prompt** button on every box.

The live HTML version is published at
<https://lusystemsbio.github.io/numericalReciPy/>.

## Repository layout

- `NN-<topic>/` — chapter folders grouped by part (e.g. `02-odes/`), each holding
  `.qmd` chapters, their figures under `images/`, and any compiled-code sources
  under `src/`.
- `index.qmd` — book preface; `references.qmd` — consolidated bibliography.
- `_quarto.yml` — book configuration (HTML); `_quarto-print.yml` — the PDF and
  EPUB formats, applied by the `print` profile.
- `references.bib` — shared bibliography for the whole book.
- `requirements.txt` / `scripts/install-packages.R` — the full Python/R package
  lists needed to render the book (see Prerequisites below).
- `scripts/` — build helpers (some run automatically during rendering) and
  authoring utilities.
- `archive/` — the pre-Quarto `.Rmd`/`.html` originals, kept as a reference for the
  proofreading pass; not part of the rendered book. See `proofread.md` for how to
  compare a chapter against its original.
- `appendix-*.qmd` — method reference, R/Python library map, and the recipe index.
- `archive/` — the original R Markdown and Jupyter source material, kept for
  reference; not part of the rendered book.

## Building the book

### Prerequisites

- **[Quarto](https://quarto.org/docs/get-started/)** (1.9 or newer).
- **R** (4.6+) with the packages listed in `scripts/install-packages.R`; install
  them all with:
  ```bash
  Rscript scripts/install-packages.R
  ```
- **Python** (3.x) on your `PATH` with the packages listed in `requirements.txt`;
  install them all with:
  ```bash
  python3 -m pip install -r requirements.txt
  ```
  R and Python chunks run in the same document via `knitr` + `reticulate`; the
  repo-root `.Rprofile` picks the first interpreter that has these packages, so a
  Homebrew or python.org install both work. To force a specific one, set
  `RETICULATE_PYTHON=/path/to/python3`.
- **macOS + python.org Python only:** run the installer's certificate step once.
  ```bash
  "/Applications/Python 3.14/Install Certificates.command"    # match your version
  ```
  The python.org build does not point Python at a CA bundle, so without this every
  HTTPS download fails. Chapter 1D.1 fetches a structure from the Protein Data Bank
  and would fail with `'NoneType' object has no attribute 'readlines'`, which is
  Biopython returning `None` for the failed download rather than a certificate error.
  Homebrew Python does not need this.
- **gfortran** — required for the pre-build step below
  (macOS: `brew install gcc`; Debian/Ubuntu: `apt-get install gfortran`).
  On macOS, Homebrew's gfortran often cannot find the Command Line Tools SDK and
  fails to link with `ld: library 'System' not found`. Point it at the SDK:
  ```bash
  export SDKROOT=$(xcrun --show-sdk-path)     # add to ~/.zshrc to make it stick
  ```
  Also make sure Homebrew's `bin` is on your `PATH` (`/opt/homebrew/bin` on Apple
  Silicon), or `gfortran` will not be found at all.
- **ImageMagick** (`convert`/`magick` on your `PATH`) — required for the
  figure-crop hook (macOS: `brew install imagemagick`; Debian/Ubuntu:
  `apt-get install imagemagick`).
- **TinyTeX** — only needed for the PDF edition: `quarto install tinytex` (one time).

If a chapter fails to render with a missing-package error, re-run the two
install commands above before troubleshooting further — package lists here
have drifted out of sync with the actual environment before.

### Platform support

Developed and rendered on macOS; Linux is expected to work unchanged. Windows is
not currently supported for two parts of the book: the chapters that call compiled
Fortran (1B, 2B, 6D) load `.so` libraries built by a POSIX shell script, where
Windows R needs `.dll`, and the parallel examples in 1B rely on Unix `fork`, which
Windows lacks for both `multiprocessing` and `parallel::mclapply`. Everything else
is platform-neutral. On Windows, WSL2 avoids both problems.

### Pre-build step (automatic)

Some chapters call compiled Fortran routines. The shared libraries are **not**
committed; they are compiled from the `.f90` sources for your platform by
`scripts/build-fortran.sh`, which Quarto runs automatically as a
[`pre-render`](https://quarto.org/docs/projects/scripts.html#pre-and-post-render)
step (configured in `_quarto.yml`). You do not need to run it by hand — it fires
on every `quarto render`. It only requires `gfortran` to be installed; if you
want to run it standalone:

```bash
sh scripts/build-fortran.sh
```

### Render

```bash
quarto render              # build the HTML book into _book/
quarto preview             # live-reloading local preview
```

The rendered site is written to `_book/` (git-ignored).

To work on a single chapter, pass its path:

```bash
quarto render 02-odes/02b-numerical-integration.qmd   # -> _book/02-odes/...html
quarto preview 02-odes/02b-numerical-integration.qmd  # same, live-reloading
```

A single-chapter render always re-executes that chapter, so it is the way to pick up
a package or data change. A whole-book `quarto render` is incremental instead:
`execute: freeze: auto` re-runs only the chapters whose source changed and loads the
rest from the `_freeze/` cache, so a full build after one edit takes seconds. Keep
`_freeze/`; to force a chapter to re-run, render that chapter on its own (or delete
its `_freeze/` entry). Either way the pre-render steps still run first.

PDF and EPUB are **not** built by default, since they are slow and rarely needed
while writing. Those two formats are configured in `_quarto-print.yml` and are
enabled by Quarto's `print` profile:

```bash
scripts/build-print.sh          # both
scripts/build-print.sh pdf      # just one
quarto render --profile print --to epub    # or call Quarto directly
```

Both are whole-book merged documents, so always render the whole project rather
than a single chapter. The PDF additionally needs TinyTeX (see Prerequisites).

## Single-language practice copies

Every chapter presents both an R and a Python implementation. For readers who want to
work in just one language, a setup script generates a single-language copy of each
dual-language chapter, written **in place** beside the original:

```bash
python3 scripts/make-lang-variants.py
```

For `02-odes/02b-numerical-integration.qmd` this writes
`02-odes/02b-numerical-integration-R.qmd` (Python implementation removed, `engine: knitr`)
and `02-odes/02b-numerical-integration-py.qmd` (R removed, `engine: jupyter`, so it needs
no R). Because each copy lives in its chapter's own folder, all of the chapter's
resources (`images/`, `src/`, data files) resolve normally, and you can open and run the
copy directly. These copies are git-ignored and are **not** part of the rendered book;
regenerate them any time by re-running the script. The online HTML edition additionally
offers a **Both / R / Python** switch at the top of each chapter to hide one language's
code while reading.

## License

The source code in this repository (the R, Python, and Fortran in the chapters,
and the build scripts) is licensed under the MIT license (`LICENSE-CODE`) and
can be reused freely.

The book's text, figures, and exercises are not covered by that license.
Copyright (c) 2021-2026 Mingyang Lu, Lu Lab for Computational Systems Biology.
All rights reserved. To use them in a course or elsewhere, please contact the
author.
