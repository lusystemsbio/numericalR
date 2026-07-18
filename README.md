# Numerical Quartorial

### *Systems Biology in R & Python* — abbreviated **NumericalQ**

#### [Mingyang Lu](https://lusystemsbio.northeastern.edu) — Lu Lab for Computational Systems Biology, Northeastern University

This book is a hands-on introduction to the numerical methods and algorithms used
in computational systems biology: modeling gene-regulatory circuits, solving
ordinary and stochastic differential equations, running simulations, performing
optimization, and analyzing high-dimensional data. Rather than treating solvers
as black boxes, each topic develops the underlying theory and then implements the
algorithm from scratch, in the context of real problems in bioengineering,
biomedical engineering, and data science.

Every method is presented **language-agnostically** and then implemented in both
**R** and **Python**, so the same concept can be studied and run in either
language. The book is built with [Quarto](https://quarto.org) and rendered to
HTML, PDF, and EPUB from a single source.

The live HTML version is published at
<https://lusystemsbio.github.io/numericalR>.

## Repository layout

- `NN-<topic>/` — chapter folders grouped by part (e.g. `02-odes/`), each holding
  `.qmd` chapters, their figures under `images/`, and any compiled-code sources
  under `src/`.
- `index.qmd` — book preface; `references.qmd` — consolidated bibliography.
- `references.bib` — shared bibliography for the whole book.
- `scripts/` — build helpers run automatically during rendering.
- `archive/` — original R Markdown / Jupyter source material retained for
  reference during the migration; not part of the rendered book.

## Building the book

### Prerequisites

- **[Quarto](https://quarto.org/docs/get-started/)** (1.9 or newer).
- **R** (4.6+) with the packages: `reticulate`, `deSolve`, `microbenchmark`,
  `ggplot2` (install via `install.packages(...)`).
- **Python** (3.x) on your `PATH` with `numpy`, `scipy`, `matplotlib`, `jupyter`,
  and `ipykernel`. R and Python chunks run in the same document via `knitr` +
  `reticulate`; the interpreter is pinned once in the repo-root `.Rprofile`.
- **gfortran** — required for the pre-build step below
  (macOS: `brew install gcc`; Debian/Ubuntu: `apt-get install gfortran`).
- **TinyTeX** for PDF output: `quarto install tinytex` (one time).

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
quarto render              # build HTML, PDF, and EPUB into _book/
quarto render --to html    # a single format
quarto preview             # live-reloading local preview
```

The rendered site is written to `_book/` (git-ignored).
