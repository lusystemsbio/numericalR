# Proofreading tracker

Second-pass proofread of each migrated chapter against its original `.Rmd`.
Work top to bottom; flip `☐` → `☑` as each chapter is done, and jot issues in Notes.

## Workflow

Open a chapter next to its original in VS Code:

```sh
scripts/proofread.sh 2a          # opens 02a-*.qmd (edit) + archive/02A.Rmd (reference)
scripts/proofread.sh 2a --diff   # or the side-by-side diff view (old | new)
```

Then press **⌘\\** to put them in two panes. Edit the `.qmd`; keep the `.Rmd` as a
read-only reference. Prefer plain side-by-side for reading rewritten prose; use
`--diff` for a quick "what changed / was anything dropped" scan. The old file is
resolved from the `.qmd`'s `old_id:` frontmatter, so split chapters (2e and 2f both
map to `02E.Rmd`) open the right source.

After editing, re-render to check: `quarto render <file>.qmd --to html` (single
chapter), or a full `quarto render` before publishing.

## Per-chapter checklist (apply to each)

- [ ] **No dropped explanations** — the original's per-figure / interpretive sentences
      ("the plot below shows…", "for k=0.2…") survive, in a *shared* location (concept
      before the two impls, or discussion after Python), not buried in only one language.
- [ ] **R/Python parity** — the two implementations are the *same* algorithm and defaults;
      genuine language differences are noted inside their own block.
- [ ] **No dropped code** — nothing from the original `.Rmd` was silently lost (e.g. a
      validation step); diff-scan to be sure.
- [ ] **Citations** — every named method / dataset / figure cites its source (`[@key]`).
- [ ] **De-AI wording** — no em-dashes for asides, no "not X but Y", no rule-of-three
      triads, no "the takeaway is"; matches the author's voice.
- [ ] **Copyright** — flag any close paraphrase / verbatim external text; rewrite in
      original words.
- [ ] **Accuracy** — numbers/claims match what the code actually produces (cf. the
      2D.3 / 8D.2 fixes found during recipe evaluation).

## Chapters

| Done | Ch | Title | New file | Original | Notes |
|:--:|:--|:--|:--|:--|:--|
| ☐ | 01A | 1A. Getting started with R and Python | `01-intro/01a-getting-started.qmd` | `archive/01A.Rmd` |  |
| ☐ | 01B | 1B. Efficient programming | `01-intro/01b-efficient-programming.qmd` | `archive/01B.Rmd` |  |
| ☐ | 01C | 1C. Numerical methods | `01-intro/01c-numerical-methods.qmd` | `archive/01C.Rmd` |  |
| ☐ | 01D | 1D. Exercises | `01-intro/01d-exercises.qmd` | `archive/01D.Rmd` |  |
| ☐ | 02A | 2A. Modeling gene circuits with rate equations | `02-odes/02a-modeling-gene-circuits.qmd` | `archive/02A.Rmd` |  |
| ☐ | 02B | 2B. Numerical integration | `02-odes/02b-numerical-integration.qmd` | `archive/02B.Rmd` |  |
| ☐ | 02C | 2C. Modeling bacterial growth | `02-odes/02c-bacterial-growth.qmd` | `archive/02C.Rmd` |  |
| ☐ | 02D | 2D. Stability and effective potential | `02-odes/02d-effective-potential.qmd` | `archive/02D.Rmd` |  |
| ☐ | 02E | 2E. Bifurcation: finding roots/steady states | `02-odes/02e-bifurcation.qmd` | `archive/02E.Rmd` | split from 02E.Rmd (with 02f) |
| ☐ | 02F | 2F. Bifurcation: finding curves | `02-odes/02f-bifurcation-curves.qmd` | `archive/02E.Rmd` | split from 02E.Rmd (with 02e) |
| ☐ | 02G | 2G. Exercises | `02-odes/02g-exercises.qmd` | `archive/02F.Rmd` |  |
| ☐ | 03A | 3A. Nullclines | `03-phase-plane/03a-nullclines.qmd` | `archive/03A.Rmd` |  |
| ☐ | 03B | 3B. Steady states in 2D | `03-phase-plane/03b-steady-states-2d.qmd` | `archive/03B.Rmd` |  |
| ☐ | 03C | 3C. Modeling a chemostat | `03-phase-plane/03c-chemostat.qmd` | `archive/03C.Rmd` |  |
| ☐ | 03D | 3D. Modeling predator-prey dynamics | `03-phase-plane/03d-predator-prey.qmd` | `archive/03D.Rmd` |  |
| ☐ | 03E | 3E. Bifurcation for two-variable systems | `03-phase-plane/03e-bifurcation-2d.qmd` | `archive/03E.Rmd` |  |
| ☐ | 03F | 3F. Separatrix | `03-phase-plane/03f-separatrix.qmd` | `archive/03F.Rmd` |  |
| ☐ | 03G | 3G. Effective potential revisited | `03-phase-plane/03g-effective-potential-2d.qmd` | `archive/03G.Rmd` |  |
| ☐ | 03H | 3H. Multi-component systems | `03-phase-plane/03h-multi-component.qmd` | `archive/03H.Rmd` |  |
| ☐ | 03I | 3I. Exercises | `03-phase-plane/03i-exercises.qmd` | `archive/03I.Rmd` |  |
| ☐ | 04A | 4A. Delayed differential equations | `04-time-delays/04a-delayed-differential-equations.qmd` | `archive/04A.Rmd` |  |
| ☐ | 04B | 4B. Modeling systems with time delays | `04-time-delays/04b-modeling-time-delays.qmd` | `archive/04B.Rmd` |  |
| ☐ | 04C | 4C. Delays from indirect interactions | `04-time-delays/04c-indirect-interactions.qmd` | `archive/04C.Rmd` | R prose existed; Python written new |
| ☐ | 04D | 4D. Exercises | `04-time-delays/04d-exercises.qmd` | `archive/04D.Rmd` |  |
| ☐ | 05A | 5A. Integrators for second-order ODEs | `05-molecular-dynamics/05a-second-order-integrators.qmd` | `archive/05A.Rmd` |  |
| ☐ | 05B | 5B. Orbital motion | `05-molecular-dynamics/05b-orbital-motion.qmd` | `archive/05B.Rmd` |  |
| ☐ | 05C | 5C. Modeling a box of particles | `05-molecular-dynamics/05c-box-of-particles.qmd` | `archive/05C.Rmd` |  |
| ☐ | 05D | 5D. Exercises | `05-molecular-dynamics/05d-exercises.qmd` | `archive/05D.Rmd` |  |
| ☐ | 06A | 6A. Random number generators | `06-stochastic/06a-random-number-generators.qmd` | `archive/06A.Rmd` |  |
| ☐ | 06B | 6B. Brownian motion | `06-stochastic/06b-brownian-motion.qmd` | `archive/06B.Rmd` |  |
| ☐ | 06C | 6C. SDE integrators | `06-stochastic/06c-sde-integrators.qmd` | `archive/06C.Rmd` |  |
| ☐ | 06D | 6D. Stochastic state transitions | `06-stochastic/06d-stochastic-transitions.qmd` | `archive/06D.Rmd` |  |
| ☐ | 06E | 6E. Exercises | `06-stochastic/06e-exercises.qmd` | `archive/06E.Rmd` |  |
| ☐ | 07A | 7A. Modeling diffusion | `07-pde/07a-modeling-diffusion.qmd` | `archive/07A.Rmd` |  |
| ☐ | 07B | 7B. Reaction-diffusion systems | `07-pde/07b-reaction-diffusion.qmd` | `archive/07B.Rmd` |  |
| ☐ | 07C | 7C. Turing instability | `07-pde/07c-turing-instability.qmd` | `archive/07C.Rmd` |  |
| ☐ | 07D | 7D. Pattern formation in Dictyostelium | `07-pde/07d-pattern-formation-dictyostelium.qmd` | `archive/07D.Rmd` |  |
| ☐ | 07E | 7E. Two-dimensional reaction-diffusion | `07-pde/07e-2d-reaction-diffusion.qmd` | `archive/07E.Rmd` | NEW chapter; no old source (07E.Rmd = the old exercises → 07f) |
| ☐ | 07F | 7F. Exercises | `07-pde/07f-exercises.qmd` | `archive/07E.Rmd` | = old 07E exercises |
| ☐ | 08A | 8A. Monte Carlo methods | `08-monte-carlo/08a-monte-carlo-methods.qmd` | `archive/08A.Rmd` |  |
| ☐ | 08B | 8B. Metropolis algorithm | `08-monte-carlo/08b-metropolis-algorithm.qmd` | `archive/08B.Rmd` |  |
| ☐ | 08C | 8C. Particles in a box: MCMC sampling | `08-monte-carlo/08c-particles-in-a-box.qmd` | `archive/08C.Rmd` |  |
| ☐ | 08D | 8D. Gillespie algorithm | `08-monte-carlo/08d-gillespie-algorithm.qmd` | `archive/08D.Rmd` |  |
| ☐ | 08E | 8E. Exercises | `08-monte-carlo/08e-exercises.qmd` | `archive/08E.Rmd` | new GRN problem 8E.2 added |
| ☐ | 09A | 9A. MCMC optimization methods | `09-optimization/09a-mcmc-optimization.qmd` | `archive/09A.Rmd` |  |
| ☐ | 09B | 9B. Dynamic programming | `09-optimization/09b-dynamic-programming.qmd` | `archive/09B.Rmd` |  |
| ☐ | 09C | 9C. Genetic algorithm | `09-optimization/09c-genetic-algorithm.qmd` | `archive/09C.Rmd` |  |
| ☐ | 09D | 9D. Exercises | `09-optimization/09d-exercises.qmd` | `archive/09D.Rmd` |  |
| ☐ | 10A | 10A. Dimensionality reduction | `10-high-dim-data/10a-dimensionality-reduction.qmd` | `archive/10A.Rmd` |  |
| ☐ | 10B | 10B. Clustering | `10-high-dim-data/10b-clustering.qmd` | `archive/10B.Rmd` |  |
| ☐ | 10C | 10C. Network algorithms | `10-high-dim-data/10c-network-algorithms.qmd` | `archive/10C.Rmd` | heavily authored from PPTs |
| ☐ | 10D | 10D. Exercises | `10-high-dim-data/10d-exercises.qmd` | `archive/10D.Rmd` |  |
