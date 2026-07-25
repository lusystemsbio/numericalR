# Recipe prompt-optimization history

Each entry records one recipe-box edit made to improve its cold-generation
reproduction rate (RFS = fraction of K=5 context-free generations that reproduce
the book's result). Newest first. RFS values are from the validated judge
(pilot: Opus-judged + manually checked).

Format: box | RFS before -> after | field changed | rationale (before -> after).

---

## 7E.3.1  (07-pde/07e-2d-reaction-diffusion.qmd) — commit f988b19
- **RFS:** 0/5 -> 1/5  (still fragile; see note)
- **Field:** Test
- **Failure class:** wrong-output (pattern morphology)
- **Before:** "... nearly uniform initial u = v = 1 with +-0.1 noise (seed 10), identical to 7E.2."
- **After:** "... noise (seed 10); integrate to `t = 100` in successive blocks of 20 (the same initial condition and run length as 7E.2)."
- **Rationale:** The Test never stated the total run time. Cold generations stopped mid-transient and produced a labyrinth instead of the book's settled array of spots (the book reaches spots only by t~=60-100). Adding the run length fixed the omission.
- **Note:** Only 1/5 after the fix. Even with correct kinetics/params/run time, the activator-inhibitor spot morphology is implementation/realization-sensitive and does not reliably reproduce from the recipe alone. Flagged as inherently fragile; a deeper fix (or accepting it as hard) is TBD.

## 9B.3.1  (09-optimization/09b-dynamic-programming.qmd) — commit f988b19
- **RFS:** 2/5 -> 5/5
- **Field:** Test
- **Failure class:** underspecified-test
- **Before:** "10 cities at fixed coordinates, the tour starting and ending at city 1."
- **After:** "10 cities with x = (0, -28.87, ...) and y = (0, 0, 43.39, ...), Euclidean distances; the tour starting and ending at city 1."
- **Rationale:** The Test said "fixed coordinates" but listed none, so a context-free AI invented its own cities and could not reproduce the book's specific tour (optimal length 193.73). Adding the actual coordinates made it reproducible.

## 6D.3.1 (06-stochastic/06d-stochastic-transitions.qmd) — RECIPE BOX REMOVED
- **RFS:** 2/5 -> (removed)
- **Rationale:** The box asked a cold AI to write + compile a Fortran MFPT routine and bind it from R and Python, then run a long (t=1e5) simulation. Two generations came back empty and one failed to compile. This is inherent difficulty (compiled code + long simulation), not a promptable gap, so per author's decision the recipe box is removed (the 6D.3 section keeps its content). It is not a good candidate for a self-contained recipe.

## 2D.3.1 (02-odes/02d-effective-potential.qmd) — PROSE FIX (not a generation failure)
- **Field:** Verification
- **Before:** "two basins (near 100 and 300 nM) split by a barrier near 200 nM"
- **After:**  "two basins (near 72 and 332 nM) split by a barrier near 171 nM"
- **Rationale:** The generations were correct; the box's numbers were wrong. This self-activating gene at k=0.15 has stable states ~71.5/331.6 and an unstable state ~170.8 (the same roots as 2E.3), not 100/300/200. Surfaced because the LLM judge failed correct generations against the inaccurate prose.

## 8D.2.1 (08-monte-carlo/08d-gillespie-algorithm.qmd) — PROSE FIX
- **Field:** Verification
- **Before:** "The standard deviation follows sqrt(x_bar) (Poisson)"
- **After:**  "grows with the mean ... below sqrt(x_bar) at this finite run length (a single trajectory undersamples) ... relative noise ~1/sqrt(x_bar)"
- **Rationale:** The book's OWN figure shows the measured SD sitting well below sqrt(x_bar) at high copy number (finite tmax undersampling); "follows sqrt(x_bar)" overstated it.
