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
