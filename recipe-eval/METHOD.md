# Recipe-box cold-generation evaluation — method

How we measure, and improve, whether each recipe box is a good enough spec that
an AI can reproduce the book's result from the recipe alone (with no book, no
web). The headline metric is the **Recipe Fidelity Score (RFS)**.

---

## 1. Goal and setup

For every recipe box, take the exact prompt its **Copy prompt** button produces,
hand it to a context-free AI, generate code, and check whether the code
reproduces the book's result. Fix the recipes that don't reproduce.

Settled parameters:
- **Language:** Python for the full sweep; R as a spot-check on a subset.
- **Generator model:** Claude Opus 4.8 (matches how the boxes were drafted; if a
  recipe fails even Opus cold, the recipe genuinely needs work).
- **Samples per box:** K = 5 (RFS is a rate, not a single draw).
- **Pilot first**, then the full 116-box sweep.

---

## 2. Generation — one-shot, context-free

Each sample: `claude -p "<prompt>" --model claude-opus-4-8 --allowedTools ""`
run in an **empty temp directory** — no book, no files, no web, and **one-shot**
(no self-execution / self-debugging). This measures the *prompt as written*, not
the model's ability to iterate to correctness. Fresh process per sample => no
leakage between samples or from the book. `generate.py` does this K times per box
and caches by file existence (so it is pausable/resumable and never regenerates).

---

## 3. Execution — a fair environment

Generated scripts run under `run.py` (parallel; each script is its own process).

**Fairness fix:** this machine has bleeding-edge numpy 2 / Python 3.14, which
*removed* APIs the AI reasonably uses (`np.trapz`, `ndarray.ptp()`,
`matplotlib.cm.get_cmap`, `boxplot(labels=)`). Counting those as failures would
measure library churn, not recipe quality. So scripts run in a dedicated
**Python-3.9 runner venv** with the common stable versions (numpy 1.26,
matplotlib 3.8, scipy, scikit-learn, networkx, umap-learn, pcurvepy2). This took
"can't run" from 27/578 down to 8/578 (~99% run).

An execution result is one of: **ran clean**, or **did not run** (real code bug,
failed compile, etc.).

---

## 4. Scoring — Recipe Fidelity Score (RFS)

**RFS(box) = (# of the box's K generations that reproduce the book's result) / K.**
Aggregate = mean RFS and the distribution of per-box RFS.

### 4.1 Compare the target quantity, wherever it lives
Do **not** compare screen output, and do **not** rely on an LLM's holistic
judgment (see the Haiku lesson below). Instead, for each box define a **target
quantity** — the key checkable object the recipe is about — and compare it to
ground truth. The target is a small object, e.g.:
- a scalar / final value (`N(100)=exp(10)`), the slope of a curve;
- a set of points (roots, 2 steady states, a tour, city order);
- a distribution (Gaussian samples, a degree distribution);
- a 2D field's morphology (Turing spots vs stripes, a spiral);
- an exact discrete structure (Boolean attractor set, an alignment).

The target lives in either the **plotted data** or the **printed output**; we
grab it from whichever holds it, so "no plot" is not a special case.

### 4.2 Where the data comes from
A matplotlib **capture shim** (`capture/sitecustomize.py`, auto-imported via
`PYTHONPATH`) records the actual arrays handed to `plot/scatter/hist/imshow/
contour/bar/...`, dumped to `$CAPTURE_OUT` as JSON. Run BOTH:
- the **book's Python reference** for the box -> ground-truth plotted data;
- each **generation** -> candidate plotted data.
For print-only targets, read the value from the generation's captured stdout.

### 4.3 Comparison rules (handle mismatched sampling carefully)
The book and the AI may save different numbers of points, so never compare arrays
element-by-element:
- **Curve `y(x)`:** interpolate BOTH series onto a common x-grid (their
  overlapping range) before computing relative / max error. (Confirmed real on
  2C.1: book plots the exact curve at 101 points, the AI at 1001.)
- **Parametric trajectory / orbit `(x(t),y(t))`:** resample by arc length onto a
  common parameterization (or a sampling-robust distance), not by index.
- **Point set (roots, steady states, tour):** order-independent set matching
  within tolerance; the count is part of the check.
- **Distribution:** Kolmogorov-Smirnov statistic (independent of sample count).
- **2D field:** pattern metric (peak/spot count, dominant wavelength via FFT),
  not a pixel diff.

### 4.4 Output type and the no-plot marks
Classify each box from its **Show** field (automatic): PLOT-output (78),
PRINT-only (14), or BOTH (24). Mark:
1. **Print-only boxes** — expect no figure; score the printed target, do not
   penalize the missing plot.
2. **Expected-plot-missing** — a plot-output box whose generation captured ZERO
   plotted series -> flag as an output failure (it did not produce the artifact
   the recipe asked for).

### 4.5 Method chosen (A + C + vision), with the alternative recorded
- **A (backbone): deterministic target-quantity comparison.** Per box: the key
  quantity + its source + tolerance; extract and compare to ground truth in code.
  No LLM for the numeric/curve/field/discrete majority. Objective, near-zero
  usage, handles no-plot by construction.
- **C (where cleaner): property checks.** Where a Verification is a property
  (energy conserved, exactly 3 roots, 2 steady states, attractor set exact,
  variance proportional to t), check the property on the generation's output
  instead of re-running the book.
- **Vision only for the ~5-8 pure-morphology boxes** (Turing spots/stripes,
  spiral) where there is no number — an LLM reads the figure.
- **Alternative B (fallback, not chosen):** an LLM as a structured *extractor*
  (reads code+stdout+figure, returns the named target values as JSON) with the
  pass/fail done deterministically in our code. Uniform for plot/print and robust
  (the LLM never judges), but costs some usage. Use only if per-box deterministic
  extractors get too fiddly.

Each verdict also carries a **failure class** — no-run / wrong-method /
wrong-params / wrong-output / numeric-off / underspecified-test — which drives
the prompt fixes.

---

## 5. Judge lessons (hard-won)
- **Ground the oracle in the book's ACTUAL output, not the verification prose.**
  On the pilot, 8D.2 scored 0/5 only because the criterion (SD ~= sqrt(x_bar))
  was stricter than the book's OWN figure, which shows SD well below sqrt(x_bar)
  at high copy number (finite tmax undersamples). The AI reproduced the book; the
  oracle was wrong. (Also flags a candidate accuracy fix to the book's wording.)
- **A weak LLM judge is unreliable.** A cheap holistic judge (Haiku) agreed with
  Opus on only 28% of pilot samples, and every disagreement was Haiku wrongly
  *failing* correct code (it failed all five 2C.1 runs that print N(100)=22026).
  This is why scoring moved to deterministic target-quantity comparison; any LLM
  is used only to *extract* or to read a genuinely visual figure, never to judge.
- **Keep the oracle in sync with recipe edits** (9B.3 stayed "fail" after its fix
  only because its criterion was not updated).
- **Some morphology targets are inherently fragile** (7E.3 spots): even with
  correct kinetics/params/run-time, cold generations often land in a labyrinth.
  Flag such boxes; don't force a false pass.

---

## 6. Optimization loop
- **Generation is one-shot** (measures the prompt).
- **Optimization is a diagnosis-driven loop, human-in-the-loop:** score -> read
  the failing generations + failure class -> edit the specific recipe field that
  caused it (Method / Test / Model / Show) -> regenerate that box -> re-score.
  Not the model rewriting its own prompt; edits are deliberate for quality.
- **Every edit is logged** with rationale and RFS before/after in
  `optimization_log.md`.
- **Baseline preserved:** all 116 pre-optimization boxes are archived in
  `original_boxes.{json,md}` for before/after diffs.

---

## 7. Files
- `prompts.py` — box -> exact Copy-prompt text (per language); `all_prompts.json`.
- `generate.py` — context-free generation (`$RECIPE_EVAL_BASE`, resumable).
- `run.py` — parallel execution (`$RECIPE_EVAL_PY` = the 3.9 runner venv).
- `capture/sitecustomize.py` — matplotlib plotted-data capture shim.
- `scorer.py` — verdict per generation (being moved from LLM-judge to the
  deterministic target-quantity design above).
- `original_boxes.{json,md}` — pre-optimization baseline.
- `optimization_log.md` — edit history + rationale.
- `gen/` — the generated scripts (cached).

---

## 8. Status
Generation complete (580/580). Execution ~99% run on the fair venv. Scoring is
being rebuilt from the failed LLM-judge to the deterministic target-quantity
method in section 4. Pilot (10 boxes, validated judge + manual): RFS 74% -> 92%
after two recipe fixes.
