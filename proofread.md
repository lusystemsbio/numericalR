# Proofreading a migrated chapter

`scripts/proofread.sh` opens a chapter of the Quarto book next to the pre-Quarto
R Markdown original it was migrated from, so the two can be compared line by line.
It is the tool for the copyright/proofreading pass described in `CLAUDE.md`.

The originals live in `archive/` (248 files) and are kept for exactly this purpose.
They are not part of the rendered book and nothing in the build depends on them.

## Prerequisite

The script drives VS Code through its command-line helper, so `code` has to be on
your `PATH`. If it is not, the script fails with `code: command not found`. To
install it, open VS Code and run <kbd>Cmd</kbd>+<kbd>Shift</kbd>+<kbd>P</kbd> →
**Shell Command: Install 'code' command in PATH**.

## Usage

```bash
scripts/proofread.sh <chapter> [--diff]
```

`<chapter>` is the chapter code in lowercase (`2a`, `2e`, `10a`), not a file path.
It can be run from anywhere; the script changes to the repo root itself.

```bash
scripts/proofread.sh 2a          # two tabs: the .qmd to edit, the .Rmd to read
scripts/proofread.sh 2a --diff   # one side-by-side diff tab instead
```

It prints what it resolved before opening anything:

```
$ scripts/proofread.sh 2a
new: 02-odes/02a-modeling-gene-circuits.qmd
old: archive/02A.Rmd  (old_id=02A)
  tip: press Cmd-\ to put them side by side
```

### The two modes

- **Default** adds both files as tabs to your VS Code window. Press
  <kbd>Cmd</kbd>+<kbd>\\</kbd> to split them into side-by-side panes. Use this when
  you are editing the `.qmd` and want the original visible as a reference.
- **`--diff`** opens a single diff editor, old on the left and new on the right.
  Use this to hunt for dropped prose or dropped code, which is the main risk the
  pass is looking for.

### Where the files open

The script passes no window flag, so VS Code uses its default and opens the files
in your **last active window**. Both files are inside the repo, so with the project
already open they appear as tabs in that window and in the existing Explorer tree.

If no VS Code window is open at all, VS Code opens a new one containing just those
two files, without the project folder as the workspace, so the Explorer will not
show the rest of the book. Open the project first to avoid that (or add
`-r`/`--reuse-window` to the `code` calls in the script).

## How the old file is found

The script does not guess from the filename. It reads the `old_id:` field in the
chapter's own frontmatter and looks for `archive/<OLD_ID>.Rmd`. That is what makes
split chapters resolve correctly:

| command | new file | original |
|---------|----------|----------|
| `scripts/proofread.sh 2e` | `02-odes/02e-bifurcation.qmd` | `archive/02E.Rmd` |
| `scripts/proofread.sh 2f` | `02-odes/02f-bifurcation-curves.qmd` | `archive/02E.Rmd` |

Both map to `02E.Rmd`, because one original chapter became two.

The single-language practice copies (`*-R.qmd`, `*-py.qmd`) are skipped, so the
script always picks the real dual-language chapter.

If a chapter has no archived original, because it was newly authored rather than
migrated, the script says so and opens just the `.qmd`:

```
  (no archived .Rmd; this chapter is newly authored)
```

## What to check while comparing

From the proofreading conventions in `CLAUDE.md`:

1. **Dropped explanations.** Per-figure sentences ("the plot below shows...", "In
   this example...") and any sentence interpreting a plot must survive. Condense or
   revise them freely, but do not drop them.
2. **Dropped code.** Diff the code chunks too, not just prose. Whole pieces have
   gone missing in migration before (2B's Fortran RK4 validation).
3. **Prose placement.** Anything language-agnostic (what a plot means, a stability
   conclusion, a benchmark takeaway) belongs in the shared concept text before the
   two implementations, or in a shared paragraph after the Python block. It must not
   sit only inside the R implementation, which would leave Python-only readers with
   bare code.
4. **R/Python parity.** Same algorithm, same default parameters, not merely a
   plausible Python translation.
5. **Citations.** Any named method, model, or dataset cites its original source.
6. **Copyright.** Flag close paraphrasing or verbatim text from outside sources for
   human review rather than quietly rewriting it.
7. **Wording.** No em-dashes in chapter prose; use a comma, colon, semicolon,
   parentheses, or a sentence break.
