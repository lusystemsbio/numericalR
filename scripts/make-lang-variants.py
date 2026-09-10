#!/usr/bin/env python3
"""Generate single-language practice copies of every dual-language chapter,
written IN PLACE next to each original .qmd (e.g. 02-odes/02b-numerical-integration-R.qmd
beside 02-odes/02b-numerical-integration.qmd). Run this once as a setup step after
installing/cloning the book; the files are for readers to practice with and are NOT part
of the rendered book (they are not listed in _quarto.yml and are git-ignored).

Because each copy sits in its chapter's own folder, every relative resource the chapter
uses (images/<id>/..., src/*.f90, data/...) resolves exactly as for the original, and
rendering a copy inside the project also picks up the book's bibliography and filters.

R-only  (<name>-R.qmd):  Python implementation heading + ```{python} chunks removed,
                         engine forced to knitr.
Py-only (<name>-py.qmd): R implementation heading + ```{r} chunks removed,
                         engine forced to jupyter (needs no R).

Prose is kept verbatim (concept text and the shared post-Python discussion), and the
retained language's ### implementation heading is kept, matching the design settled
with the author. Display-only listings are treated by language just like executable
chunks: ```r / ```python (no braces) are dropped from the other language's variant,
while shared ones (```bash / ```sh / ```fortran, or a plain fence) are kept in both. Eligibility is decided purely by content: any page carrying both
a `### R implementation` and a `### Python implementation` heading gets variants, including
exercises pages (1D, 2G); pages with only one language, or none, get none. Idempotent and fast.
"""
import os
import re
import glob

FENCE = re.compile(r"^(`{3,})(.*)$")
# an R/Python code fence, executable OR display-only: ```{r} ```{python} ```r ```python
# (```bash / ```sh / ```fortran are shared infrastructure and match nothing here -> kept).
LANG = re.compile(r"^\{?\.?\s*(r|python)\b")
IMPL = re.compile(r"^###\s+(R|Python)\s+implementation\s*\{\.impl\}\s*$", re.I)


def strip(text, keep):
    """Return `text` with the non-`keep` language's impl heading + code removed.
    keep is 'r' or 'python'."""
    drop = "python" if keep == "r" else "r"
    lines = text.split("\n")

    # split frontmatter (leading --- ... ---) from body
    fm, body = [], lines
    if lines and lines[0].strip() == "---":
        for i in range(1, len(lines)):
            if lines[i].strip() == "---":
                fm, body = lines[: i + 1], lines[i + 1:]
                break

    out = []
    in_fence = False
    fence_len = 0
    dropping = False  # inside a to-be-removed executable chunk
    for ln in body:
        m = FENCE.match(ln)
        if m and not in_fence:
            in_fence = True
            fence_len = len(m.group(1))
            info = m.group(2).strip()
            lm = LANG.match(info)
            dropping = bool(lm and lm.group(1) == drop)
            if not dropping:
                out.append(ln)
            continue
        if in_fence:
            # a bare fence of >= opening length closes the block
            if re.match(r"^`{%d,}\s*$" % fence_len, ln):
                if not dropping:
                    out.append(ln)
                in_fence = False
                dropping = False
                continue
            if not dropping:
                out.append(ln)
            continue
        # outside any fence: drop the other language's impl heading line
        im = IMPL.match(ln)
        if im and im.group(1).lower().startswith(drop[0]) and \
           (im.group(1).lower() == "python") == (drop == "python"):
            continue
        out.append(ln)

    # collapse 3+ blank lines (left by removals) down to one
    joined = "\n".join(out)
    joined = re.sub(r"\n{3,}", "\n\n", joined).strip("\n")

    fm = set_engine(fm, keep)
    note = ("<!-- Auto-generated %s-only version of this chapter; the %s implementation "
            "has been removed. Regenerate via scripts/make-lang-variants.py. -->"
            % (keep.upper() if keep == "r" else "Python", "Python" if keep == "r" else "R"))
    return "\n".join(fm) + "\n\n" + note + "\n\n" + joined + "\n"


def set_engine(fm, keep):
    """Insert engine directives into the frontmatter line list (before closing ---)."""
    adds = ["engine: knitr"] if keep == "r" else ["engine: jupyter", "jupyter: python3"]
    if len(fm) >= 2 and fm[0].strip() == "---" and fm[-1].strip() == "---":
        keys = {a.split(":")[0].strip() for a in adds}
        kept = [l for l in fm[1:-1] if l.split(":")[0].strip() not in keys]
        return ["---"] + kept + adds + ["---"]
    return ["---"] + adds + ["---"]  # no frontmatter -> create one


def main():
    proj = os.environ.get("QUARTO_PROJECT_DIR") or os.getcwd()

    n = 0
    for src in sorted(glob.glob(os.path.join(proj, "[0-9][0-9]-*", "*.qmd"))):
        base = os.path.basename(src)
        if base.endswith("-R.qmd") or base.endswith("-py.qmd"):
            continue  # never regenerate from a previously generated variant
        text = open(src, encoding="utf-8").read()
        if not (re.search(r"^###\s+R\s+implementation", text, re.M) and
                re.search(r"^###\s+Python\s+implementation", text, re.M)):
            continue  # no paired R/Python implementations -> no variants
        stem = os.path.splitext(src)[0]  # write in place, beside the original
        open(stem + "-R.qmd", "w", encoding="utf-8").write(strip(text, "r"))
        open(stem + "-py.qmd", "w", encoding="utf-8").write(strip(text, "python"))
        n += 1
    print("make-lang-variants: wrote in-place R/Python practice copies for %d chapters" % n)


if __name__ == "__main__":
    main()
