#!/usr/bin/env python3
"""Quarto pre-render hook: scan every chapter for numbered Recipe callouts and
(re)generate recipe-index-generated.qmd, a table fragment that the Method and library
reference chapter (reference-tables.qmd) pulls in with an {{< include >}} shortcode.
Each row links a recipe's number to its section and shows its Objective. The fragment
is git-ignored and rebuilt here, so the index stays in sync as recipe boxes change.
Run scripts/number-recipes.py first so the callout titles carry their numbers.
"""
import glob
import os
import re

ROOT = os.environ.get("QUARTO_PROJECT_DIR") or os.getcwd()
HEAD = re.compile(r"^## (.+)$", re.M)
CALL = re.compile(r'::: \{\.callout-note title="(Recipe[^"]*)"\}\n(.*?)\n:::', re.S)


def slug(title):
    """Reproduce Quarto's (pandoc's) section-anchor slug for a heading title.
    Pandoc keeps Unicode alphanumerics (so a Greek 'π' survives in the anchor),
    drops other punctuation except . _ -, turns spaces into hyphens, lowercases,
    and removes everything up to the first ASCII letter (the section number)."""
    s = title.lower()
    s = re.sub(r"^[^a-z]+", "", s)                              # drop the leading "8a.1 " number
    s = s.replace(" ", "-")
    s = "".join(c for c in s if c.isalnum() or c in "._-")      # keep unicode letters (e.g. π)
    s = re.sub(r"-+", "-", s).strip("-")
    return s


def main():
    rows = []
    total = 0
    for f in sorted(glob.glob(os.path.join(ROOT, "[0-9][0-9]-*", "*.qmd"))):
        if f.endswith(("-R.qmd", "-py.qmd")):
            continue
        text = open(f, encoding="utf-8").read()
        heads = [(m.start(), m.group(1).strip()) for m in HEAD.finditer(text)]
        rel = os.path.relpath(f, ROOT)
        for m in CALL.finditer(text):
            label, body = m.group(1), m.group(2)
            sec = None
            for pos, title in heads:
                if pos < m.start():
                    sec = title
                else:
                    break
            if sec is None:
                continue
            number = label.replace("Recipe ", "")          # "2B.1.1"
            parts = sec.split(None, 1)
            name = parts[1] if len(parts) > 1 else sec       # section title without its number
            rows.append(f"[{number}]({rel}#{slug(sec)}) | {name}")
            total += 1

    # Pack two recipes per line (four columns) so the index stays compact.
    lines = []
    for i in range(0, len(rows), 2):
        left = rows[i]
        right = rows[i + 1] if i + 1 < len(rows) else " | "
        lines.append(f"| {left} | {right} |")
    frag = ["| Recipe | Method / application | Recipe | Method / application |",
            "|--------|------------------------------------------|--------|------------------------------------------|",
            *lines, ""]
    dest = os.path.join(ROOT, "recipe-index-generated.qmd")
    open(dest, "w", encoding="utf-8").write("\n".join(frag) + "\n")
    print(f"make-recipe-index: indexed {total} recipe boxes")


if __name__ == "__main__":
    main()
