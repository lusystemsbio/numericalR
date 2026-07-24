#!/usr/bin/env python3
"""Quarto pre-render hook: scan every chapter for Recipe callouts and (re)generate
recipe-index.qmd, an auto-maintained index of all recipe boxes. Each entry links to
the section the box sits in and shows its Objective, so the index stays in sync as
recipe boxes are added or removed. recipe-index.qmd is git-ignored and rebuilt here.
"""
import glob
import os
import re

ROOT = os.environ.get("QUARTO_PROJECT_DIR") or os.getcwd()


def slug(title):
    """Reproduce Quarto's section-anchor slug for a heading title."""
    s = title.lower()
    s = re.sub(r"^[^a-z]+", "", s)          # drop leading non-letters (e.g. the "2")
    s = s.replace(" ", "-")
    s = re.sub(r"[^a-z0-9._-]", "", s)       # drop other punctuation (: ' , etc.)
    s = re.sub(r"-+", "-", s).strip("-")
    return s


def chapter_title(text, path):
    m = re.search(r'^title:\s*"?(.*?)"?\s*$', text, re.M)
    return m.group(1) if m else os.path.basename(path)


def main():
    files = sorted(
        f for f in glob.glob(os.path.join(ROOT, "[0-9][0-9]-*", "*.qmd"))
        if not f.endswith(("-R.qmd", "-py.qmd"))
    )

    chapters = []   # (chapter_title, relpath, [(sec_title, anchor, objective), ...])
    total = 0
    heading_re = re.compile(r"^## (.+)$", re.M)
    for f in files:
        text = open(f, encoding="utf-8").read()
        # positions of every level-2 heading
        heads = [(m.start(), m.group(1).strip()) for m in heading_re.finditer(text)]
        entries = []
        for cm in re.finditer(
            r'::: \{\.callout-note title="Recipe"\}\n(.*?)\n:::', text, re.S
        ):
            # enclosing section = last ## heading before this callout
            sec = None
            for pos, title in heads:
                if pos < cm.start():
                    sec = title
                else:
                    break
            if sec is None:
                continue
            om = re.search(r"\*\*Objective\.\*\*\s*(.*?)(?:\n\n|\Z)", cm.group(1), re.S)
            obj = " ".join(om.group(1).split()) if om else ""
            obj = obj.split(". ")[0].rstrip(".") + "." if obj else ""
            entries.append((sec, slug(sec), obj))
            total += 1
        if entries:
            rel = os.path.relpath(f, ROOT)
            chapters.append((chapter_title(text, f), rel, entries))

    out = ['---', 'title: "Recipe index"', '---', '',
           "Every recipe box in the book, grouped by chapter and linked to its section. "
           "Each box states an objective, a method, a test case, and how to verify the "
           "result, and can be turned into a prompt with its **Copy prompt** button "
           "(see [Part 1A](01-intro/01a-getting-started.qmd)).", '']
    for title, rel, entries in chapters:
        out.append(f"**{title}**")
        out.append("")
        for sec, anchor, obj in entries:
            link = f"[{sec}]({rel}#{anchor})"
            out.append(f"- {link}" + (f" --- {obj}" if obj else ""))
        out.append("")

    dest = os.path.join(ROOT, "recipe-index.qmd")
    open(dest, "w", encoding="utf-8").write("\n".join(out) + "\n")
    print(f"make-recipe-index: indexed {total} recipe boxes across {len(chapters)} chapters")


if __name__ == "__main__":
    main()
