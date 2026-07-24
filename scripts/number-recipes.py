#!/usr/bin/env python3
"""Number the Recipe callouts in the chapter sources: each becomes
title="Recipe <section>.<k>", where <section> is the number of the enclosing
level-2 section (e.g. 2B.1) and <k> counts recipes within that section (usually 1).
Idempotent: re-run after adding boxes or renumbering sections. This edits the .qmd
sources in place, so run it manually and commit the result (not a render hook)."""
import glob
import os
import re

ROOT = os.environ.get("QUARTO_PROJECT_DIR") or os.getcwd()
HEAD = re.compile(r"^## (\S+)", re.M)
CALL = re.compile(r'::: \{\.callout-note title="Recipe[^"]*"\}')


def main():
    total = 0
    for f in sorted(glob.glob(os.path.join(ROOT, "[0-9][0-9]-*", "*.qmd"))):
        if f.endswith(("-R.qmd", "-py.qmd")):
            continue
        text = open(f, encoding="utf-8").read()
        heads = [(m.start(), m.group(1)) for m in HEAD.finditer(text)]
        out, last, counter = [], 0, {}
        for m in CALL.finditer(text):
            sec = None
            for pos, tok in heads:
                if pos < m.start():
                    sec = tok
                else:
                    break
            if sec is None:
                continue
            counter[sec] = counter.get(sec, 0) + 1
            out.append(text[last:m.start()])
            out.append(f'::: {{.callout-note title="Recipe {sec}.{counter[sec]}"}}')
            last = m.end()
            total += 1
        if out:
            out.append(text[last:])
            open(f, "w", encoding="utf-8").write("".join(out))
    print(f"number-recipes: numbered {total} recipe boxes")


if __name__ == "__main__":
    main()
