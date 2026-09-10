#!/usr/bin/env python3
"""Rewrap each Markdown paragraph onto a single source line.

The book's prose is authored one-paragraph-per-line (a blank line ends a
paragraph, so source line breaks inside one are cosmetic and Pandoc collapses
them to spaces). This script normalizes chapters that were still hard-wrapped,
which changes the source only, never the rendered output.

Joined:    plain paragraphs, list items with their continuation lines,
           blockquotes.
Verbatim:  YAML frontmatter, fenced code blocks, raw HTML blocks, $$..$$ math
           (single- and multi-line), \\begin{..}..\\end{..} LaTeX blocks,
           headings, ::: div/callout fences, pipe tables, thematic breaks,
           and 4-space-indented code.

Usage:  python3 scripts/unwrap-paragraphs.py <file.qmd> [more.qmd ...]
"""
import re
import sys

FENCE    = re.compile(r'^\s*(`{3,}|~{3,})')
HTMLOPEN = re.compile(r'^<([A-Za-z][\w-]*)\b')
DIV      = re.compile(r'^:{3,}')
HEAD     = re.compile(r'^#{1,6}\s')
LISTIT   = re.compile(r'^(\s*)([-*+]|\d+[.)]|\(\d+\))\s+')
QUOTE    = re.compile(r'^\s*>\s?')
TABLE    = re.compile(r'^\s*\|')
RULE     = re.compile(r'^\s*([-*_])(\s*\1){2,}\s*$')
CODEIND  = re.compile(r'^ {4,}\S')
TEXOPEN  = re.compile(r'^\s*\\begin\{([a-zA-Z*]+)\}')


def is_break(l):
    """A line that must begin its own block (never absorbed into a join)."""
    return (not l.strip() or FENCE.match(l) or DIV.match(l) or HEAD.match(l)
            or TABLE.match(l) or RULE.match(l) or HTMLOPEN.match(l)
            or CODEIND.match(l) or TEXOPEN.match(l) or '$$' in l)


def unwrap(text):
    lines = text.split('\n')
    out, i, n = [], 0, len(lines)

    if lines and lines[0].strip() == '---':          # YAML frontmatter
        out.append(lines[0]); i = 1
        while i < n and lines[i].strip() != '---':
            out.append(lines[i]); i += 1
        if i < n:
            out.append(lines[i]); i += 1

    def passthrough_until(pred, include_first=True):
        """Copy lines verbatim until pred(line) is true (that line included)."""
        nonlocal i
        if include_first:
            out.append(lines[i].rstrip()); i += 1
        while i < n and not pred(lines[i]):
            out.append(lines[i].rstrip()); i += 1
        if i < n:
            out.append(lines[i].rstrip()); i += 1

    while i < n:
        line = lines[i]

        m = FENCE.match(line)                         # ``` / ~~~ code block
        if m:
            mk = m.group(1)
            close = re.compile(r'^\s*' + re.escape(mk[0]) + '{' + str(len(mk)) + r',}\s*$')
            passthrough_until(close.match)
            continue

        m = HTMLOPEN.match(line)                      # raw HTML block
        if m:
            closer = '</%s>' % m.group(1)
            if closer in line:
                out.append(line.rstrip()); i += 1
            else:
                passthrough_until(lambda l, c=closer: l.lstrip().startswith(c))
            continue

        if '$$' in line:                              # display math
            if line.count('$$') >= 2:
                out.append(line.rstrip()); i += 1     # opens and closes here
            else:
                passthrough_until(lambda l: '$$' in l)
            continue

        m = TEXOPEN.match(line)                       # \begin{..} .. \end{..}
        if m:
            closer = r'\end{%s}' % m.group(1)
            if closer in line:
                out.append(line.rstrip()); i += 1
            else:
                passthrough_until(lambda l, c=closer: c in l)
            continue

        if (not line.strip() or DIV.match(line) or HEAD.match(line)
                or TABLE.match(line) or RULE.match(line) or CODEIND.match(line)):
            out.append(line.rstrip()); i += 1
            continue

        if QUOTE.match(line):                         # blockquote -> one line
            indent = re.match(r'^\s*', line).group(0)
            parts = []
            while i < n and QUOTE.match(lines[i]):
                parts.append(QUOTE.sub('', lines[i], count=1).strip()); i += 1
            out.append((indent + '> ' + ' '.join(p for p in parts if p)).rstrip())
            continue

        # list item, or plain paragraph -> join continuation lines
        parts = [line.rstrip()]; i += 1
        while i < n and not is_break(lines[i]) and not LISTIT.match(lines[i]) \
                and not QUOTE.match(lines[i]):
            parts.append(lines[i].strip()); i += 1
        out.append(' '.join(p if k == 0 else p.strip() for k, p in enumerate(parts)))

    return '\n'.join(out)


if __name__ == '__main__':
    for path in sys.argv[1:]:
        src = open(path, encoding='utf-8').read()
        new = unwrap(src)
        if new != src:
            open(path, 'w', encoding='utf-8').write(new)
            print('rewrapped %s' % path)
