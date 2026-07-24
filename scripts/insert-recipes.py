#!/usr/bin/env python3
"""Insert recipe boxes: for each (file, section-code -> body), place a
callout right before the first '### R implementation {.impl}' after the
'## <code> ...' heading. Title is 'Recipe <code>'; run number-recipes.py
afterward to normalize to 'Recipe <code>.<k>'. Idempotent-ish: refuses if a
Recipe box already exists in that section."""
import re, sys, json
def insert(path, boxes):
    lines = open(path).read().split("\n")
    # work from bottom to top so indices stay valid
    items = []
    for code, body in boxes.items():
        # find heading '## <code> '
        h = next((i for i,l in enumerate(lines) if re.match(r'## '+re.escape(code)+r'(\s|$)', l)), None)
        if h is None: print(f"  {path}: heading {code} NOT FOUND"); continue
        # next section heading (to bound the search)
        nxt = next((i for i in range(h+1,len(lines)) if re.match(r'## \S', lines[i])), len(lines))
        # bail if a recipe already there
        if any('callout-note title="Recipe' in lines[i] for i in range(h,nxt)):
            print(f"  {path}: {code} already has a Recipe box, skipping"); continue
        impl = next((i for i in range(h+1,nxt) if lines[i].strip().startswith('### R implementation')), None)
        if impl is None: print(f"  {path}: {code} has no R implementation block, skipping"); continue
        items.append((impl, code, body))
    for impl, code, body in sorted(items, reverse=True):
        block = [f'::: {{.callout-note title="Recipe {code}"}}', body.strip("\n"), ':::', '']
        lines[impl:impl] = block
        print(f"  {path}: inserted Recipe {code}")
    open(path,"w").write("\n".join(lines))

if __name__ == "__main__":
    data = json.load(sys.stdin)
    for path, boxes in data.items():
        insert(path, boxes)
