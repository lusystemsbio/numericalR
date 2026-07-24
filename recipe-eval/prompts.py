#!/usr/bin/env python3
"""Turn a recipe box into the exact context-free prompt its Copy-prompt button
produces (mirrors recipe-copy.html), for a chosen language."""
import re, glob, json, sys

REPO = "/Users/lvmy/neu/teaching/numericalR/numericalR"
CALL = re.compile(r'::: \{\.callout-note title="Recipe ([^"]*)"\}\n(.*?)\n:::', re.S)

def _fields(body):
    f = {}
    for name in ["Objective","Model","Method","Test","Show","Verification"]:
        m = re.search(r'\*\*'+name+r'\.\*\*\s*(.*?)(?=\n\n\*\*|\Z)', body, re.S)
        if m:
            f[name] = " ".join(m.group(1).split()).replace("`","")  # copy button drops backticks
    return f

def all_boxes():
    boxes = {}
    for fn in sorted(glob.glob(REPO + "/[0-9][0-9]-*/*.qmd")):
        if fn.endswith(("-R.qmd","-py.qmd")): continue
        t = open(fn).read()
        for m in CALL.finditer(t):
            boxes[m.group(1)] = _fields(m.group(2))
    return boxes

def prompt_for(fields, lang="Python"):
    frag = lambda s: (s or "").strip().rstrip(".").strip()
    lc = lambda s: (s[0].lower()+s[1:]) if s else s
    p  = f"In {lang}, " + lc(frag(fields["Objective"])) + ". "
    if fields.get("Model"):  p += "The model is " + lc(frag(fields["Model"])) + ". "
    p += "Use " + lc(frag(fields["Method"])) + ". "
    p += "Test it on " + frag(fields["Test"]) + ". "
    if fields.get("Show"):   p += "Produce " + lc(frag(fields["Show"])) + ". "
    p += "As a separate check, confirm that " + lc(frag(fields["Verification"])) + ". "
    p += ("Implement the method explicitly with short comments rather than calling a "
          "routine that does it in one step, and explain in one sentence why that "
          "check confirms the result.")
    return p

if __name__ == "__main__":
    boxes = all_boxes()
    pilot = ["2C.1.1","2E.3.1","5A.4.1","6A.4.1","3A.1.1",
             "8D.2.1","9B.3.1","10C.5.1","7E.3.1","10B.1.1"]
    out = {n: prompt_for(boxes[n], "Python") for n in pilot}
    json.dump(out, open(sys.argv[1] if len(sys.argv)>1 else "/dev/stdout","w"), indent=1)
    for n in pilot:
        print(f"\n===== {n} =====\n{out[n]}")
