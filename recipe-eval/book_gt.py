#!/usr/bin/env python3
"""Capture the book's reference plotted data per box (ground truth).

For each chapter that has recipe boxes, extract its Python chunks in order,
tag each with its enclosing '## <sec>' heading, and run them (in the book's own
interpreter) under the capture shim, so every plotted series is recorded with
its section. Result: book_gt.json = {section_code: [plotted series...]}.
Run from the repo root (book code reads data files by relative path)."""
import glob, json, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CAPTURE = os.path.join(HERE, "capture")
PY_BOOK = os.environ.get("PY_BOOK", "/usr/local/bin/python3")   # book's own env
HEAD = re.compile(r'^## (\d+[A-Z]\.\d+)', re.M)
PYBLK = re.compile(r'```\{python[^}]*\}\n(.*?)```', re.S)

def chapter_blocks(path):
    """[(section_code, python_code)] in document order."""
    t = open(path).read()
    heads = [(m.start(), m.group(1)) for m in HEAD.finditer(t)]
    def sec_at(pos):
        s = None
        for hp, code in heads:
            if hp < pos: s = code
            else: break
        return s
    return [(sec_at(m.start()), m.group(1)) for m in PYBLK.finditer(t)]

def run_chapter(path, capout):
    blocks = [(s, c) for s, c in chapter_blocks(path) if s]
    if not blocks: return {}
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(blocks, f); bf = f.name
    env = dict(os.environ, PYTHONPATH=CAPTURE, CAPTURE_OUT=capout, MPLBACKEND="Agg")
    try:
        subprocess.run([PY_BOOK, os.path.join(HERE, "gt_driver.py"), bf],
                       cwd=REPO, env=env, capture_output=True, text=True, timeout=400)
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT {os.path.basename(path)}")
    series = json.load(open(capout)) if os.path.exists(capout) else []
    by_sec = {}
    for s in series:
        by_sec.setdefault(s.get("sec"), []).append(s)
    return by_sec

def main():
    files = sorted(glob.glob(os.path.join(REPO, "[0-9][0-9]-*", "*.qmd")))
    files = [f for f in files if not f.endswith(("-R.qmd", "-py.qmd"))]
    only = sys.argv[1] if len(sys.argv) > 1 else None   # optional folder filter e.g. 02-
    gt = {}
    tmp = tempfile.mkdtemp()
    for f in files:
        if only and only not in f: continue
        if not HEAD.search(open(f).read()): continue
        capout = os.path.join(tmp, os.path.basename(f) + ".json")
        by_sec = run_chapter(f, capout)
        for sec, series in by_sec.items():
            gt.setdefault(sec, []).extend(series)
        n = sum(len(v) for v in by_sec.values())
        print(f"  {os.path.basename(f)}: {n} series across {len(by_sec)} sections", flush=True)
    json.dump(gt, open(os.path.join(HERE, "book_gt.json"), "w"))
    print(f"book ground truth: {len(gt)} sections captured -> book_gt.json")

if __name__ == "__main__":
    main()
