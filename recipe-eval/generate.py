#!/usr/bin/env python3
"""Context-free code generation for recipe-box evaluation.

For each box prompt and each of K samples, call the Claude CLI in print mode
inside an EMPTY working directory with ALL tools disabled -> a fresh, one-shot,
context-free generation (no book, no files, no web, no self-execution). Extract
the single python code block and save it as gen/<box>_s<k>.py.
"""
import json, os, re, subprocess, sys, tempfile, concurrent.futures as cf

BASE = os.environ.get("RECIPE_EVAL_BASE", "/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval")
GEN = os.path.join(BASE, "gen"); os.makedirs(GEN, exist_ok=True)
MODEL = "claude-opus-4-8"
K = 5
WORKERS = 6

WRAPPER = (
 'Respond with ONLY one complete Python script inside a single triple-backtick '
 'python code block, and no prose outside it. Put `import matplotlib` then '
 '`matplotlib.use("Agg")` before importing pyplot, and save any figure with '
 'plt.savefig("{fig}") instead of plt.show(). Print every numerical result you '
 'compute, each on its own line with a clear label.\n\nUSER REQUEST:\n{prompt}'
)

def extract_code(md):
    m = re.search(r"```(?:python)?\s*\n(.*?)```", md, re.S)
    return m.group(1) if m else md

def one(box, k, prompt):
    py = os.path.join(GEN, f"{box}_s{k}.py")
    fig = os.path.join(GEN, f"{box}_s{k}.png")
    if os.path.exists(py):
        return box, k, "cached"
    wrapper = WRAPPER.format(fig=fig, prompt=prompt)
    with tempfile.TemporaryDirectory() as d:
        try:
            r = subprocess.run(
                ["claude","-p",wrapper,"--model",MODEL,"--allowedTools","",
                 "--output-format","text"],
                cwd=d, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=300)
            code = extract_code(r.stdout)
            open(py,"w").write(code)
            open(py+".raw.md","w").write(r.stdout)
            return box, k, ("ok" if code.strip() else "empty")
        except subprocess.TimeoutExpired:
            return box, k, "timeout"
        except Exception as e:
            return box, k, f"err:{e}"

def main():
    prompts = json.load(open(sys.argv[1]))
    jobs = [(b, k, p) for b, p in prompts.items() for k in range(1, K+1)]
    print(f"generating {len(jobs)} scripts ({len(prompts)} boxes x {K}) with {MODEL}", flush=True)
    status = {}
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for box, k, st in ex.map(lambda a: one(*a), jobs):
            status[f"{box}_s{k}"] = st
            print(f"  {box}_s{k}: {st}", flush=True)
    json.dump(status, open(os.path.join(BASE,"gen_status.json"),"w"), indent=1)
    print("done:", {v: sum(1 for x in status.values() if x==v) for v in set(status.values())})

if __name__ == "__main__":
    main()
